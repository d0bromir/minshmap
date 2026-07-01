"""minSHmap - the whole sketch-based read mapper in ~150 lines, nothing else.

The single, pedagogical version: one sketcher (ntHash), no max_matches, no
second-best/mapq - we just keep the single best hit. `minshmap.cpp` is a
line-for-line port of this file (same algorithm, same output) that additionally
offers `-j` for output-identical parallelism; keep the two in lockstep.

  1. SKETCH  - FracMinHash: keep the small fraction `hfrac` of k-mers, chosen by
               a rolling, strand-canonical ntHash (a k-mer and its rev-comp get
               the same hash, so either strand is found with one index).
  2. INDEX   - reference: hash -> [(segment, position, strand), ...].
  3. MAP     - sketch the read; rank its k-mers rarest-first; scatter the rarest
               ones into overlapping reference blocks; for each
               block compute containment while pruning with the SEED HEURISTIC

                   sh = 1 - (seeds_used - matches) / m   (upper bound it can reach)

               the moment sh < theta the block is provably hopeless -> skip it.

Usage:  python minshmap.py reference.fa reads.fa -k 15 -r 0.05 -t 0.9
"""
from __future__ import annotations
import argparse

MASK64 = (1 << 64) - 1
_LUT_FW = {"A": 0x3C8BFBB395C60474, "C": 0x3193C18562A02B4C,
           "G": 0x20323ED082572324, "T": 0x295549F54BE24456}
_LUT_RC = {"A": _LUT_FW["T"], "C": _LUT_FW["G"], "G": _LUT_FW["C"], "T": _LUT_FW["A"]}


def _rotl(x, r): r &= 63; return ((x << r) | (x >> (64 - r))) & MASK64
def _rotr(x, r): r &= 63; return ((x >> r) | (x << (64 - r))) & MASK64


def sketch(seq, k, hfrac):
    """FracMinHash via a rolling ntHash. Yields (pos, hash, strand).      O(len)"""
    n = len(seq)
    if n < k:
        return
    thr = int(hfrac * (MASK64 + 1))
    fw = [_LUT_FW.get(c, 0) for c in seq]
    rc = [_LUT_RC.get(c, 0) for c in seq]
    h_fw = h_rc = 0
    for t in range(k):                       # hash of the first window
        h_fw ^= _rotl(fw[t], k - 1 - t)
        h_rc ^= _rotl(rc[t], t)
    for i in range(n - k + 1):
        h = h_fw if h_fw <= h_rc else h_rc   # canonical: min over the two strands
        if h <= thr:
            yield i, h, 0 if h_fw <= h_rc else 1
        if i + k < n:                        # roll the window one base
            h_fw = _rotl(h_fw, 1) ^ _rotl(fw[i], k) ^ fw[i + k]
            h_rc = _rotr(h_rc, 1) ^ _rotr(rc[i], 1) ^ _rotl(rc[i + k], k - 1)


def build_index(refs, k, hfrac):
    """Reference sketch inverted: hash -> [(segm_id, pos, strand), ...]."""
    index, names, lengths = {}, [], []
    for sid, (name, seq) in enumerate(refs):
        names.append(name); lengths.append(len(seq))
        for pos, h, strand in sketch(seq, k, hfrac):
            index.setdefault(h, []).append((sid, pos, strand))
    return index, names, lengths


def map_read(seq, index, k, hfrac, theta):
    """Best reference block for the read: (segm, t_start, t_end, score, codir) or None."""
    sk = list(sketch(seq, k, hfrac))
    m = len(sk)                              # informative k-mers in the read
    if m == 0:
        return None
    B = max(len(seq), 1)                     # each candidate block is read-length-wide

    # One seed per distinct read k-mer: remember its reference hits and read strand.
    seeds = {}
    for _, h, strand in sk:
        if h not in seeds:
            seeds[h] = (index.get(h, ()), strand)
    order = sorted(seeds, key=lambda h: len(seeds[h][0]))   # rarest k-mers first

    # Candidate blocks (pre-#3 reference): sorted(cand) ascending, prune vs theta.
    S = int((1 - theta) * m) + 1
    cand = set()
    for h in order[:S]:
        for sid, pos, _ in seeds[h][0]:
            b = pos // B
            cand.add((sid, b))
            if b:
                cand.add((sid, b - 1))

    best = None
    for sid, b in sorted(cand):
        lo, hi = b * B, (b + 2) * B
        used = matches = codir = 0
        r_min = r_max = -1
        for h in order:
            used += 1
            hits, rstrand = seeds[h]
            for s2, pos, hstrand in hits:
                if s2 == sid and lo <= pos < hi:
                    matches += 1
                    codir += 1 if hstrand == rstrand else -1
                    r_min = pos if r_min < 0 else min(r_min, pos)
                    r_max = max(r_max, pos)
                    break
            if 1 - (used - matches) / m < theta:
                break
        else:
            score = matches / m
            if score >= theta and (best is None or score > best[3]):
                best = (sid, r_min, r_max + k, score, codir)
    return best


def read_fasta(path):
    """Yield (name, sequence) per FASTA record."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def main():
    p = argparse.ArgumentParser(description="minSHmap - minimal sketch mapper")
    p.add_argument("reference"); p.add_argument("reads")
    p.add_argument("-k", type=int, default=15)
    p.add_argument("-r", "--hfrac", type=float, default=0.05)
    p.add_argument("-t", "--theta", type=float, default=0.9)
    a = p.parse_args()
    index, names, lengths = build_index(read_fasta(a.reference), a.k, a.hfrac)
    for name, seq in read_fasta(a.reads):
        mp = map_read(seq, index, a.k, a.hfrac, a.theta)
        if mp is None:
            continue
        sid, ts, te, score, codir = mp
        nmatch = round(score * len(seq))
        print("\t".join(str(x) for x in (
            name, len(seq), 0, len(seq), "+" if codir >= 0 else "-",
            names[sid], lengths[sid], ts, te, nmatch, te - ts, 60)))


if __name__ == "__main__":
    main()
