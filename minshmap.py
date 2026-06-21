"""minSHmap - a minimalistic sketch-based read mapper.

A short, readable re-implementation of `shmap` in the spirit of `minSH`:
the goal is to be as simple to *understand* as possible, not as fast as possible.

The whole idea in three steps:

  1. SKETCH   - turn each sequence into a small set of (position, hash, strand)
                k-mers via FracMinHash (keep only k-mers whose hash is small).
  2. INDEX    - for the reference, remember where every sketch k-mer occurs:
                hash -> [(segment, position, strand), ...].
  3. MAP      - for each read, sketch it, then look for the reference *window*
                that shares the most read k-mers. To avoid scanning every
                window we use the SEED HEURISTIC (the "SH" in shmap):

                    sh = 1 - (seeds_used - matches) / m

                `sh` is an upper bound on the containment a window can still
                reach while we add read k-mers (seeds) one by one, rarest
                first. If `sh` ever drops below the homology threshold, the
                window can never be good enough and we prune it immediately.

This is the same seed heuristic that `minSH` uses to guide A* alignment, here
re-used to prune candidate *locations* during read mapping.

Usage:
    python minshmap.py reference.fa reads.fa      # map reads, print PAF
    python minshmap.py --demo                     # generate data and self-test
"""

from __future__ import annotations

import argparse
import hashlib
import math
import random
from collections import namedtuple

# --------------------------------------------------------------------------- #
# Sketching: FracMinHash
# --------------------------------------------------------------------------- #
#
# A sketcher turns a sequence into a small, strand-canonical set of k-mers:
#
#     sketch(seq, k, hfrac) -> [(pos, h, strand), ...]
#
# It keeps a k-mer only if its 64-bit hash is among the smallest fraction
# `hfrac` of all hashes (FracMinHash). The hash is *canonical*: a k-mer and its
# reverse complement get the same value, so a read still matches when it comes
# from the opposite strand. `strand` is 0 (forward) or 1 (reverse complement).
#
# Three interchangeable implementations are provided to study the speed/clarity
# trade-off (all selected with `--hash`):
#
#   naive   - hash every k-mer string from scratch.            O(len * k)
#   poly    - Rabin-Karp polynomial rolling hash.              O(len)
#   nthash  - rotation-based rolling hash, mirrors shmap.      O(len)

MASK64 = (1 << 64) - 1
_COMP = str.maketrans("ACGT", "TGCA")
_VAL = {"A": 0, "C": 1, "G": 2, "T": 3}


def revcomp(s: str) -> str:
    """Reverse complement of a DNA string."""
    return s.translate(_COMP)[::-1]


def _threshold(hfrac: float) -> int:
    """Largest hash value still kept by FracMinHash."""
    return int(hfrac * (MASK64 + 1))


# ---- 1. Naive: hash each k-mer string independently ----------------------- #

def _h64(kmer: str) -> int:
    """Deterministic 64-bit hash of a k-mer (process-independent)."""
    return int.from_bytes(hashlib.blake2b(kmer.encode(), digest_size=8).digest(), "little")


def sketch_naive(seq: str, k: int, hfrac: float):
    """Reference FracMinHash: re-hash every k-mer from scratch.        O(len*k)"""
    thr = _threshold(hfrac)
    sketch = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        h_fw = _h64(kmer)
        h_rc = _h64(revcomp(kmer))
        h = min(h_fw, h_rc)  # canonical: same hash for a k-mer and its rev-comp
        if h <= thr:
            sketch.append((i, h, 0 if h_fw <= h_rc else 1))
    return sketch


# ---- 2. Polynomial (Rabin-Karp) rolling hash ------------------------------ #
#
# A k-mer is read as a base-`_B` number. Sliding the window one step costs O(1):
# drop the outgoing base's contribution, shift, add the incoming base. The
# reverse-complement hash is rolled in parallel (the modular inverse `_INVB`
# lets us drop its least-significant base without recomputing).

_M = (1 << 61) - 1        # a Mersenne prime keeps the hashes well spread
_B = 0x9E3779B1           # base / multiplier
_INVB = pow(_B, _M - 2, _M)


def sketch_poly(seq: str, k: int, hfrac: float):
    """FracMinHash with a Rabin-Karp rolling hash.                       O(len)"""
    n = len(seq)
    if n < k:
        return []
    thr = int(hfrac * _M)
    v = [_VAL.get(c, 0) for c in seq]
    bk1 = pow(_B, k - 1, _M)        # weight of the most/least significant base

    h_fw = 0
    h_rc = 0
    for t in range(k):              # hash of the first window
        h_fw = (h_fw * _B + v[t]) % _M
        h_rc = (h_rc + (3 - v[t]) * pow(_B, t, _M)) % _M

    sketch = []
    for i in range(n - k + 1):
        h = min(h_fw, h_rc)
        if h <= thr:
            sketch.append((i, h, 0 if h_fw <= h_rc else 1))
        if i + k < n:               # roll the window from i to i+1
            out, inc = v[i], v[i + k]
            h_fw = ((h_fw - out * bk1) * _B + inc) % _M
            h_rc = ((h_rc - (3 - out)) * _INVB + (3 - inc) * bk1) % _M
    return sketch


# ---- 3. ntHash-style rotation rolling hash (as in shmap) ------------------ #
#
# Each base maps to a random 64-bit word; the k-mer hash is the XOR of those
# words cyclically rotated by their offset. Rotations are cheap and invertible,
# so the window rolls in O(1) without any modular inverse. This mirrors the
# hashing in shmap's FracMinHash (LUT_fw / LUT_rc with rotl/rotr).

_LUT_FW = {
    "A": 0x3C8BFBB395C60474, "C": 0x3193C18562A02B4C,
    "G": 0x20323ED082572324, "T": 0x295549F54BE24456,
}
_LUT_RC = {"A": _LUT_FW["T"], "C": _LUT_FW["G"], "G": _LUT_FW["C"], "T": _LUT_FW["A"]}


def _rotl(x: int, r: int) -> int:
    r &= 63
    return ((x << r) | (x >> (64 - r))) & MASK64


def _rotr(x: int, r: int) -> int:
    r &= 63
    return ((x >> r) | (x << (64 - r))) & MASK64


def sketch_nthash(seq: str, k: int, hfrac: float):
    """FracMinHash with a rotation-based (ntHash-like) rolling hash.     O(len)"""
    n = len(seq)
    if n < k:
        return []
    thr = _threshold(hfrac)
    fw = [_LUT_FW.get(c, 0) for c in seq]
    rc = [_LUT_RC.get(c, 0) for c in seq]

    h_fw = h_rc = 0
    for t in range(k):              # hash of the first window
        h_fw ^= _rotl(fw[t], k - 1 - t)
        h_rc ^= _rotl(rc[t], t)

    sketch = []
    for i in range(n - k + 1):
        h = h_fw if h_fw <= h_rc else h_rc   # canonical: min of the two strands
        if h <= thr:
            sketch.append((i, h, 0 if h_fw <= h_rc else 1))
        if i + k < n:               # roll the window from i to i+1
            h_fw = _rotl(h_fw, 1) ^ _rotl(fw[i], k) ^ fw[i + k]
            h_rc = _rotr(h_rc, 1) ^ _rotr(rc[i], 1) ^ _rotl(rc[i + k], k - 1)
    return sketch


SKETCHERS = {"naive": sketch_naive, "poly": sketch_poly, "nthash": sketch_nthash}


# --------------------------------------------------------------------------- #
# Indexing the reference
# --------------------------------------------------------------------------- #

Segment = namedtuple("Segment", "name length")


def build_index(refs, k: int, hfrac: float, sketch=sketch_naive):
    """Index reference segments: hash -> [(segm_id, pos, strand), ...].

    `refs` is an iterable of (name, sequence).                          O(total)
    """
    index = {}
    segments = []
    for sid, (name, seq) in enumerate(refs):
        segments.append(Segment(name, len(seq)))
        for pos, h, strand in sketch(seq, k, hfrac):
            index.setdefault(h, []).append((sid, pos, strand))
    return index, segments


# --------------------------------------------------------------------------- #
# Mapping a read
# --------------------------------------------------------------------------- #

Mapping = namedtuple("Mapping", "segm_id t_start t_end score codir")


def _containment_with_sh(seeds, rstrand, m, sid, lo, hi, thr):
    """Containment of the read in reference window [lo, hi) of segment `sid`.

    Adds read k-mers (seeds) rarest first and tracks the seed heuristic

        sh = 1 - (seeds_used - matches) / m

    which upper-bounds the containment still reachable. Returns None as soon as
    `sh < thr` (the window is pruned), otherwise the exact containment.
    """
    seeds_used = matches = 0
    codir = r_min = r_max = 0
    r_min = math.inf
    r_max = -1
    for n_hits, occ_in_p, h, hits in seeds:
        seeds_used += occ_in_p
        matched = 0
        for s2, pos, hstrand in hits:
            if s2 == sid and lo <= pos < hi:
                matched += 1
                codir += 1 if hstrand == rstrand[h] else -1
                r_min = min(r_min, pos)
                r_max = max(r_max, pos)
        matches += min(matched, occ_in_p)  # a seed can match its window at most occ_in_p times
        if 1 - (seeds_used - matches) / m < thr:  # seed heuristic: prune hopeless window
            return None
    return matches / m, codir, r_min, r_max


def map_read(seq, index, segments, k, hfrac, theta, min_diff, sketch=sketch_naive,
             max_matches=0):
    """Find the best (and second best) reference location for one read.

    `max_matches` (0 = off) drops k-mers that occur more than that many times in
    the reference: such repetitive k-mers are uninformative and, on real genomes
    with satellite arrays, otherwise dominate the work (see shmap's
    erase_frequent_kmers / max_matches). Returns (best, second) as Mapping or None.
    """
    sk = sketch(seq, k, hfrac)
    m = len(sk)
    if m == 0:
        return None, None
    halflen = max(len(seq), 1)  # buckets are read-length windows over the reference

    # Group read k-mers by hash: where they occur in the read and their strand.
    occ, rstrand = {}, {}
    for pos, h, strand in sk:
        occ.setdefault(h, []).append(pos)
        rstrand.setdefault(h, strand)

    # One "seed" per distinct read k-mer, annotated with its rarity in the index.
    seeds = []
    for h, positions in occ.items():
        hits = index.get(h, ())
        seeds.append((len(hits), len(positions), h, hits))

    # Drop over-frequent (uninformative) k-mers; m becomes the informative count.
    if max_matches:
        seeds = [s for s in seeds if s[0] <= max_matches]
        m = sum(occ_in_p for _, occ_in_p, _, _ in seeds)
        if m == 0:
            return None, None

    seeds.sort(key=lambda s: s[0])  # rarest first: strongest, cheapest evidence


    # ---- Seed candidate windows using the rarest seeds ----------------------
    # Any read-length region of the reference lies fully inside some window
    # because consecutive windows overlap (we add each hit to bucket b and b-1).
    theta2 = theta - min_diff
    S = int((1 - theta2) * m) + 1  # enough seeds that any homology shows >=1 match
    candidates = set()
    used = 0
    for n_hits, occ_in_p, h, hits in seeds:
        if used >= S:
            break
        if n_hits == 0:
            continue
        for sid, pos, _ in hits:
            b = pos // halflen
            candidates.add((sid, b))
            if b > 0:
                candidates.add((sid, b - 1))
        used += 1

    # ---- Refine each candidate window, pruning with the seed heuristic ------
    found = []
    for sid, b in candidates:
        lo, hi = b * halflen, (b + 2) * halflen
        res = _containment_with_sh(seeds, rstrand, m, sid, lo, hi, theta)
        if res is None:
            continue
        score, codir, r_min, r_max = res
        if score >= theta:
            found.append(Mapping(sid, r_min, r_max + k, score, codir))

    if not found:
        return None, None
    found.sort(key=lambda mp: mp.score, reverse=True)
    best = found[0]
    # Second best from a clearly different region (for mapping quality).
    second = next((mp for mp in found[1:] if mp.segm_id != best.segm_id
                   or mp.t_start >= best.t_end or mp.t_end <= best.t_start), None)
    return best, second


def mapq(best, second, theta, min_diff):
    """Mapping quality (60 high / 0 ambiguous), minimap2-style."""
    score2 = second.score if second else (theta - min_diff)
    frac = 1 - score2 / best.score
    return 60 if frac > min_diff else 0


def to_paf(query, qlen, best, second, segments, theta, min_diff):
    """Format one mapping as a PAF line."""
    seg = segments[best.segm_id]
    strand = "+" if best.codir >= 0 else "-"
    nmatch = round(best.score * qlen)
    return "\t".join(str(x) for x in (
        query, qlen, 0, qlen, strand,
        seg.name, seg.length, best.t_start, best.t_end,
        nmatch, best.t_end - best.t_start, mapq(best, second, theta, min_diff),
    ))


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """Yield (name, sequence) for every record in a FASTA file."""
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


def map_reads(ref_path, reads_path, k, hfrac, theta, min_diff, sketch=sketch_naive,
              out=print, max_matches=0):
    """Index the reference and print a PAF line per mapped read."""
    index, segments = build_index(read_fasta(ref_path), k, hfrac, sketch)
    for name, seq in read_fasta(reads_path):
        best, second = map_read(seq, index, segments, k, hfrac, theta, min_diff, sketch,
                                max_matches)
        if best is not None:
            out(to_paf(name, len(seq), best, second, segments, theta, min_diff))


# --------------------------------------------------------------------------- #
# Demo / self-test
# --------------------------------------------------------------------------- #

_NUC = "ACGT"


def _random_seq(n, rng):
    return "".join(rng.choice(_NUC) for _ in range(n))


def _mutate(s, e, rng):
    """Apply substitutions/insertions/deletions at rate `e` (like minSH)."""
    out = []
    for c in s:
        if rng.random() > e:
            out.append(c)
        elif rng.random() < 1 / 3:
            out.append(rng.choice(_NUC))            # substitution
        elif rng.random() < 1 / 2:
            out.append(rng.choice(_NUC) + c)        # insertion
        # else: deletion
    return "".join(out)


def demo(k=11, hfrac=0.2, theta=0.5, min_diff=0.02, error=0.04, seed=1, sketch=sketch_nthash):
    """Generate a reference + reads with known truth and check mapping.

    Note on parameters: a read k-mer survives a mutated read only if all k of
    its bases are error-free, so the expected containment is about (1-e)^k.
    The threshold `theta` must sit below that, hence smaller k / higher theta
    trade off for noisy reads (here e=4%: (1-0.04)^11 ~ 0.64 > theta=0.5).
    """
    rng = random.Random(seed)
    ref = _random_seq(5000, rng)
    segments = [("chr_demo", len(ref))]
    index, segments = build_index([("chr_demo", ref)], k, hfrac, sketch)

    correct = 0
    total = 20
    for i in range(total):
        start = rng.randint(0, len(ref) - 400)
        true = ref[start : start + 300]
        read = _mutate(true, error, rng)
        if rng.random() < 0.5:
            read = revcomp(read)
        best, second = map_read(read, index, segments, k, hfrac, theta, min_diff, sketch)
        ok = best is not None and best.t_start <= start + 150 <= best.t_end
        correct += ok
        loc = f"{best.t_start}-{best.t_end} score={best.score:.2f}" if best else "UNMAPPED"
        print(f"read{i:02d} true~{start:5d}  ->  {loc}  {'OK' if ok else 'MISS'}")
    print(f"\nCorrectly placed {correct}/{total} reads.")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main():
    p = argparse.ArgumentParser(description="minSHmap - minimalistic sketch-based read mapper")
    p.add_argument("reference", nargs="?", help="reference FASTA")
    p.add_argument("reads", nargs="?", help="reads FASTA")
    p.add_argument("-k", type=int, default=15, help="k-mer length (default 15)")
    p.add_argument("-r", "--hfrac", type=float, default=0.05, help="FracMinHash ratio (default 0.05)")
    p.add_argument("-t", "--theta", type=float, default=0.9, help="homology threshold (default 0.9)")
    p.add_argument("-d", "--min-diff", type=float, default=0.02, help="best vs 2nd-best margin (default 0.02)")
    p.add_argument("-M", "--max-matches", type=int, default=0,
                   help="ignore k-mers occurring > this many times in the reference (0 = off)")
    p.add_argument("--hash", choices=SKETCHERS, default="nthash", help="sketch hash (default nthash)")
    p.add_argument("--demo", action="store_true", help="run a synthetic self-test")
    args = p.parse_args()

    if args.demo:
        demo(sketch=SKETCHERS[args.hash])
    elif args.reference and args.reads:
        map_reads(args.reference, args.reads, args.k, args.hfrac, args.theta,
                  args.min_diff, SKETCHERS[args.hash], max_matches=args.max_matches)
    else:
        p.error("provide REFERENCE and READS, or use --demo")


if __name__ == "__main__":
    main()
