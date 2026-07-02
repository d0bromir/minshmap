"""minSHmap - the whole sketch-based read mapper, nothing else.

Pedagogical version: one sketcher, a single best hit plus a phi-FREE mapping quality.
Three stages:

  1. SKETCH - canonical (w, k)-minimizers from the `minimizer-iter` library
              (rust-seq / Igor Martayan), exposed to Python by the tiny PyO3 wrapper
              in `minimizer_ext/`. We do NOT hand-roll the sketch: the rolling hash,
              sliding-window minimum and canonicalization all live in that library, so
              a read and its reverse complement select the same minimizers (one index
              serves both strands). `w` must be ODD (the canonical scheme uses an odd
              window to break forward/reverse ties).
  2. INDEX  - reference minimizers inverted: hash -> [(segment, position, strand), ...].
  3. MAP    - sketch the read, rank minimizers rarest-first, scatter them into
              overlapping blocks, score containment, prune with the SEED HEURISTIC
              sh = 1 - (seeds_used - matches) / m (best it can reach): sh < theta -> skip.
              Then a mapping quality mapq in {0, 60}: 60 iff the best mapping is
              confidently unique. Uniqueness uses the phi-FREE rule (Def. 5'): an
              *alternative* is any candidate whose reference interval is DISJOINT from
              the best one. Overlap is measured in reference coordinates, so the
              reverse-complement image of the best coincides with it and is not a
              spurious alternative -- no maximal-overlap parameter phi is needed.
              mapq = 60 iff every disjoint alternative is weaker by more than delta.

Setup once:  cd minimizer_ext && maturin build --release -i <python> && pip install <wheel>
Usage:       python minshmap.py reference.fa reads.fa -k 15 -w 11 -t 0.9 -d 0.15
"""
import argparse
from bisect import bisect_left

from Bio import SeqIO
from minimizer_ext import canonical_minimizers


def minimizers(seq, k, w):
    """Canonical (w, k)-minimizers straight from the minimizer-iter library: yields
    (pos, hash, strand) per emitted minimizer (`w` must be odd). The only thing we do
    here is reshape the library's tuple - it returns the strand as a bool, the mapper
    wants 0/1 - so none of the sketch logic is ours."""
    for pos, h, strand in canonical_minimizers(seq, k, w):
        yield pos, h, 1 if strand else 0


def build_index(refs, k, w):
    """Reference minimizers inverted: hash -> [(segm_id, pos, strand), ...]. Segments are
    indexed in order and positions ascend within each, so every hit list comes out sorted
    by (segm_id, pos) -- which is exactly what _score_block's binary search relies on."""
    index, names, lengths = {}, [], []
    for sid, (name, seq) in enumerate(refs):
        names.append(name); lengths.append(len(seq))
        for pos, h, strand in minimizers(seq, k, w):
            index.setdefault(h, []).append((sid, pos, strand))
    return index, names, lengths


def _seeds_rarest_first(sk, index):
    """One seed per distinct read minimizer: {hash: [ref_hits, read_strand, mult]}, mult =
    its multiplicity in the read (for weighted containment); hashes ordered rarest-first."""
    seeds = {}
    for _, h, strand in sk:                  # tally per-minimizer multiplicity (strand = first seen)
        seeds.setdefault(h, [index.get(h, ()), strand, 0])[2] += 1
    return seeds, sorted(seeds, key=lambda h: len(seeds[h][0]))


def _candidate_blocks(seeds, order, m, theta, B):
    """Drop the rarest seeds' hits into overlapping blocks (b and b-1) so any read-sized
    region falls wholly inside one block; return the block keys by vote count (most first)."""
    S = int((1 - theta) * m) + 1             # this many rare seeds are enough to seed
    cand = {}
    for h in order[:S]:
        for sid, pos, _ in seeds[h][0]:
            b = pos // B
            cand[(sid, b)] = cand.get((sid, b), 0) + 1
            if b:
                cand[(sid, b - 1)] = cand.get((sid, b - 1), 0) + 1
    return sorted(cand, key=lambda key: (-cand[key], key))      # votes desc, key asc


def seed_heuristic(used, matches, m):
    """Upper bound on the containment a block can still reach, in [0, 1] (compare to theta)."""
    return 1 - (used - matches) / m


def _score_block(seeds, order, sid, lo, hi, m, target):
    """Add seeds rarest-first to block [lo, hi); stop as soon as the seed heuristic proves the
    block cannot reach `target`. Returns (score, codir, r_min, r_max), or None if pruned.
    Each seed's smallest-pos hit in the block is found by BINARY SEARCH (hit lists are sorted
    by (sid, pos)), so a frequent minimizer costs O(log hits) instead of O(hits) -- identical
    result, same algorithm as the C++ port (minshmap.cpp::score_block)."""
    used = matches = codir = 0; r_min = r_max = -1
    for h in order:
        hits, rstrand, cnt = seeds[h]         # cnt = read multiplicity (weighted containment)
        used += cnt
        a = bisect_left(hits, (sid, lo))     # first hit with (s2, pos) >= (sid, lo)
        if a < len(hits) and hits[a][0] == sid and hits[a][1] < hi:
            pos, hstrand = hits[a][1], hits[a][2]
            matches += cnt
            codir += 1 if hstrand == rstrand else -1
            r_min = pos if r_min < 0 else min(r_min, pos)
            r_max = max(r_max, pos)
        if seed_heuristic(used, matches, m) < target:
            return None                      # provably hopeless -> skip
    return matches / m, codir, r_min, r_max  # weighted containment = matched occurrences / m


def _disjoint(a, b):
    """True iff mappings a, b = (sid, t_start, t_end) cover disjoint reference intervals."""
    return a[0] != b[0] or a[2] <= b[1] or b[2] <= a[1]


def map_read(seq, index, k, w, theta, delta=0.15):
    """Best reference block for the read plus a phi-free mapq:
    (segm, t_start, t_end, score, codir, mapq) or None. The placement (everything but
    mapq) is exactly the single best block; mapq is 60 iff every *alternative* mapping
    -- any scored candidate whose reference interval is DISJOINT from the best (Def. 5')
    -- is weaker than the best by more than delta (Def. 6'). No max-overlap phi."""
    sk = list(minimizers(seq, k, w))
    m = len(sk)                              # informative minimizers in the read
    if m == 0:
        return None
    B = max(len(seq), 1)                     # each candidate block is read-length-wide
    seeds, order = _seeds_rarest_first(sk, index)
    blocks = _candidate_blocks(seeds, order, m, theta, B)
    best = best_key = None
    cands = []                               # scored blocks that may be best OR a near-best alternative
    for key in blocks:                       # score promising blocks first
        sid, b = key
        # keep blocks within delta of the best so a disjoint second-best survives the prune
        target = max(theta, best[3] - delta) if best else theta
        res = _score_block(seeds, order, sid, b * B, (b + 2) * B, m, target)
        if res is None:
            continue
        score, codir, r_min, r_max = res
        if score < theta:
            continue
        mapping = (sid, r_min, r_max + k, score, codir)
        cands.append(mapping)
        if best is None or score > best[3] or (score == best[3] and key < best_key):  # tie-break on key
            best, best_key = mapping, key
    if best is None:
        return None
    second = 0.0                             # strongest alternative DISJOINT from the best
    for mp in cands:
        if mp is not best and _disjoint(mp, best):
            second = max(second, mp[3])
    mapq = 60 if second < best[3] - delta else 0
    return (*best, mapq)


def fasta(path):
    """(name, uppercase sequence) per record from Biopython; upper because the encoder only knows ACGT."""
    return ((rec.id, str(rec.seq).upper()) for rec in SeqIO.parse(path, "fasta"))


def main():
    p = argparse.ArgumentParser(description="minSHmap - minimal sketch mapper")
    p.add_argument("reference"); p.add_argument("reads")
    p.add_argument("-k", type=int, default=15)
    p.add_argument("-w", "--window", type=int, default=11)   # must be odd (canonical minimizers)
    p.add_argument("-t", "--theta", type=float, default=0.9)
    p.add_argument("-d", "--delta", type=float, default=0.15)  # mapq similarity margin (Def. 6')
    a = p.parse_args()
    if a.window % 2 == 0:
        p.error("window (-w) must be odd for canonical minimizers")
    index, names, lengths = build_index(fasta(a.reference), a.k, a.window)
    for name, seq in fasta(a.reads):
        mp = map_read(seq, index, a.k, a.window, a.theta, a.delta)
        if mp is None:
            continue
        sid, ts, te, score, codir, mapq = mp
        nmatch = round(score * len(seq))
        print("\t".join(str(x) for x in (
            name, len(seq), 0, len(seq), "+" if codir >= 0 else "-",
            names[sid], lengths[sid], ts, te, nmatch, te - ts, mapq)))


if __name__ == "__main__":
    main()
