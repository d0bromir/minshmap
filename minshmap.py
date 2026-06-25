"""minSHmap - the whole sketch-based read mapper in ~150 lines, nothing else.

The single, pedagogical version: one sketcher, no max_matches, no second-best/mapq
- we just keep the single best hit. `minshmap.cpp` is a port of this file (same
algorithm, same output) that additionally offers `-j` for output-identical
parallelism and a few cache/layout optimizations; keep the two in lockstep.

  1. SKETCH  - (w, k)-minimizers: slide a window of `w` consecutive k-mers and keep
               the one with the smallest hash (a sparse, well-spread sample of the
               sequence). Each base is 2 bits; forward and reverse-complement codes
               roll in O(1); the canonical (smaller) code is hashed, so either
               strand is found with one index.
  2. INDEX   - reference: hash -> [(segment, position, strand), ...].
  3. MAP     - sketch the read; rank its minimizers rarest-first; scatter the rarest
               ones into overlapping reference windows ("buckets"); for each
               window compute containment while pruning with the SEED HEURISTIC

                   sh = 1 - (seeds_used - matches) / m   (upper bound it can reach)

               the moment sh < theta the window is provably hopeless -> skip it.

Usage:  python minshmap.py reference.fa reads.fa -k 15 -w 10 -t 0.9
"""
from __future__ import annotations
import argparse

from Bio import SeqIO

MASK64 = (1 << 64) - 1
_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}     # 2 bits per base; complement(b) = 3 - b


def _mix(x):
    """SplitMix64 finalizer: a packed k-mer code -> a well-mixed 64-bit hash."""
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return x ^ (x >> 31)


def minimizers(seq, k, w):
    """(w, k)-minimizers via a rolling 2-bit canonical code. Yields (pos, hash, strand):
    for each window of w consecutive k-mers keep the smallest-hash one (leftmost on
    ties); each position is emitted once, in increasing order. k<=32. O(n*w)."""
    n = len(seq)
    if n < k:
        return
    mask = (1 << (2 * k)) - 1
    sh = 2 * (k - 1)
    fw = rc = 0
    for t in range(k):                       # pack the first window
        b = _CODE.get(seq[t], 0)
        fw = ((fw << 2) | b) & mask
        rc = (rc >> 2) | ((3 - b) << sh)
    kmers, i = [], 0                         # (hash, pos, strand) for every k-mer
    while True:
        c = fw if fw <= rc else rc           # canonical k-mer (min over both strands)
        kmers.append((_mix(c), i, 0 if fw <= rc else 1))
        if i + k >= n:
            break
        b = _CODE.get(seq[i + k], 0)         # roll the window one base
        fw = ((fw << 2) | b) & mask
        rc = (rc >> 2) | ((3 - b) << sh)
        i += 1
    ww, last = min(w, len(kmers)), -1
    for j in range(len(kmers) - ww + 1):     # slide the window, emit each new minimizer
        best = min(range(j, j + ww), key=lambda t: kmers[t][0])   # leftmost smallest hash
        if best != last:
            h, pos, strand = kmers[best]
            yield pos, h, strand
            last = best


def build_index(refs, k, w):
    """Reference minimizers inverted: hash -> [(segm_id, pos, strand), ...]."""
    index, names, lengths = {}, [], []
    for sid, (name, seq) in enumerate(refs):
        names.append(name); lengths.append(len(seq))
        for pos, h, strand in minimizers(seq, k, w):
            index.setdefault(h, []).append((sid, pos, strand))
    return index, names, lengths


def _seeds_rarest_first(sk, index):
    """One seed per distinct read minimizer: {hash: (ref_hits, read_strand)}, and the
    hashes ordered rarest (fewest reference hits) first."""
    seeds = {}
    for _, h, strand in sk:
        if h not in seeds:
            seeds[h] = (index.get(h, ()), strand)
    return seeds, sorted(seeds, key=lambda h: len(seeds[h][0]))


def _candidate_windows(seeds, order, m, theta, W):
    """Drop the rarest seeds' hits into overlapping buckets (b and b-1) so any read-sized
    region falls wholly inside one bucket; return the window keys by vote count (most first)."""
    S = int((1 - theta) * m) + 1             # this many rare seeds are enough to seed
    cand = {}
    for h in order[:S]:
        for sid, pos, _ in seeds[h][0]:
            b = pos // W
            cand[(sid, b)] = cand.get((sid, b), 0) + 1
            if b:
                cand[(sid, b - 1)] = cand.get((sid, b - 1), 0) + 1
    return sorted(cand, key=lambda key: (-cand[key], key))      # votes desc, key asc


def seed_heuristic(used, matches, m):
    """Upper bound on the containment a window can still reach, in [0, 1] (compare to theta):
    even if every remaining seed matched, the score cannot exceed this."""
    return 1 - (used - matches) / m


def _score_window(seeds, order, sid, lo, hi, m, target):
    """Add seeds rarest-first to window [lo, hi); stop as soon as the seed heuristic proves the
    window cannot reach `target`. Returns (score, codir, r_min, r_max), or None if pruned."""
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
                break                        # a read minimizer counts in a window once
        if seed_heuristic(used, matches, m) < target:
            return None                      # provably hopeless -> skip
    return matches / m, codir, r_min, r_max  # containment = fraction of read minimizers hit


def map_read(seq, index, k, w, theta):
    """Best reference window for the read: (segm, t_start, t_end, score, codir) or None.
    Stages run in order: sketch -> rank seeds -> candidate windows -> score each window."""
    sk = list(minimizers(seq, k, w))
    m = len(sk)                              # informative minimizers in the read
    if m == 0:
        return None
    W = max(len(seq), 1)                     # each candidate window is read-length-wide
    seeds, order = _seeds_rarest_first(sk, index)
    windows = _candidate_windows(seeds, order, m, theta, W)
    best = best_key = None
    for key in windows:                      # score promising windows first
        sid, b = key
        target = max(theta, best[3]) if best else theta   # can't beat best -> prune harder
        res = _score_window(seeds, order, sid, b * W, (b + 2) * W, m, target)
        if res is None:
            continue
        score, codir, r_min, r_max = res
        if score >= theta and (best is None or score > best[3]
                               or (score == best[3] and key < best_key)):  # tie-break on key
            best = (sid, r_min, r_max + k, score, codir)
            best_key = key
    return best


def read_fasta(path):
    """Yield (name, uppercase sequence) per FASTA record, parsed by Biopython."""
    for rec in SeqIO.parse(path, "fasta"):
        yield rec.id, str(rec.seq).upper()


def main():
    p = argparse.ArgumentParser(description="minSHmap - minimal sketch mapper")
    p.add_argument("reference"); p.add_argument("reads")
    p.add_argument("-k", type=int, default=15)
    p.add_argument("-w", "--window", type=int, default=10)
    p.add_argument("-t", "--theta", type=float, default=0.9)
    a = p.parse_args()
    index, names, lengths = build_index(read_fasta(a.reference), a.k, a.window)
    for name, seq in read_fasta(a.reads):
        mp = map_read(seq, index, a.k, a.window, a.theta)
        if mp is None:
            continue
        sid, ts, te, score, codir = mp
        nmatch = round(score * len(seq))
        print("\t".join(str(x) for x in (
            name, len(seq), 0, len(seq), "+" if codir >= 0 else "-",
            names[sid], lengths[sid], ts, te, nmatch, te - ts, 60)))


if __name__ == "__main__":
    main()
