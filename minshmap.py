"""minSHmap - the whole sketch-based read mapper, nothing else.

Pedagogical version: one sketcher, a single best hit (no max_matches / second-best /
mapq). Three stages:

  1. SKETCH - canonical (w, k)-minimizers from the `minimizer-iter` library
              (rust-seq / Igor Martayan), exposed to Python by the tiny PyO3 wrapper
              in `minimizer_ext/`. We do NOT hand-roll the sketch: the rolling hash,
              sliding-window minimum and canonicalization all live in that library, so
              a read and its reverse complement select the same minimizers (one index
              serves both strands). `w` must be ODD (the canonical scheme uses an odd
              window to break forward/reverse ties).
  2. INDEX  - reference minimizers inverted: hash -> [(segment, position, strand), ...].
  3. MAP    - sketch the read, rank minimizers rarest-first, scatter them into
              overlapping windows, score containment, prune with the SEED HEURISTIC
              sh = 1 - (seeds_used - matches) / m (best it can reach): sh < theta -> skip.

Setup once:  cd minimizer_ext && maturin build --release -i <python> && pip install <wheel>
Usage:       python minshmap.py reference.fa reads.fa -k 15 -w 11 -t 0.9
"""
import argparse

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
        if h not in seeds:                   # first occurrence wins; skip duplicate lookups
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
    """Upper bound on the containment a window can still reach, in [0, 1] (compare to theta)."""
    return 1 - (used - matches) / m


def _score_window(seeds, order, sid, lo, hi, m, target):
    """Add seeds rarest-first to window [lo, hi); stop as soon as the seed heuristic proves the
    window cannot reach `target`. Returns (score, codir, r_min, r_max), or None if pruned."""
    used = matches = codir = 0; r_min = r_max = -1
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


def fasta(path):
    """(name, uppercase sequence) per record from Biopython; upper because the encoder only knows ACGT."""
    return ((rec.id, str(rec.seq).upper()) for rec in SeqIO.parse(path, "fasta"))


def main():
    p = argparse.ArgumentParser(description="minSHmap - minimal sketch mapper")
    p.add_argument("reference"); p.add_argument("reads")
    p.add_argument("-k", type=int, default=15)
    p.add_argument("-w", "--window", type=int, default=11)   # must be odd (canonical minimizers)
    p.add_argument("-t", "--theta", type=float, default=0.9)
    a = p.parse_args()
    if a.window % 2 == 0:
        p.error("window (-w) must be odd for canonical minimizers")
    index, names, lengths = build_index(fasta(a.reference), a.k, a.window)
    for name, seq in fasta(a.reads):
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
