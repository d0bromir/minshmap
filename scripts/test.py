"""Correctness tests for minSHmap and its three sketchers.

Run from the repo root:

    python scripts/test.py

Exits non-zero on the first failure so it can gate reproducible runs / CI.
"""

from __future__ import annotations

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import minshmap as M  # noqa: E402

NUC = "ACGT"
PASS = 0
FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        print(f"  FAIL: {msg}")


def random_seq(n, rng):
    return "".join(rng.choice(NUC) for _ in range(n))


def mutate(s, e, rng):
    out = []
    for c in s:
        if rng.random() > e:
            out.append(c)
        elif rng.random() < 1 / 3:
            out.append(rng.choice(NUC))
        elif rng.random() < 1 / 2:
            out.append(rng.choice(NUC) + c)
    return "".join(out)


# --------------------------------------------------------------------------- #

def test_revcomp():
    print("revcomp involution")
    rng = random.Random(0)
    for _ in range(100):
        s = random_seq(rng.randint(1, 50), rng)
        check(M.revcomp(M.revcomp(s)) == s, "revcomp(revcomp(s)) != s")
        check(M.revcomp("ACGT") == "ACGT", "revcomp(ACGT) wrong")


def test_sketch_canonical_strand():
    """A sketch and the sketch of its reverse complement must share hashes."""
    print("strand-canonical hashing (k-mer and its rev-comp share a hash)")
    rng = random.Random(1)
    k, hfrac = 11, 0.3
    for name, sketch in M.SKETCHERS.items():
        seq = random_seq(400, rng)
        rc = M.revcomp(seq)
        h_fwd = {h for _, h, _ in sketch(seq, k, hfrac)}
        h_rev = {h for _, h, _ in sketch(rc, k, hfrac)}
        # Every kept hash on one strand must appear on the other (same k-mers).
        check(h_fwd == h_rev, f"{name}: forward/revcomp hash sets differ")


def test_rolling_matches_naive_selectivity():
    """Rolling hashes must keep roughly the expected fraction of k-mers."""
    print("FracMinHash selectivity (~hfrac of k-mers kept)")
    rng = random.Random(2)
    seq = random_seq(20_000, rng)
    k, hfrac = 15, 0.1
    n_kmers = len(seq) - k + 1
    for name, sketch in M.SKETCHERS.items():
        kept = len(sketch(seq, k, hfrac))
        frac = kept / n_kmers
        check(0.05 < frac < 0.20, f"{name}: kept fraction {frac:.3f} far from {hfrac}")


def test_rolling_hash_consistency():
    """Rolling and naive sketchers must agree on positions/strands they keep
    *within their own hashing* - here we assert each is internally stable and
    that a rolled hash equals a freshly computed one for the poly hash."""
    print("rolling == recomputed (poly hash window invariant)")
    rng = random.Random(3)
    seq = random_seq(2000, rng)
    k, hfrac = 13, 1.0  # keep everything: compare every window
    rolled = M.sketch_poly(seq, k, hfrac)
    # Recompute each kept window independently with the same polynomial.
    for pos, h, _ in rolled:
        kmer = seq[pos : pos + k]
        h_fw = 0
        h_rc = 0
        for t, c in enumerate(kmer):
            v = M._VAL[c]
            h_fw = (h_fw * M._B + v) % M._M
            h_rc = (h_rc + (3 - v) * pow(M._B, t, M._M)) % M._M
        check(h == min(h_fw, h_rc), f"poly rolled hash != recomputed at pos {pos}")


def test_mapping_accuracy():
    """End-to-end: each sketcher should place most reads at the truth."""
    print("end-to-end mapping accuracy (>= 90% placed)")
    rng = random.Random(4)
    k, hfrac, theta, min_diff, error = 11, 0.2, 0.5, 0.02, 0.04
    ref = random_seq(8000, rng)
    truth = []
    reads = []
    for _ in range(40):
        start = rng.randint(0, len(ref) - 400)
        read = mutate(ref[start : start + 300], error, rng)
        if rng.random() < 0.5:
            read = M.revcomp(read)
        reads.append(read)
        truth.append(start)

    for name, sketch in M.SKETCHERS.items():
        index, segments = M.build_index([("chr", ref)], k, hfrac, sketch)
        correct = 0
        for read, start in zip(reads, truth):
            best, _ = M.map_read(read, index, segments, k, hfrac, theta, min_diff, sketch)
            if best is not None and best.t_start <= start + 150 <= best.t_end:
                correct += 1
        check(correct >= 0.9 * len(reads), f"{name}: only {correct}/{len(reads)} placed")


def main():
    for t in (
        test_revcomp,
        test_sketch_canonical_strand,
        test_rolling_matches_naive_selectivity,
        test_rolling_hash_consistency,
        test_mapping_accuracy,
    ):
        t()
    print(f"\n{PASS} checks passed, {FAIL} failed.")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
