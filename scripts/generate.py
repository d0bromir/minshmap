"""Generate a deterministic reference + reads FASTA pair for testing/benchmarks.

Truth is encoded in each read header so mappings can be scored:

    >read000 segm=chr1 pos=1234 strand=+ len=300

Run from the repo root:

    python scripts/generate.py                 # writes data/ref.fa, data/reads.fa
    python scripts/generate.py --reads 2000    # bigger read set
"""

from __future__ import annotations

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from minshmap import revcomp  # noqa: E402

NUC = "ACGT"


def random_seq(n, rng):
    return "".join(rng.choice(NUC) for _ in range(n))


def mutate(s, e, rng):
    """Apply substitution/insertion/deletion errors at per-base rate `e`."""
    out = []
    for c in s:
        if rng.random() > e:
            out.append(c)
        elif rng.random() < 1 / 3:
            out.append(rng.choice(NUC))          # substitution
        elif rng.random() < 1 / 2:
            out.append(rng.choice(NUC) + c)      # insertion
        # else: deletion
    return "".join(out)


def write_fasta(path, records, width=70):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i : i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description="generate synthetic reference + reads")
    ap.add_argument("--ref-len", type=int, default=50_000, help="reference length per segment")
    ap.add_argument("--segments", type=int, default=2, help="number of reference segments")
    ap.add_argument("--reads", type=int, default=500, help="number of reads")
    ap.add_argument("--read-len", type=int, default=300, help="read length")
    ap.add_argument("--error", type=float, default=0.05, help="per-base error rate")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed (reproducible)")
    ap.add_argument("--outdir", default="data", help="output directory")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    segments = [(f"chr{s + 1}", random_seq(args.ref_len, rng)) for s in range(args.segments)]
    write_fasta(os.path.join(args.outdir, "ref.fa"), segments)

    reads = []
    for i in range(args.reads):
        sid = rng.randrange(args.segments)
        name, seq = segments[sid]
        start = rng.randint(0, len(seq) - args.read_len)
        sub = seq[start : start + args.read_len]
        strand = "+"
        read = mutate(sub, args.error, rng)
        if rng.random() < 0.5:
            read = revcomp(read)
            strand = "-"
        header = f"read{i:04d} segm={name} pos={start} strand={strand} len={len(read)}"
        reads.append((header, read))
    write_fasta(os.path.join(args.outdir, "reads.fa"), reads)

    print(f"Wrote {args.segments} segment(s) x {args.ref_len} bp -> {args.outdir}/ref.fa")
    print(f"Wrote {args.reads} reads (len {args.read_len}, error {args.error:.0%}) -> {args.outdir}/reads.fa")


if __name__ == "__main__":
    main()
