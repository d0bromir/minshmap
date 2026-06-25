"""Benchmark the *Python* minSHmap mapper and compare it with minSH (Python).

The two tools do different jobs, so this is not a head-to-head race:
  - minSHmap (minshmap.py) MAPS a read -> best reference window (a PAF location).
  - minSH (minsh.astar) ALIGNS two sequences -> exact edit distance via seed-A*.

So we use them in tandem, the standard way to compare a mapper with an aligner:
minSHmap places each read, then minSH aligns the read against that placed window
to score the placement (edit distance / identity). We report each tool's own
throughput plus minSH's identity, all on the same Python interpreter.

Modes:
  synthetic : data/ref.fa + data/reads.fa (short ~300 bp reads, 5% error). Both
              tools finish on every read; reads carry truth in the header
              (segm=/pos=/strand=) so we also report minSHmap mapping accuracy.
  chr21     : data_rw/chr21.fa + data_rw/{hifi,ont,clr}.fa (real long reads).
              Pure-Python mapping is ~2 reads/s and minSH's A* cannot finish
              16 kb reads, so this mode is for small subsets / a spot-check.

Usage:
  python 09_bench_py_vs_minsh.py synthetic
  python 09_bench_py_vs_minsh.py chr21 --datasets hifi ont clr --max-reads 100
"""
import argparse
import csv
import math
import os
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)
WORKSPACE = os.path.dirname(MINSHMAP_ROOT)
sys.path.insert(0, MINSHMAP_ROOT)            # for `import minshmap`
sys.path.insert(0, os.path.join(WORKSPACE, "minSH"))  # for `import minsh`

import minshmap as M                          # noqa: E402
from minsh.astar import build_seedh           # noqa: E402
from minsh.banded import align_banded         # noqa: E402

RESULTS = os.path.join(HERE, "results_rw")
DATA_RW = os.path.join(HERE, "data_rw")
SMALL = os.path.join(MINSHMAP_ROOT, "data")

LEN_CAP = 12000          # minSH A* is intractable on longer reads -> skip for alignment
BUDGET_FRAC = 0.40       # minSH banded budget = 40% edits of the read length

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def kfor(n):
    return max(3, math.ceil(math.log(max(n, 4), 4)))


def parse_truth(header):
    """read0001 segm=chr1 pos=36202 strand=- len=301 -> dict, or {} if absent."""
    out = {}
    for tok in header.split()[1:]:
        if "=" in tok:
            key, val = tok.split("=", 1)
            out[key] = val
    return out


def run_minshmap(refs, reads, k, w, theta):
    """Map every read. Returns (per-read placements, timings, counts)."""
    t0 = time.perf_counter()
    index, names, lengths = M.build_index(refs, k, w)
    idx_s = time.perf_counter() - t0
    name2id = {n: i for i, n in enumerate(names)}

    placements, t0 = [], time.perf_counter()
    for rid, seq in reads:
        mp = M.map_read(seq, index, k, w, theta)
        placements.append((rid, seq, mp))
    map_s = time.perf_counter() - t0
    return placements, names, lengths, name2id, idx_s, map_s


def score_with_minsh(placements, ref_by_id):
    """Align each mapped read against its placed locus with banded minSH.
    The placed interval is the minimizer span, slightly shorter than the read, so we
    align against a read-length reference window anchored at the placement start --
    the identity then reflects true read divergence, not the window's clipped ends.
    Returns identities, edit%, alignment wall time, and #skipped/#over-budget."""
    idents, edit_pcts, n_skip, n_fail = [], [], 0, 0
    t0 = time.perf_counter()
    for _, seq, mp in placements:
        if mp is None:
            continue
        sid, ts, _te, _score, codir = mp
        B = seq if codir >= 0 else revcomp(seq)
        ts = max(0, ts)
        A = ref_by_id[sid][ts: ts + len(B)]  # read-length reference window
        if not A or not B:
            continue
        if len(B) > LEN_CAP:
            n_skip += 1
            continue
        k = kfor(len(A))
        h = build_seedh(A, B, k)
        res = align_banded(A, B, h, max_cost=math.ceil(BUDGET_FRAC * len(B)),
                           return_stats=True)
        if res is None:                      # diverges beyond budget
            n_fail += 1
            continue
        _g, dist, _cells = res
        denom = max(len(A), len(B))
        idents.append(1 - dist / denom)
        edit_pcts.append(100 * dist / denom)
    align_s = time.perf_counter() - t0
    return idents, edit_pcts, align_s, n_skip, n_fail


def read_truth(path):
    """first-token -> {segm,pos,strand,len} from FASTA headers (synthetic reads)."""
    truth = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                toks = line[1:].split()
                truth[toks[0]] = parse_truth(line[1:])
    return truth


def mapping_accuracy(placements, names, lengths, truth_map):
    """Placement precision: of mapped reads that carry truth, how many landed at the
    true locus (segm + strand + within one read length of the true start).
    Returns (precision, n_mapped_with_truth)."""
    n_mapped_truth = correct = 0
    for rid, seq, mp in placements:
        t = truth_map.get(rid, {})
        if not t or "pos" not in t or mp is None:
            continue
        n_mapped_truth += 1
        if mp is None:
            continue
        sid, ts, te, _s, codir = mp
        ok_seg = names[sid] == t.get("segm")
        ok_str = ("+" if codir >= 0 else "-") == t.get("strand")
        tol = len(seq)                       # within one read length of the true start
        ok_pos = abs(ts - int(t["pos"])) <= tol
        if ok_seg and ok_str and ok_pos:
            correct += 1
    return (correct / n_mapped_truth) if n_mapped_truth else None, n_mapped_truth


def summarize(label, placements, names, lengths, ref_by_id, idx_s, map_s, truth_map):
    n = len(placements)
    mapped = [p for p in placements if p[2] is not None]
    nm = len(mapped)
    idents, edit_pcts, align_s, n_skip, n_fail = score_with_minsh(mapped, ref_by_id)
    acc, n_truth = mapping_accuracy(placements, names, lengths, truth_map)

    print(f"\n=== {label} ===")
    print(f"  reads                  : {n}")
    print(f"  minSHmap index build   : {idx_s:8.2f} s")
    print(f"  minSHmap map           : {map_s:8.2f} s  ({n / map_s:7.1f} reads/s)")
    print(f"  minSHmap mapped        : {nm} ({100 * nm / n:.1f}%)")
    if acc is not None:
        print(f"  minSHmap precision     : {100 * acc:.1f}%  (mapped reads at true locus, of {n_truth})")
    n_aln = len(idents)
    if n_aln:
        print(f"  minSH aligned          : {n_aln} reads"
              f"{f' (+{n_skip} too long, +{n_fail} over budget)' if (n_skip or n_fail) else ''}")
        print(f"  minSH align time       : {align_s:8.2f} s  ({n_aln / align_s:7.2f} reads/s)")
        print(f"  minSH median identity  : {statistics.median(idents):.4f}")
        print(f"  minSH median edit rate : {statistics.median(edit_pcts):.2f}%")
    else:
        print(f"  minSH aligned          : 0 reads"
              f" ({n_skip} too long, {n_fail} over budget) - long-read A* intractable")

    return {
        "label": label, "reads": n,
        "minshmap_index_s": round(idx_s, 3), "minshmap_map_s": round(map_s, 3),
        "minshmap_reads_per_s": round(n / map_s, 2),
        "minshmap_mapped": nm, "minshmap_mapped_frac": round(nm / n, 4),
        "minshmap_precision": round(acc, 4) if acc is not None else "",
        "minsh_aligned": n_aln, "minsh_skipped_long": n_skip, "minsh_over_budget": n_fail,
        "minsh_align_s": round(align_s, 3),
        "minsh_reads_per_s": round(n_aln / align_s, 4) if n_aln else "",
        "minsh_median_identity": round(statistics.median(idents), 4) if idents else "",
        "minsh_median_edit_pct": round(statistics.median(edit_pcts), 3) if edit_pcts else "",
    }


def load_fasta_list(path):
    return list(M.fasta(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["synthetic", "chr21"])
    ap.add_argument("--datasets", nargs="+", default=["hifi", "ont", "clr"])
    ap.add_argument("--max-reads", type=int, default=0, help="0 = no cap")
    ap.add_argument("-k", type=int, default=15)
    ap.add_argument("-w", "--window", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(RESULTS, exist_ok=True)
    rows = []

    if args.mode == "synthetic":
        refs = load_fasta_list(os.path.join(SMALL, "ref.fa"))
        ref_by_id = {i: seq for i, (_n, seq) in enumerate(refs)}
        reads = load_fasta_list(os.path.join(SMALL, "reads.fa"))
        if args.max_reads:
            reads = reads[: args.max_reads]
        theta = 0.5                          # same operating point as the parity test
        truth_map = read_truth(os.path.join(SMALL, "reads.fa"))
        placements, names, lengths, _n2i, idx_s, map_s = run_minshmap(
            refs, reads, args.k, args.window, theta)
        rows.append(summarize(f"synthetic (k={args.k} w={args.window} t={theta})",
                              placements, names, lengths, ref_by_id, idx_s, map_s, truth_map))
    else:
        theta_preset = {"hifi": 0.20, "ont": 0.15, "clr": 0.18}
        refs = load_fasta_list(os.path.join(DATA_RW, "chr21.fa"))
        ref_by_id = {i: seq for i, (_n, seq) in enumerate(refs)}
        for ds in args.datasets:
            reads = load_fasta_list(os.path.join(DATA_RW, f"{ds}.fa"))
            if args.max_reads:
                reads = reads[: args.max_reads]
            theta = theta_preset.get(ds, 0.2)
            placements, names, lengths, _n2i, idx_s, map_s = run_minshmap(
                refs, reads, args.k, args.window, theta)
            rows.append(summarize(f"chr21/{ds} (k={args.k} w={args.window} t={theta})",
                                  placements, names, lengths, ref_by_id, idx_s, map_s, {}))

    stamp = time.strftime("%Y%m%d-%H%M%S")
    out = os.path.join(RESULTS, f"py_vs_minsh_{args.mode}_{stamp}.csv")
    with open(out, "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wri.writeheader()
        wri.writerows(rows)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
