"""Benchmark the two minSHmap implementations head-to-head: Python vs C++.

Both do the exact same job (map a read -> best reference window, PAF) and now use
the SAME minimizer library (minimizer-iter), so this is a fair apples-to-apples
race. minSH is intentionally NOT here: it solves a different task (global pairwise
alignment vs our semi-global mapping), so comparing the two is not meaningful.

For each dataset we report, per implementation: index build time, map throughput
(reads/s), how many reads mapped, and -- when reads carry truth in the header
(segm=/pos=/strand=) -- placement precision (mapped reads landing at the true locus).

Modes:
  synthetic : data/ref.fa + data/reads.fa (short ~300 bp reads, 5% error). Fast.
  chr21     : data_rw/chr21.fa + data_rw/{hifi,ont,clr}.fa (real long reads).
              Pure-Python mapping is ~2 reads/s; use --max-reads for a subset.

Usage:
  python 09_bench_py_vs_cpp.py synthetic
  python 09_bench_py_vs_cpp.py chr21 --datasets hifi ont clr --max-reads 200
"""
import argparse
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)
sys.path.insert(0, MINSHMAP_ROOT)            # for `import minshmap`

import minshmap as M                          # noqa: E402

RESULTS = os.path.join(HERE, "results_rw")
DATA_RW = os.path.join(HERE, "data_rw")
SMALL = os.path.join(MINSHMAP_ROOT, "data")
CPP_BIN = os.path.join(MINSHMAP_ROOT, "minshmap.exe" if os.name == "nt" else "minshmap_linux")


def parse_truth(header):
    """read0001 segm=chr1 pos=36202 strand=- len=301 -> dict, or {} if absent."""
    out = {}
    for tok in header.split()[1:]:
        if "=" in tok:
            key, val = tok.split("=", 1)
            out[key] = val
    return out


def read_truth(path):
    truth = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                toks = line[1:].split()
                truth[toks[0]] = parse_truth(line[1:])
    return truth


def bench_py(ref_path, reads_path, k, w, theta, max_reads):
    """Python minSHmap: returns (index_s, map_s, n, mapped, placements-by-id)."""
    refs = list(M.fasta(ref_path))
    reads = list(M.fasta(reads_path))
    if max_reads:
        reads = reads[:max_reads]
    t0 = time.perf_counter()
    index, names, lengths = M.build_index(refs, k, w)
    idx_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    placements = {rid: M.map_read(seq, index, k, w, theta) for rid, seq in reads}
    map_s = time.perf_counter() - t0
    mapped = sum(1 for mp in placements.values() if mp is not None)
    return idx_s, map_s, len(reads), mapped, names, placements


def py_precision(placements, names, reads_path, max_reads):
    """Of mapped reads with header truth, fraction at the true locus (seg+strand+pos)."""
    truth = read_truth(reads_path)
    correct = total = 0
    for i, (rid, mp) in enumerate(placements.items()):
        if max_reads and i >= max_reads:
            break
        t = truth.get(rid, {})
        if not t or "pos" not in t or mp is None:
            continue
        total += 1
        sid, ts, _te, _s, codir = mp
        if names[sid] == t.get("segm") and ("+" if codir >= 0 else "-") == t.get("strand"):
            correct += 1 if abs(ts - int(t["pos"])) <= 1000 else 0
    return (correct / total) if total else None, total


def bench_cpp(ref_path, reads_path, k, w, theta, max_reads):
    """C++ minSHmap: returns (map_s, mapped). Index+map not separable via CLI -> one wall."""
    cap = reads_path
    if max_reads:
        cap = os.path.join(RESULTS, "_cap.fa")
        with open(reads_path) as fi, open(cap, "w") as fo:
            n = 0
            for line in fi:
                if line.startswith(">") and (n := n + 1) > max_reads:
                    break
                fo.write(line)
    t0 = time.perf_counter()
    out = subprocess.run([CPP_BIN, ref_path, cap, "-k", str(k), "-w", str(w), "-t", str(theta)],
                         capture_output=True, text=True)
    wall = time.perf_counter() - t0
    mapped = len({l.split("\t")[0] for l in out.stdout.splitlines() if len(l.split("\t")) >= 9})
    if cap != reads_path:
        os.remove(cap)
    return wall, mapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["synthetic", "chr21"])
    ap.add_argument("--datasets", nargs="+", default=["hifi", "ont", "clr"])
    ap.add_argument("--max-reads", type=int, default=0, help="0 = no cap")
    ap.add_argument("-k", type=int, default=15)
    ap.add_argument("-w", "--window", type=int, default=11)   # odd: canonical minimizers
    args = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)

    if args.mode == "synthetic":
        jobs = [("synthetic", os.path.join(SMALL, "ref.fa"), os.path.join(SMALL, "reads.fa"), 0.5)]
    else:
        ref = os.path.join(DATA_RW, "chr21.fa")
        theta = {"hifi": 0.20, "ont": 0.15, "clr": 0.18}
        jobs = [(ds, ref, os.path.join(DATA_RW, f"{ds}.fa"), theta.get(ds, 0.2)) for ds in args.datasets]

    rows = []
    for label, ref, reads, theta in jobs:
        if not os.path.exists(reads):
            print(f"skip {label}: {reads} missing", flush=True)
            continue
        idx_s, map_s, n, py_mapped, names, plc = bench_py(ref, reads, args.k, args.window, theta, args.max_reads)
        prec, n_truth = py_precision(plc, names, reads, args.max_reads)
        cpp_wall, cpp_mapped = (None, None)
        if os.path.exists(CPP_BIN):
            cpp_wall, cpp_mapped = bench_cpp(ref, reads, args.k, args.window, theta, args.max_reads)
        print(f"\n=== {label} (k={args.k} w={args.window} t={theta}, {n} reads) ===")
        print(f"  py  index {idx_s:7.2f}s  map {map_s:7.2f}s ({n/map_s:7.1f} r/s)  mapped {py_mapped}"
              + (f"  precision {100*prec:.1f}% of {n_truth}" if prec is not None else ""))
        if cpp_wall:
            print(f"  cpp total {cpp_wall:7.2f}s ({n/cpp_wall:7.1f} r/s)  mapped {cpp_mapped}")
        rows.append({"dataset": label, "reads": n, "theta": theta,
                     "py_index_s": round(idx_s, 3), "py_map_s": round(map_s, 3),
                     "py_reads_per_s": round(n/map_s, 2), "py_mapped": py_mapped,
                     "py_precision": round(prec, 4) if prec is not None else "",
                     "cpp_wall_s": round(cpp_wall, 3) if cpp_wall else "",
                     "cpp_reads_per_s": round(n/cpp_wall, 2) if cpp_wall else "",
                     "cpp_mapped": cpp_mapped if cpp_mapped is not None else ""})

    stamp = time.strftime("%Y%m%d-%H%M%S")
    out = os.path.join(RESULTS, f"py_vs_cpp_{args.mode}_{stamp}.csv")
    with open(out, "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wri.writeheader()
        wri.writerows(rows)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
