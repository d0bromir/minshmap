#!/usr/bin/env python3
"""Real-world benchmark: shmap (FracMinHash) vs minSHmap C++ vs minSHmap Python,
all run inside WSL/Linux so timings share one environment.

Two SCOPES, following the rule "map reads against the reference they came from":
  * wholegenome : the reads (hifi/ont/clr) are whole-genome HG002 WGS, so they are
                  mapped against the WHOLE human genome (hs1.fa = T2T-CHM13v2.0).
                  This is the correct pairing and reaches ~100% for accurate reads.
                  Python is skipped here (a whole-genome index does not fit in RAM).
  * chromosome  : from the whole-genome run we take the reads that truly belong to
                  chr21 (cpp mapped them to chr21), then map THOSE reads against
                  chr21 only. This is the chromosome<->chromosome pairing, so it too
                  should map ~all of them. All three mappers participate.

Equal sensitivity: minSHmap keeps ~2/(w+1) of the k-mers as (w,k)-minimizers, so
shmap's FracMinHash ratio -r is set to 2/(w+1); both share -k and -t; shmap uses
the Containment metric to match minSHmap's containment score.

Per run we report: total time, INDEX time and MAP time separately, index resident
memory and peak memory, reads_in, mapped, map%, mean mapq, and map throughput.
minSHmap prints these on stderr when MINSHMAP_BENCH is set (stdout PAF stays
byte-identical); shmap (never modified) is wrapped in /usr/bin/time -v for peak RSS.

Run from WSL:
    PYTHONPATH=~/pylib python3 realworld/11_bench_3way.py --scope both
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                       # minshmap/
DATA = HERE / "data_rw"
RESULTS = HERE / "results_rw"
SHMAP = (ROOT.parent / "shmap" / "release" / "shmap").resolve()
CPP = (ROOT / "minshmap_linux").resolve()
PYSCRIPT = (ROOT / "minshmap.py").resolve()
PYLIB = Path(os.path.expanduser("~/pylib"))
CHR_REF = DATA / "chr21.fa"
WG_REF = DATA / "hs1.fa"
CHR_NAME = "chr21"                        # segment name inside hs1.fa / chr21.fa

# Per-dataset homology threshold (long-read error rates differ).
THETA = {"hifi": 0.20, "ont": 0.15, "clr": 0.18}

_BENCH_RE = re.compile(
    r"index_s=(\S+)\s+map_s=(\S+)\s+reads=(\S+)\s+mapped=(\S+)\s+"
    r"index_rss_mb=(\S+)\s+peak_rss_mb=(\S+)")
_TIME_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")


def subset_reads(src: Path, dst: Path, max_reads: int) -> int:
    """Write a UNIFORM sample of `max_reads` records across the whole file (0 = all).
    Even spacing matters because SRA read files are often position-grouped, so the
    first-N reads may miss whole chromosomes; a strided sample stays representative.
    Returns the number of records written."""
    total = 0
    with open(src) as f:
        for line in f:
            if line.startswith(">"):
                total += 1
    if max_reads <= 0 or max_reads >= total:
        stride, want = 1, total
    else:
        stride, want = total // max_reads, max_reads
    idx, written, keep = -1, 0, False
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                idx += 1
                keep = (idx % stride == 0) and written < want
                if keep:
                    written += 1
            if keep:
                fo.write(line)
    return written


def extract_reads(src: Path, dst: Path, names: set) -> int:
    """Write the records of src whose name (first header token) is in `names`."""
    n = 0
    keep = False
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                nm = line[1:].split()[0]
                keep = nm in names
                if keep:
                    n += 1
            if keep:
                fo.write(line)
    return n


def parse_paf(stdout: str, seg: str = None):
    """(#mapped unique queries, mean mapq, {names mapped to `seg`}). Column 12 is
    mapq for both tools; shmap appends extra tags after it (ignored). Column 6 is
    the target segment name."""
    names, mapqs, on_seg = set(), [], set()
    for line in stdout.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        names.add(f[0])
        try:
            mapqs.append(int(f[11]))
        except ValueError:
            pass
        if seg is not None and f[5] == seg:
            on_seg.add(f[0])
    mean_mapq = sum(mapqs) / len(mapqs) if mapqs else 0.0
    return len(names), mean_mapq, on_seg


def run_mapper(mapper, ref, reads_file, k, w, theta, density, threads):
    """Run one mapper; return dict(stdout, total_s, index_s, map_s, index_rss_mb,
    peak_rss_mb). index/map split comes from minSHmap's [bench] stderr; for shmap
    only the peak RSS (via /usr/bin/time -v) is available (index/map split = None)."""
    env = dict(os.environ)
    if mapper == "shmap":
        env["OMP_NUM_THREADS"] = str(threads)
        cmd = ["/usr/bin/time", "-v", str(SHMAP), "-p", str(reads_file), "-s", str(ref),
               "-k", str(k), "-r", f"{density:.6f}", "-t", str(theta), "-m", "Containment"]
    elif mapper == "cpp":
        env["MINSHMAP_BENCH"] = "1"
        cmd = [str(CPP), str(ref), str(reads_file),
               "-k", str(k), "-w", str(w), "-t", str(theta), "-j", str(threads)]
    else:  # py
        env["MINSHMAP_BENCH"] = "1"
        env["PYTHONPATH"] = str(PYLIB) + os.pathsep + env.get("PYTHONPATH", "")
        cmd = [sys.executable, str(PYSCRIPT), str(ref), str(reads_file),
               "-k", str(k), "-w", str(w), "-t", str(theta)]

    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    total_s = time.perf_counter() - t0
    if r.returncode != 0:
        sys.stderr.write(f"  ! {mapper} rc={r.returncode}: {r.stderr.strip()[:300]}\n")

    res = {"stdout": r.stdout, "total_s": total_s,
           "index_s": None, "map_s": None, "index_rss_mb": None, "peak_rss_mb": None}
    if mapper == "shmap":
        mm = _TIME_RSS_RE.search(r.stderr)
        if mm:
            res["peak_rss_mb"] = int(mm.group(1)) / 1024.0
    else:
        mm = _BENCH_RE.search(r.stderr)
        if mm:
            res["index_s"] = float(mm.group(1))
            res["map_s"] = float(mm.group(2))
            res["index_rss_mb"] = float(mm.group(5))
            res["peak_rss_mb"] = float(mm.group(6))
    return res


def make_row(scope, ds, mapper, ref, reads_in, res, density, theta, k, w, threads):
    mapped, mmapq, _ = parse_paf(res["stdout"])
    pct = 100.0 * mapped / reads_in if reads_in else 0.0
    map_s = res["map_s"] if res["map_s"] is not None else res["total_s"]
    rps = reads_in / map_s if map_s and map_s > 0 else 0.0
    return {
        "scope": scope, "dataset": ds, "mapper": mapper, "ref": Path(ref).name,
        "reads_in": reads_in, "mapped": mapped, "map_pct": round(pct, 2),
        "mean_mapq": round(mmapq, 1),
        "total_s": round(res["total_s"], 2),
        "index_s": round(res["index_s"], 2) if res["index_s"] is not None else "",
        "map_s": round(res["map_s"], 2) if res["map_s"] is not None else "",
        "index_rss_mb": round(res["index_rss_mb"], 1) if res["index_rss_mb"] is not None else "",
        "peak_rss_mb": round(res["peak_rss_mb"], 1) if res["peak_rss_mb"] is not None else "",
        "map_reads_per_s": round(rps, 1),
        "k": k, "w": w, "density_or_r": f"{density:.4f}", "theta": theta, "threads": threads,
    }


def print_row(r):
    isec = r["index_s"] if r["index_s"] != "" else "-"
    msec = r["map_s"] if r["map_s"] != "" else "-"
    imem = r["index_rss_mb"] if r["index_rss_mb"] != "" else "-"
    pmem = r["peak_rss_mb"] if r["peak_rss_mb"] != "" else "-"
    print(f"    {r['mapper']:5s} mapped={r['mapped']:<6d} ({r['map_pct']:5.2f}%) "
          f"mapq={r['mean_mapq']:5.1f}  total={r['total_s']:8.2f}s "
          f"index={isec}s map={msec}s  idx_mem={imem}MB peak={pmem}MB")


def main():
    ap = argparse.ArgumentParser(description="real-world benchmark (chromosome + whole genome)")
    ap.add_argument("--scope", choices=["chromosome", "wholegenome", "both"], default="both")
    ap.add_argument("--datasets", nargs="+", default=["hifi", "ont", "clr"])
    ap.add_argument("-k", type=int, default=15)
    ap.add_argument("-w", type=int, default=31, help="minimizer window (odd)")
    ap.add_argument("--max-reads", type=int, default=2000,
                    help="whole-genome reads per dataset to feed in (0 = all)")
    ap.add_argument("--threads", type=int, default=1,
                    help="cpp -j / shmap OMP_NUM_THREADS (py is single-threaded)")
    a = ap.parse_args()
    if a.w % 2 == 0:
        ap.error("-w must be odd")

    density = 2.0 / (a.w + 1)            # FracMinHash ratio matching minimizer density
    RESULTS.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    csv_path = RESULTS / f"bench_{a.scope}_{stamp}.csv"
    need_wg = a.scope in ("wholegenome", "both")
    need_chr = a.scope in ("chromosome", "both")

    print(f"scope={a.scope}  k={a.k}  w={a.w}  shmap -r={density:.4f}  "
          f"threads={a.threads}  max_reads={a.max_reads}")
    print(f"binaries : shmap={SHMAP.exists()}  cpp={CPP.exists()}  "
          f"py_lib={(PYLIB / 'minimizer_ext.so').exists()}  "
          f"wg_ref={WG_REF.exists()}  chr_ref={CHR_REF.exists()}\n")

    rows = []
    for ds in a.datasets:
        src = DATA / f"{ds}.fa"
        if not src.exists():
            print(f"[skip] {ds}: {src} not found")
            continue
        theta = THETA.get(ds, 0.15)
        subset = RESULTS / f"_subset_{ds}.fa"
        n_wg = subset_reads(src, subset, a.max_reads)
        print(f"### {ds}  (theta={theta}, whole-genome reads={n_wg})")

        # ---- whole-genome scope (also identifies the genuine chr21 reads) ----
        # cpp always runs (its PAF identifies the chr21 reads); shmap only when the
        # whole-genome scope is actually reported (it is by far the slowest step).
        chr21_names = None
        if need_wg or need_chr:
            wg_mappers = ["shmap", "cpp"] if need_wg else ["cpp"]
            for m in wg_mappers:
                res = run_mapper(m, WG_REF, subset, a.k, a.w, theta, density, a.threads)
                if m == "cpp":
                    _, _, chr21_names = parse_paf(res["stdout"], seg=CHR_NAME)
                if need_wg:
                    row = make_row("wholegenome", ds, m, WG_REF, n_wg, res,
                                   density, theta, a.k, a.w, a.threads)
                    rows.append(row)
                    print_row(row)

        # ---- chromosome scope: the true chr21 reads mapped against chr21 only ----
        if need_chr:
            chr_subset = RESULTS / f"_chr21_{ds}.fa"
            n_chr = extract_reads(subset, chr_subset, chr21_names or set())
            if n_chr == 0:
                print(f"    -> no {CHR_NAME} reads in this subset "
                      f"(SRA reads are often position-grouped; skipping chromosome scope)")
            else:
                print(f"    -> {n_chr} reads belong to {CHR_NAME}; mapping them to {CHR_NAME} only")
                for m in ("shmap", "cpp", "py"):
                    res = run_mapper(m, CHR_REF, chr_subset, a.k, a.w, theta, density, a.threads)
                    row = make_row("chromosome", ds, m, CHR_REF, n_chr, res,
                                   density, theta, a.k, a.w, a.threads)
                    rows.append(row)
                    print_row(row)
        print()

    fields = ["scope", "dataset", "mapper", "ref", "reads_in", "mapped", "map_pct",
              "mean_mapq", "total_s", "index_s", "map_s", "index_rss_mb", "peak_rss_mb",
              "map_reads_per_s", "k", "w", "density_or_r", "theta", "threads"]
    with open(csv_path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        wr.writerows(rows)
    print(f"CSV -> {csv_path}")

    print("\n| scope | dataset | mapper | reads_in | mapped | map% | mapq | total_s | index_s | map_s | idx_MB | peak_MB |")
    print("|-------|---------|--------|---------:|-------:|-----:|-----:|--------:|--------:|------:|-------:|--------:|")
    for r in rows:
        print(f"| {r['scope']} | {r['dataset']} | {r['mapper']} | {r['reads_in']} | "
              f"{r['mapped']} | {r['map_pct']} | {r['mean_mapq']} | {r['total_s']} | "
              f"{r['index_s']} | {r['map_s']} | {r['index_rss_mb']} | {r['peak_rss_mb']} |")


if __name__ == "__main__":
    main()
