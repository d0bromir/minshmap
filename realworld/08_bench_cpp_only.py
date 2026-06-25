#!/usr/bin/env python3
"""Benchmark ONLY the minSHmap C++ binary, pre-filling the other tools.

Per the current policy we no longer re-run shmap / minSH (they are expensive and
their numbers are stable): this script freshly times only `minshmap_linux` on
chr21 and copies the `shmap` and `minsh` rows verbatim from the existing baseline
CSV (results_rw/realworld.csv). Each output row carries a `source` column:
`measured` (this run) or `baseline` (carried over).

minSHmap now has a single sketcher (ntHash) and no max_matches, matching the
pedagogical minshmap.py, so there is one `minshmap-cpp` row per dataset. Same
read subsets and core parameters as the baseline: the first MAX_READS reads of
each data_rw/{hifi,ont,clr}.fa.

Run inside WSL:
  python3 08_bench_cpp_only.py
Writes results_rw/realworld_cpp_<YYYYMMDD-HHMMSS>.csv.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)                 # minshmap/
DATA = os.path.join(HERE, "data_rw")
RESULTS = os.path.join(HERE, "results_rw")
CPP_BIN = os.path.join(MINSHMAP_ROOT, "minshmap_linux")
REF = os.path.join(DATA, "chr21.fa")
BASELINE_CSV = os.path.join(RESULTS, "realworld.csv")

DATASETS = [("hifi", "PacBio HiFi"), ("ont", "Oxford Nanopore"), ("clr", "PacBio CLR")]
MAX_READS = 2000          # same subset size as the baseline (realworld.csv)

# Core mapping parameters (identical to minshmap.py / minshmap.cpp on the CLI).
# (w,k)-minimizer containment runs lower than the old FracMinHash (one error changes a
# whole window's minimum), so the acceptance threshold theta must track the read error
# rate. Like minimap2's map-hifi/map-ont/map-pb presets we use a PER-DATASET theta:
# higher for low-error HiFi, lower for noisy ONT/CLR where few k-mers survive intact
# ((1-e)^k ~ 0.15 ONT, ~0.09 CLR at k=15). Values picked to stay AT/UNDER shmap's truth
# (hifi 17 / ont 32 / clr 2) -> no false-positive explosion (theta<=0.05 blows CLR to 48+).
PARAMS = dict(k=15, window=10)
THETA = {"hifi": 0.20, "ont": 0.15, "clr": 0.18}   # hifi 16, ont 28, clr 2 (<= shmap truth)
PREFILL_TOOLS = ("shmap", "minsh")   # rows copied from the baseline, not re-run


def read_reads(path, limit):
    """Yield up to `limit` (name, seq) records from a FASTA file."""
    name, chunks, n = None, [], 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                    n += 1
                    if n >= limit:
                        return
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None and n < limit:
        yield name, "".join(chunks)


def parse_paf_mapped(text):
    """Count distinct reads with at least one PAF placement line."""
    mapped = set()
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 9:
            mapped.add(f[0])
    return len(mapped)


def load_prefill(path):
    """Read the baseline CSV -> list of shmap / minsh rows (carried over verbatim)."""
    if not os.path.exists(path):
        sys.exit(f"Missing baseline CSV: {path}")
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["tool"] in PREFILL_TOOLS:
                rows.append(row)
    return rows


def main():
    if not os.path.exists(CPP_BIN):
        sys.exit(f"Missing C++ binary: {CPP_BIN} (build it in WSL first)")
    if not os.path.exists(REF):
        sys.exit("Missing chr21.fa. Run 02_download_ref.sh")
    os.makedirs(RESULTS, exist_ok=True)

    prefill_rows = load_prefill(BASELINE_CSV)

    rows = []     # assembled output rows (measured + baseline)
    for label, platform in DATASETS:
        reads_path = os.path.join(DATA, f"{label}.fa")
        if not os.path.exists(reads_path):
            print(f"skip {label}: {reads_path} missing", flush=True)
            continue

        # Reproduce the baseline read subset exactly (first MAX_READS reads).
        reads = list(read_reads(reads_path, MAX_READS))
        n = len(reads)
        cap_path = os.path.join(RESULTS, f"_cap_{label}.fa")
        with open(cap_path, "w") as f:
            for nm, sq in reads:
                f.write(f">{nm}\n{sq}\n")

        kc = ["-k", str(PARAMS["k"]), "-w", str(PARAMS["window"]), "-t", str(THETA[label])]
        print(f"  minshmap-cpp mapping {label} ({n} reads, theta={THETA[label]}) ...", flush=True)
        t0 = time.perf_counter()
        out = subprocess.run([CPP_BIN, REF, cap_path] + kc, capture_output=True, text=True)
        sec = time.perf_counter() - t0
        mapped = parse_paf_mapped(out.stdout)
        rows.append({
            "dataset": label, "platform": platform, "tool": "minshmap-cpp",
            "mapped": mapped,
            "mapped_frac": round(mapped / n, 4) if n else 0,
            "wall_sec": round(sec, 3),
            "reads_per_sec": round(n / sec) if sec else "",
            "agree_with_shmap": "",   # not recomputed (shmap is not re-run)
            "identity": "",
            "source": "measured",
        })

        # Carry over the pre-filled rows for this dataset (shmap, minsh).
        for r in prefill_rows:
            if r["dataset"] == label:
                rows.append({**r, "source": "baseline"})

        os.remove(cap_path)

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_path = os.path.join(RESULTS, f"realworld_cpp_{stamp}.csv")
    cols = ["dataset", "platform", "tool", "mapped", "mapped_frac", "wall_sec",
            "reads_per_sec", "agree_with_shmap", "identity", "source"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote {out_path}")
    print(f"  measured rows: {sum(1 for r in rows if r['source'] == 'measured')}")
    print(f"  baseline rows: {sum(1 for r in rows if r['source'] == 'baseline')}")


if __name__ == "__main__":
    main()
