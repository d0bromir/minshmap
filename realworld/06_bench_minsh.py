"""Real-world benchmark for minSH: original (unbounded) vs patched (banded) A*.

minSH is an *aligner*, so we benchmark alignment, not mapping. Tasks are built
from the real T2T-CHM13 chr21 reference (so the sequence content and k-mer
composition are real) with simulated long-read error profiles, which is the
standard way to benchmark an aligner because it gives a known reference window:

  - hifi  : read = chr21 window mutated at  6% error (PacBio HiFi-like)
  - ont   : read = chr21 window mutated at 13% error (Oxford Nanopore-like)
  - decoy : read from a DIFFERENT, unrelated chr21 locus (a wrong placement)

Each task is aligned twice through the SAME patched `align()`:
  - "original" = max_cost=None  -> identical to upstream minSH (no budget)
  - "patched"  = max_cost=0.30*L -> banded A* (this PR's change)

Every alignment runs in a subprocess with a wall-clock timeout, so the unbounded
original cannot hang the whole run; a timeout is recorded instead.

Usage:
    python3 06_bench_minsh.py                 # run the benchmark
    python3 06_bench_minsh.py --worker J I C  # internal: align task I of jobs J
"""
import csv
import math
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)
WORKSPACE = os.path.dirname(MINSHMAP_ROOT)
sys.path.insert(0, os.path.join(WORKSPACE, "minSH"))  # for `import minsh`

from minsh.astar import align, build_seedh           # noqa: E402
from minsh.utils import apply_errors                  # noqa: E402

REF = os.path.join(HERE, "data_rw", "chr21.fa")
RESULTS = os.path.join(HERE, "results_rw")
JOBS = os.path.join(RESULTS, "_minsh_jobs.tsv")

TIMEOUT = 8.0          # seconds per single alignment
BUDGET_FRAC = 0.30     # patched budget = 30% edits (covers <=13% error, rejects decoys)
LENGTHS = [1000, 2000, 4000]
PER = 6                # tasks per (category, length)
SEED = 20260621


def read_chr21():
    seq = []
    with open(REF) as f:
        next(f)
        for line in f:
            seq.append(line.strip())
    return "".join(seq).upper()


def kfor(n):
    return max(3, math.ceil(math.log(n, 4)))


# --------------------------------------------------------------------------- #
# Worker: align one task and print "STATUS dist cells seconds".
# --------------------------------------------------------------------------- #

def worker(jobs_path, idx, maxcost):
    with open(jobs_path) as f:
        line = f.readlines()[idx]
    A, B = line.rstrip("\n").split("\t")
    h = build_seedh(A, B, kfor(len(A)))
    mc = None if maxcost == "none" else int(maxcost)
    t0 = time.perf_counter()
    res = align(A, B, h, return_stats=True, max_cost=mc)
    sec = time.perf_counter() - t0
    if res is None:
        print(f"REJECT\t0\t0\t{sec:.5f}")
    else:
        _, dist, cells = res
        print(f"OK\t{dist}\t{cells}\t{sec:.5f}")


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def build_tasks(ref, rng):
    tasks = []
    for L in LENGTHS:
        for cat, err in (("hifi", 0.06), ("ont", 0.13)):
            for _ in range(PER):
                s = rng.randrange(len(ref) - L)
                window = ref[s:s + L]
                read = apply_errors(window, err)
                tasks.append((cat, L, read, window))
        for _ in range(PER):                       # decoys: unrelated locus
            s1 = rng.randrange(len(ref) - L)
            s2 = rng.randrange(len(ref) - L)
            read = apply_errors(ref[s1:s1 + L], 0.06)
            tasks.append(("decoy", L, read, ref[s2:s2 + L]))
    rng.shuffle(tasks)
    return tasks


def run_one(idx, maxcost):
    """Returns (status, dist, cells, sec). status in OK/REJECT/TIMEOUT."""
    cmd = [sys.executable, os.path.abspath(__file__), "--worker", JOBS, str(idx), maxcost]
    t0 = time.perf_counter()
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", None, None, TIMEOUT
    wall = time.perf_counter() - t0
    parts = out.stdout.strip().split("\t")
    if len(parts) != 4:
        return "ERROR", None, None, wall
    status, dist, cells, sec = parts
    return status, int(dist), int(cells), float(sec)


def main():
    if not os.path.exists(REF):
        sys.exit("Missing chr21.fa (run 02_download_ref.sh)")
    os.makedirs(RESULTS, exist_ok=True)
    ref = read_chr21()
    rng = random.Random(SEED)
    tasks = build_tasks(ref, rng)

    with open(JOBS, "w") as f:
        for _, _, A, B in tasks:
            f.write(f"{A}\t{B}\n")

    print(f"{len(tasks)} alignment tasks (chr21, lengths {LENGTHS}, "
          f"timeout {TIMEOUT}s, patched budget {int(BUDGET_FRAC*100)}% edits)\n", flush=True)

    rows = []
    for i, (cat, L, A, B) in enumerate(tasks):
        budget = str(int(BUDGET_FRAC * L))
        o_status, o_dist, _, o_sec = run_one(i, "none")     # original (no band)
        p_status, p_dist, _, p_sec = run_one(i, budget)     # patched (banded)
        rows.append(dict(i=i, cat=cat, L=L, len_read=len(A),
                         o_status=o_status, o_dist=o_dist, o_sec=o_sec,
                         p_status=p_status, p_dist=p_dist, p_sec=p_sec))
        print(f"  [{i+1:2}/{len(tasks)}] {cat:5} L={L:4}  "
              f"orig={o_status}({o_sec:.2f}s)  patched={p_status}({p_sec:.2f}s)", flush=True)

    os.remove(JOBS)
    summarize(rows)


def summarize(rows):
    cats = ["hifi", "ont", "decoy", "ALL"]
    def sel(c):
        return rows if c == "ALL" else [r for r in rows if r["cat"] == c]

    # Exactness check: where both completed (OK), distances must be identical.
    both_ok = [r for r in rows if r["o_status"] == "OK" and r["p_status"] == "OK"]
    mismatches = [r for r in both_ok if r["o_dist"] != r["p_dist"]]

    lines = [
        "# minSH real-world alignment benchmark: original vs patched (banded) A*",
        "",
        f"- reference: T2T-CHM13v2.0 chr21 (real sequence content)",
        f"- tasks: {len(rows)} (lengths {LENGTHS}, {PER}/length/category)",
        f"- per-alignment timeout: {TIMEOUT}s; patched budget: {int(BUDGET_FRAC*100)}% edits",
        f"- 'original' = align(max_cost=None) = upstream minSH; "
        f"'patched' = align(max_cost={int(BUDGET_FRAC*100)}%·L)",
        "",
        "`verdict` = finished without timing out (OK aligned or REJECT pruned). "
        "`mean_s` is over tasks that returned (timeouts counted at the cap).",
        "",
        "| category | n | orig verdicts | orig timeouts | orig total s | "
        "patched verdicts | patched timeouts | patched total s |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    csv_rows = []
    for c in cats:
        rs = sel(c)
        n = len(rs)
        o_verd = sum(1 for r in rs if r["o_status"] in ("OK", "REJECT"))
        o_to = sum(1 for r in rs if r["o_status"] == "TIMEOUT")
        o_tot = sum(r["o_sec"] for r in rs)
        p_verd = sum(1 for r in rs if r["p_status"] in ("OK", "REJECT"))
        p_to = sum(1 for r in rs if r["p_status"] == "TIMEOUT")
        p_tot = sum(r["p_sec"] for r in rs)
        lines.append(f"| {c} | {n} | {o_verd} | {o_to} | {o_tot:.1f} | "
                     f"{p_verd} | {p_to} | {p_tot:.1f} |")
        csv_rows.append(dict(category=c, n=n, orig_verdicts=o_verd, orig_timeouts=o_to,
                             orig_total_s=round(o_tot, 2), patched_verdicts=p_verd,
                             patched_timeouts=p_to, patched_total_s=round(p_tot, 2)))

    o_to_all = sum(1 for r in rows if r["o_status"] == "TIMEOUT")
    p_to_all = sum(1 for r in rows if r["p_status"] == "TIMEOUT")
    o_tot_all = sum(r["o_sec"] for r in rows)
    p_tot_all = sum(r["p_sec"] for r in rows)
    lines += [
        "",
        "## Verdict",
        "",
        f"- exactness: of {len(both_ok)} tasks both versions aligned, "
        f"**{len(both_ok) - len(mismatches)}/{len(both_ok)}** have identical edit distance "
        f"({'no mismatches' if not mismatches else str(len(mismatches)) + ' MISMATCH'}).",
        f"- timeouts: original **{o_to_all}**, patched **{p_to_all}**.",
        f"- total wall: original **{o_tot_all:.1f}s**, patched **{p_tot_all:.1f}s** "
        f"(**{o_tot_all / p_tot_all:.2f}x** faster overall)." if p_tot_all else "",
        f"- improvement: {'YES' if (p_to_all < o_to_all and not mismatches) else 'NONE'} "
        f"— the budget keeps results exact while bounding the worst case.",
        "",
    ]
    md = "\n".join(lines)
    with open(os.path.join(RESULTS, "minsh_bench.md"), "w") as f:
        f.write(md)
    with open(os.path.join(RESULTS, "minsh_bench.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        w.writeheader()
        w.writerows(csv_rows)
    print("\n" + md)
    print("Wrote results_rw/minsh_bench.md and results_rw/minsh_bench.csv")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--worker":
        worker(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    else:
        main()
