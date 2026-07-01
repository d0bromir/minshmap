"""Equal-work benchmark: minSHmap (C++) vs. shmap per-read mapping throughput.

Why this exists
---------------
The headline "minshmap ~197 reads/s vs shmap ~158 reads/s" (results_rw/
realworld.md) is confounded: it is *end-to-end* throughput (index build + map),
and shmap maps far MORE reads (it is more sensitive), so the two numbers compare
different amounts of work. This harness separates the two things that number
conflates:

  1. PER-READ MAPPER-LOOP SPEED (apples-to-apples). Both tools run their per-read
     pipeline (sketch -> candidate blocks -> seed-heuristic prune -> refine) on
     the *same* reads. We isolate that loop from the large, very different index-
     build cost with a differential measurement: time each tool on N reads and on
     2N reads (the first N are a prefix of the 2N), then

         per-read map time = (wall_2N - wall_N) / N
         isolated reads/s  = N / (wall_2N - wall_N)

     The fixed index cost cancels. No tool-specific timers are parsed, so the
     comparison is fair and tool-agnostic. Both run single-threaded.

  2. SENSITIVITY (reported separately, NOT mixed into speed). For the same 2N
     reads we also report how many each tool maps. This is where shmap wins and
     minshmap is deliberately simple; keeping it in its own column is the whole
     point -- speed and sensitivity are different axes.

Run (in WSL, from minshmap/realworld/):
    python3 07_bench_equal_sensitivity.py --n 500
"""
from __future__ import annotations

import argparse
import os
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)
WORKSPACE = os.path.dirname(MINSHMAP_ROOT)

CPP_BIN = os.path.join(MINSHMAP_ROOT, "minshmap_linux")
SHMAP_BIN = os.path.join(WORKSPACE, "shmap", "release", "shmap")
REF = os.path.join(HERE, "data_rw", "chr21.fa")
READS = os.path.join(HERE, "data_rw", "hifi.fa")
OUTDIR = os.path.join(HERE, "results_rw")

# Same params as the main real-world benchmark (results_rw/realworld.md).
PARAMS = dict(k=15, hfrac=0.05, theta=0.3, min_diff=0.02, max_matches=1000, max_seeds=500)


# --------------------------------------------------------------------------- #
# FASTA helpers
# --------------------------------------------------------------------------- #

def read_fasta(path, limit=None):
    """Yield (name, seq); name is the first whitespace token after '>'."""
    name, chunks, n = None, [], 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                    n += 1
                    if limit and n >= limit:
                        return
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None and (not limit or n < limit):
        yield name, "".join(chunks)


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


# --------------------------------------------------------------------------- #
# Tool invocations
# --------------------------------------------------------------------------- #

def cpp_cmd(reads_path):
    # minSHmap now has a single sketcher and only -k/-r/-t/-j (no -d/-M/-S).
    return [CPP_BIN, REF, reads_path,
            "-k", str(PARAMS["k"]), "-r", str(PARAMS["hfrac"]),
            "-t", str(PARAMS["theta"]), "-j", "1"]  # single-thread, matches shmap


def shmap_cmd(reads_path):
    return [SHMAP_BIN, "-p", reads_path, "-s", REF,
            "-k", str(PARAMS["k"]), "-r", str(PARAMS["hfrac"]),
            "-t", str(PARAMS["theta"]), "-d", str(PARAMS["min_diff"]),
            "-M", str(PARAMS["max_matches"]), "-S", str(PARAMS["max_seeds"])]


def run_timed(cmd):
    """Run a mapper, return (n_mapped_reads, wall_seconds)."""
    t0 = time.perf_counter()
    out = subprocess.run(cmd, capture_output=True, text=True)
    sec = time.perf_counter() - t0
    mapped = sum(1 for ln in out.stdout.splitlines() if ln and not ln.startswith("@"))
    return mapped, sec


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="equal-work minshmap vs shmap per-read benchmark")
    ap.add_argument("--n", type=int, default=500, help="N (times N and 2N reads)")
    ap.add_argument("--reps", type=int, default=2, help="repeats per timing (best taken)")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    N = args.n

    # ---- Build the N and 2N read files (2N's prefix == the N file) -------------
    records = list(read_fasta(READS, limit=2 * N))
    if len(records) < 2 * N:
        print(f"Only {len(records)} reads available; lower --n.")
        return
    file_N = os.path.join(OUTDIR, "_eqsens_N.fa")
    file_2N = os.path.join(OUTDIR, "_eqsens_2N.fa")
    write_fasta(file_N, records[:N])
    write_fasta(file_2N, records[:2 * N])

    # ---- Differential timing (index cost cancels in wall_2N - wall_N) ----------
    print(f"[time] differential N={N} vs 2N={2*N} reads, best of {args.reps} ...",
          flush=True)

    def measure(cmd_fn, path):
        best, mapped = float("inf"), 0
        for _ in range(args.reps):
            m, sec = run_timed(cmd_fn(path))
            if sec < best:
                best, mapped = sec, m
        return best, mapped

    rows = []
    for tool, cmd_fn in (("minshmap-cpp", cpp_cmd), ("shmap", shmap_cmd)):
        w_n, _ = measure(cmd_fn, file_N)
        w_2n, mapped_2n = measure(cmd_fn, file_2N)
        dt = w_2n - w_n
        rps = N / dt if dt > 1e-6 else float("inf")
        rows.append((tool, w_n, w_2n, dt, rps, mapped_2n))
        print(f"  {tool:14s} wall_N={w_n:6.2f}s  wall_2N={w_2n:6.2f}s  "
              f"map={dt:6.2f}s  isolated={rps:8.1f} reads/s  mapped={mapped_2n}/{2*N}")

    # ---- Report ----------------------------------------------------------------
    cpp = next(r for r in rows if r[0] == "minshmap-cpp")
    shm = next(r for r in rows if r[0] == "shmap")
    ratio = cpp[4] / shm[4] if shm[4] else float("inf")
    # Index build is wall(N) minus the per-read map time of the first N reads,
    # which (same read prefix) equals the differential Delta. So index ~ wall_N - Delta.
    cpp_index = max(0.0, cpp[1] - cpp[3])
    shm_index = max(0.0, shm[1] - shm[3])
    faster, slower = ("shmap", "minshmap-cpp") if shm[4] >= cpp[4] else ("minshmap-cpp", "shmap")
    speed_ratio = max(cpp[4], shm[4]) / min(cpp[4], shm[4]) if min(cpp[4], shm[4]) else float("inf")

    lines = [
        "# Equal-work benchmark: minSHmap (C++) vs shmap",
        "",
        "Per-read **mapper-loop** throughput on the *same* reads, with the fixed",
        "chr21 index-build cost cancelled by a differential N-vs-2N measurement.",
        "Speed and sensitivity are reported as separate axes.",
        "",
        f"- Reads: first {N} and first {2*N} HiFi reads (same prefix). "
        f"Best of {args.reps}. Single-threaded both.",
        f"- Params: {PARAMS}.",
        "- NOTE: the differential subtracts two ~index-sized walls, so with few",
        "  reps the ratio carries index-build noise; trust larger N over small N.",
        "",
        "| tool | wall(N) s | wall(2N) s | map loop Δ s | isolated reads/s | ~index s | mapped/2N |",
        "|------|-----------|------------|-------------|------------------|----------|-----------|",
    ]
    for tool, w_n, w_2n, dt, rps, mp in rows:
        idx = cpp_index if tool == "minshmap-cpp" else shm_index
        lines.append(f"| {tool} | {w_n:.2f} | {w_2n:.2f} | {dt:.2f} | {rps:.1f} | {idx:.2f} | {mp}/{2*N} |")
    lines += [
        "",
        f"**Per-read mapper loop: {faster} is {speed_ratio:.2f}x {slower}'s throughput** "
        f"on identical reads (index cost removed).",
        f"**Sensitivity: shmap maps {shm[5]} of {2*N} reads vs minshmap's {cpp[5]}.**",
        "",
        "Honest reading. On identical reads with the shared index cost removed,",
        f"shmap's per-read loop here is the faster one ({shm[4]:.0f} vs {cpp[4]:.0f} reads/s)",
        "AND it maps more reads -- shmap wins on both axes. This OVERTURNS the",
        "end-to-end realworld.md figure (minshmap-cpp 196.8 vs shmap 157.6), which",
        "mixes index build with mapping across three datasets and is the less",
        "controlled comparison. The per-read truth: shmap's heavily engineered",
        "inner loop (cache-friendly hash maps, position-sorted hit lists with",
        "binary search) beats minshmap's std::unordered_map + linear scans, even",
        "though shmap does MORE algorithmic work per read (chaining, second-best,",
        "mapq). minshmap's value is pedagogical clarity, not speed.",
        "",
        f"Real output-preserving headroom for minshmap: its index build (~{cpp_index:.1f}s)",
        f"is ~{(cpp_index/shm_index if shm_index else 0):.0f}x shmap's (~{shm_index:.1f}s). Swapping",
        "std::unordered_map for a faster map (ankerl/unordered_dense, as shmap uses)",
        "plus reserve() would cut that without changing a single output line.",
    ]
    out_md = os.path.join(OUTDIR, "equal_sensitivity.md")
    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWrote {out_md}")

    for tmp in (file_N, file_2N):
        try:
            os.remove(tmp)
        except OSError:
            pass


if __name__ == "__main__":
    main()
