"""Cross-language benchmark: minSHmap in Python vs. C++ (same algorithm, same data).

For each sketcher (naive / poly / nthash) it measures, in both languages, the
reference indexing throughput, the read-mapping throughput, and the mapping
accuracy against the truth encoded in the read headers. The Python numbers reuse
`bench.bench_one`; the C++ numbers come from `minshmap.exe --report` (a TSV line).

The C++ binary is (re)built automatically with g++ when it is missing or older
than the source. Results are written to results/bench_langs.json and .md.

Run from the minshmap/ folder (after data/ exists):

    python scripts/generate.py        # once, to (re)create data/
    python scripts/bench_langs.py
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import minshmap as M  # noqa: E402
import bench  # reuse Python measurement + FASTA-with-truth reader  # noqa: E402

DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SRC = os.path.join(ROOT, "minshmap.cpp")
EXE = os.path.join(ROOT, "minshmap.exe" if os.name == "nt" else "minshmap")

# Same parameters as the Python-only benchmark, so the two are comparable.
PARAMS = bench.PARAMS
HASHES = list(M.SKETCHERS.keys())


def ensure_binary():
    """Build minshmap.exe with g++ if it is missing or stale. Returns path or None."""
    gpp = shutil.which("g++")
    if gpp is None:
        if os.path.exists(EXE):
            return EXE
        print("g++ not found and no prebuilt binary; skipping C++ benchmark.")
        return None
    if os.path.exists(EXE) and os.path.getmtime(EXE) >= os.path.getmtime(SRC):
        return EXE
    print("building minshmap.cpp ...", flush=True)
    cmd = [gpp, "-O3", "-std=c++17", "-march=native", "-o", EXE, SRC]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stderr)
        return None
    return EXE


def bench_cpp(exe, hash_name, ref_path, reads_path, params):
    """Run the C++ binary in --report mode and parse its TSV line."""
    cmd = [exe, ref_path, reads_path, "--hash", hash_name,
           "-k", str(params["k"]), "-r", str(params["hfrac"]),
           "-t", str(params["theta"]), "-d", str(params["min_diff"]), "--report"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    # hash index_sec ref_bp_per_s map_sec reads_per_s mapped_frac accuracy
    f = out.split("\t")
    return {
        "sketcher": f[0],
        "index_sec": float(f[1]),
        "ref_kmers_per_sec": int(f[2]),
        "map_sec": float(f[3]),
        "reads_per_sec": int(f[4]),
        "mapped_frac": float(f[5]),
        "accuracy": float(f[6]),
    }


def to_markdown(meta, rows):
    cols = [
        ("sketcher", "sketcher"), ("lang", "lang"),
        ("index_sec", "index s"), ("ref_kmers_per_sec", "ref bp/s"),
        ("map_sec", "map s"), ("reads_per_sec", "reads/s"),
        ("mapped_frac", "mapped"), ("accuracy", "accuracy"),
    ]
    lines = [
        "# minSHmap cross-language benchmark (Python vs. C++)",
        "",
        f"- generated: {meta['timestamp']}",
        f"- reference: {meta['ref_segments']} segment(s), {meta['ref_bp']:,} bp",
        f"- reads: {meta['n_reads']} x ~{meta['read_len']} bp",
        f"- params: {meta['params']}",
        f"- python: {meta['python']}",
        f"- compiler: {meta['compiler']}",
        "",
        "| " + " | ".join(h for _, h in cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(str(r.get(k, "")) for k, _ in cols) + " |")

    # Speedup summary (C++ map throughput / Python map throughput).
    by = {(r["sketcher"], r["lang"]): r for r in rows}
    speed = []
    for h in HASHES:
        p, c = by.get((h, "python")), by.get((h, "cpp"))
        if p and c and p["reads_per_sec"]:
            speed.append(f"  - {h}: {c['reads_per_sec'] / p['reads_per_sec']:.1f}x")
    if speed:
        lines += ["", "C++ mapping speedup over Python (reads/s):", *speed]
    lines.append("")
    return "\n".join(lines)


def main():
    ref_path = os.path.join(DATA, "ref.fa")
    reads_path = os.path.join(DATA, "reads.fa")
    if not (os.path.exists(ref_path) and os.path.exists(reads_path)):
        sys.exit("Missing data/. Run: python scripts/generate.py")

    refs = list(M.read_fasta(ref_path))
    reads = list(bench.read_fasta_full(reads_path))
    os.makedirs(RESULTS, exist_ok=True)

    exe = ensure_binary()
    rows = []
    for name in HASHES:
        print(f"python  {name} ...", flush=True)
        py = bench.bench_one(name, M.SKETCHERS[name], refs, reads, PARAMS)
        py["lang"] = "python"
        rows.append(py)
        if exe:
            print(f"cpp     {name} ...", flush=True)
            cpp = bench_cpp(exe, name, ref_path, reads_path, PARAMS)
            cpp["lang"] = "cpp"
            rows.append(cpp)

    compiler = "n/a"
    if exe and shutil.which("g++"):
        compiler = subprocess.run(["g++", "--version"], capture_output=True, text=True).stdout.splitlines()[0]

    meta = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": PARAMS,
        "ref_segments": len(refs),
        "ref_bp": sum(len(s) for _, s in refs),
        "n_reads": len(reads),
        "read_len": len(reads[0][1]) if reads else 0,
        "python": sys.version.split()[0],
        "compiler": compiler,
    }
    payload = {"meta": meta, "results": rows}

    with open(os.path.join(RESULTS, "bench_langs.json"), "w") as f:
        json.dump(payload, f, indent=2)
    md = to_markdown(meta, rows)
    with open(os.path.join(RESULTS, "bench_langs.md"), "w") as f:
        f.write(md)

    print("\n" + md)
    print("Wrote results/bench_langs.json and results/bench_langs.md")


if __name__ == "__main__":
    main()
