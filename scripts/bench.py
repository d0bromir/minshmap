"""Benchmark the three minSHmap sketchers on a generated reference + reads.

Measures, per sketcher (naive / poly / nthash):
  - sketch throughput  (reference k-mers/sec during indexing)
  - index build time
  - read mapping time
  - mapped fraction and accuracy vs. the truth in the read headers

Reproducible: fixed data (scripts/generate.py) + fixed parameters here.
Writes structured results to results/bench.json and results/bench.md.

Run from the repo root:

    python scripts/generate.py        # once, to (re)create data/
    python scripts/bench.py
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import minshmap as M  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")

# Fixed benchmark parameters (documented for reproducibility).
# k=11 at 5% error gives expected containment (1-0.05)^11 ~= 0.57 > theta=0.4,
# so well-placed reads clear the threshold and the comparison is meaningful.
PARAMS = dict(k=11, hfrac=0.1, theta=0.4, min_diff=0.02)


def parse_truth(header):
    """Extract (segm, pos, strand) from a generated read header."""
    fields = dict(tok.split("=") for tok in header.split() if "=" in tok)
    return fields.get("segm"), int(fields["pos"]), fields.get("strand")


def read_fasta_full(path):
    """Like minshmap.read_fasta but keeps the *whole* header line (truth lives
    after the read name, which the mapper's reader would otherwise drop)."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def bench_one(name, sketch, refs, reads, params):
    k, hfrac = params["k"], params["hfrac"]
    theta, min_diff = params["theta"], params["min_diff"]

    # --- indexing (includes reference sketching) ---
    ref_bp = sum(len(s) for _, s in refs)
    t0 = time.perf_counter()
    index, segments = M.build_index(refs, k, hfrac, sketch)
    t_index = time.perf_counter() - t0
    indexed_kmers = sum(len(v) for v in index.values())

    seg_id = {seg.name: i for i, seg in enumerate(segments)}

    # --- mapping ---
    mapped = correct = 0
    t0 = time.perf_counter()
    for header, seq in reads:
        best, second = M.map_read(seq, index, segments, k, hfrac, theta, min_diff, sketch)
        if best is None:
            continue
        mapped += 1
        tname, tpos, _ = parse_truth(header)
        if best.segm_id == seg_id.get(tname) and best.t_start <= tpos + len(seq) // 2 <= best.t_end:
            correct += 1
    t_map = time.perf_counter() - t0

    n = len(reads)
    return {
        "sketcher": name,
        "index_sec": round(t_index, 4),
        "ref_kmers_per_sec": round(ref_bp / t_index) if t_index else None,
        "indexed_kmers": indexed_kmers,
        "map_sec": round(t_map, 4),
        "reads_per_sec": round(n / t_map) if t_map else None,
        "mapped": mapped,
        "mapped_frac": round(mapped / n, 3),
        "correct": correct,
        "accuracy": round(correct / n, 3),
    }


def to_markdown(meta, rows):
    cols = [
        ("sketcher", "sketcher"), ("index_sec", "index s"),
        ("ref_kmers_per_sec", "ref bp/s"), ("map_sec", "map s"),
        ("reads_per_sec", "reads/s"), ("mapped_frac", "mapped"),
        ("accuracy", "accuracy"),
    ]
    lines = [
        "# minSHmap sketcher benchmark",
        "",
        f"- generated: {meta['timestamp']}",
        f"- reference: {meta['ref_segments']} segment(s), {meta['ref_bp']:,} bp",
        f"- reads: {meta['n_reads']} x ~{meta['read_len']} bp",
        f"- params: {meta['params']}",
        "",
        "| " + " | ".join(h for _, h in cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(str(r[k]) for k, _ in cols) + " |")
    lines.append("")
    return "\n".join(lines)


def main():
    ref_path = os.path.join(DATA, "ref.fa")
    reads_path = os.path.join(DATA, "reads.fa")
    if not (os.path.exists(ref_path) and os.path.exists(reads_path)):
        sys.exit("Missing data/. Run: python scripts/generate.py")

    refs = list(M.read_fasta(ref_path))
    reads = list(read_fasta_full(reads_path))
    os.makedirs(RESULTS, exist_ok=True)

    rows = []
    for name, sketch in M.SKETCHERS.items():
        print(f"benchmarking {name} ...", flush=True)
        rows.append(bench_one(name, sketch, refs, reads, PARAMS))

    meta = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": PARAMS,
        "ref_segments": len(refs),
        "ref_bp": sum(len(s) for _, s in refs),
        "n_reads": len(reads),
        "read_len": len(reads[0][1]) if reads else 0,
    }
    payload = {"meta": meta, "results": rows}

    with open(os.path.join(RESULTS, "bench.json"), "w") as f:
        json.dump(payload, f, indent=2)
    md = to_markdown(meta, rows)
    with open(os.path.join(RESULTS, "bench.md"), "w") as f:
        f.write(md)

    print("\n" + md)
    print(f"Wrote results/bench.json and results/bench.md")


if __name__ == "__main__":
    main()
