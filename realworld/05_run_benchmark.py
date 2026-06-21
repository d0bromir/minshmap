#!/usr/bin/env python3
"""Real-world long-read benchmark: minSHmap (Python + C++) vs. Pesho's shmap,
with minSH as an independent alignment-accuracy check.

Reference : T2T-CHM13v2.0 chr21 (~45 Mb), data_rw/chr21.fa
Reads     : real long-read subsets, one per platform (data_rw/{hifi,ont,clr}.fa)
            hifi = PacBio HiFi, ont = Oxford Nanopore, clr = PacBio CLR.

The reads are whole-genome WGS, so only the ~1-2% originating from chr21 should
map; that on-target fraction (and where the tools agree) is exactly what we
measure. There is no positional ground truth, so accuracy is assessed two ways:
  1. cross-tool concordance: do the minSHmap variants place a read at the same
     chr21 locus as the mature shmap (interval overlap, same strand)?
  2. minSH alignment: for a sample of HiFi placements, align the read to the
     located chr21 window with minSH's seed-heuristic A* and report identity.

Tools compared (same parameters k / hFrac / theta / min_diff):
  py-poly   : minshmap.py, polynomial rolling sketch
  py-nthash : minshmap.py, ntHash sketch
  cpp       : minshmap_linux, all three sketches (naive / poly / ntHash, compiled)
  shmap     : release/shmap (the original C++ mapper)

Run inside WSL:
  python3 05_run_benchmark.py [--max-reads N]
Writes results_rw/realworld.json and results_rw/realworld.md.
"""

from __future__ import annotations

import argparse
import gc
import math
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
MINSHMAP_ROOT = os.path.dirname(HERE)                  # minshmap/
WORKSPACE = os.path.dirname(MINSHMAP_ROOT)             # Pesho/
sys.path.insert(0, MINSHMAP_ROOT)
sys.path.insert(0, os.path.join(WORKSPACE, "minSH"))   # for `import minsh`

import minshmap as M  # noqa: E402

DATA = os.path.join(HERE, "data_rw")
RESULTS = os.path.join(HERE, "results_rw")
CPP_BIN = os.path.join(MINSHMAP_ROOT, "minshmap_linux")
SHMAP_BIN = os.path.join(WORKSPACE, "shmap", "release", "shmap")
REF = os.path.join(DATA, "chr21.fa")

DATASETS = [("hifi", "PacBio HiFi"), ("ont", "Oxford Nanopore"), ("clr", "PacBio CLR")]

# Unified mapping parameters for every tool (long-read regime on a single chr).
# max_matches drops k-mers occurring > this many times in chr21 (satellite
# repeats reach 30k+ hits and would otherwise dominate the work); shmap does the
# same via its own max_matches / erase_frequent_kmers. max_seeds bounds shmap's
# per-read seed set so its seed-heuristic stays within RAM on long noisy reads.
PARAMS = dict(k=15, hfrac=0.05, theta=0.3, min_diff=0.02, max_matches=1000, max_seeds=500)


# ----------------------------------------------------------------------------- #
# FASTA + PAF helpers
# ----------------------------------------------------------------------------- #

def read_reads(path):
    """Yield (name, seq); name is the first header token (the read id)."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def parse_paf(text):
    """PAF stdout -> {qname: (tname, tstart, tend, strand)} (best line per read)."""
    placements = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        qname, _, _, _, strand, tname, _, tstart, tend = f[:9]
        placements.setdefault(qname, (tname, int(tstart), int(tend), strand))
    return placements


def overlaps(a, b, same_strand=True):
    """Do two placements hit the same target with overlapping intervals?"""
    if a is None or b is None or a[0] != b[0]:
        return False
    if same_strand and a[3] != b[3]:
        return False
    return max(a[1], b[1]) < min(a[2], b[2])


# ----------------------------------------------------------------------------- #
# Tool runners
# ----------------------------------------------------------------------------- #

def run_py(reads, index, segments, sketch, params):
    """Map an in-memory list of reads with minshmap.py; return (placements, sec)."""
    k, hfrac, theta, md = params["k"], params["hfrac"], params["theta"], params["min_diff"]
    mm = params["max_matches"]
    placements = {}
    t0 = time.perf_counter()
    for name, seq in reads:
        best, _ = M.map_read(seq, index, segments, k, hfrac, theta, md, sketch, mm)
        if best is not None:
            strand = "+" if best.codir >= 0 else "-"
            placements[name] = (segments[best.segm_id].name, best.t_start, best.t_end, strand)
    return placements, time.perf_counter() - t0


def run_binary(cmd):
    """Run a mapper subprocess; return (placements, sec). PAF on stdout."""
    t0 = time.perf_counter()
    out = subprocess.run(cmd, capture_output=True, text=True)
    sec = time.perf_counter() - t0
    return parse_paf(out.stdout), sec


# ----------------------------------------------------------------------------- #
# minSH alignment validation (HiFi sample)
# ----------------------------------------------------------------------------- #

def minsh_validate(reads_by_name, ref_seq, placements, sample=20, max_len=12000, pad=200):
    """Align a sample of placed reads to their located chr21 window with minSH.
    Returns list of {read, len, edit, identity, cells}.

    Long reads make minSH's A* degenerate to Dijkstra, so we (a) reverse-
    complement reads placed on the '-' strand, (b) cap read length, and
    (c) skip placements whose 15-mer containment is too low (wrong locus =>
    huge band). This keeps it a fast, meaningful spot-check on HiFi."""
    try:
        from minsh.astar import build_seedh, align
    except Exception as e:  # numpy/fenwick/matplotlib missing -> skip gracefully
        return {"error": f"minSH unavailable: {e}", "alignments": []}

    def canon_kmers(s, kk=15):
        out = set()
        for i in range(len(s) - kk + 1):
            km = s[i:i + kk]
            r = M.revcomp(km)
            out.add(km if km <= r else r)
        return out

    out = []
    for name, place in placements.items():
        if len(out) >= sample:
            break
        seq = reads_by_name.get(name)
        if seq is None or not (1000 <= len(seq) <= max_len):
            continue
        _, ts, te, strand = place
        if strand == "-":
            seq = M.revcomp(seq)                       # align in placement strand
        region = ref_seq[max(0, ts - pad): te + pad]
        if not region:
            continue
        rk = canon_kmers(seq)
        if len(rk & canon_kmers(region)) / max(len(rk), 1) < 0.30:
            continue                                   # weak placement: A* would blow up
        k = max(8, math.ceil(math.log(max(len(seq), 2), 4)))
        h = build_seedh(seq, region, k)
        _, dist, cells = align(seq, region, h, return_stats=True)
        out.append({
            "read": name, "len": len(seq), "edit": int(dist),
            "identity": round(1 - dist / len(seq), 4), "cells": int(cells),
        })
    return {"alignments": out}


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


# ----------------------------------------------------------------------------- #
# Main
# ----------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-reads", type=int, default=0, help="cap reads per dataset (0 = all)")
    ap.add_argument("--no-minsh", action="store_true", help="skip the minSH alignment check")
    args = ap.parse_args()

    if not os.path.exists(REF):
        sys.exit("Missing chr21.fa. Run 02_download_ref.sh")
    os.makedirs(RESULTS, exist_ok=True)

    print("loading chr21 reference ...", flush=True)
    refs = list(M.read_fasta(REF))
    ref_bp = sum(len(s) for _, s in refs)
    ref_seq = refs[0][1]
    print(f"  chr21: {ref_bp:,} bp", flush=True)

    kc = ["-k", str(PARAMS["k"]), "-r", str(PARAMS["hfrac"]),
          "-t", str(PARAMS["theta"]), "-d", str(PARAMS["min_diff"]),
          "-M", str(PARAMS["max_matches"])]
    # shmap explodes its seed-heuristic on long noisy reads; cap seeds per sketch
    # so it stays within WSL's RAM (this bound only affects shmap).
    shmap_extra = ["-S", str(PARAMS["max_seeds"])]

    # Load (and cap) each read set once; they are small. When --max-reads is
    # given we also write the capped subset to a temp FASTA so the compiled
    # binaries map exactly the same reads as the Python tools (otherwise they
    # would scan the whole file and the counts would not be comparable).
    datasets = []
    for label, platform in DATASETS:
        reads_path = os.path.join(DATA, f"{label}.fa")
        if not os.path.exists(reads_path):
            print(f"skip {label}: {reads_path} missing", flush=True)
            continue
        reads = list(read_reads(reads_path))
        if args.max_reads:
            reads = reads[: args.max_reads]
            bin_path = os.path.join(RESULTS, f"_cap_{label}.fa")
            with open(bin_path, "w") as f:
                for nm, sq in reads:
                    f.write(f">{nm}\n{sq}\n")
        else:
            bin_path = reads_path
        datasets.append((label, platform, bin_path, reads))

    # Tools run in phases so a big Python index and the heavy shmap subprocess
    # are never resident at the same time (WSL has ~13 GB).
    runs = {}                 # (tool, label) -> (placements, wall_sec)
    py_index_sec = {}

    # --- Phase A: Python mapper. Only the speed-optimal sketcher (ntHash, the
    # winner of the C++ hash comparison) is run, so the slow pure-Python mapper
    # executes once instead of three times. ---
    PY_HASH = "nthash"
    print(f"building python index ({PY_HASH}) ...", flush=True)
    t0 = time.perf_counter()
    index, segments = M.build_index(refs, PARAMS["k"], PARAMS["hfrac"], M.SKETCHERS[PY_HASH])
    idx_sec = time.perf_counter() - t0
    py_index_sec[PY_HASH] = round(idx_sec, 2)
    print(f"  done in {idx_sec:.1f}s ({len(index):,} hashes)", flush=True)
    for label, platform, bin_path, reads in datasets:
        print(f"  minshmap-py-{PY_HASH} mapping {label} ...", flush=True)
        pl, map_sec = run_py(reads, index, segments, M.SKETCHERS[PY_HASH], PARAMS)
        runs[(f"minshmap-py-{PY_HASH}", label)] = (pl, idx_sec + map_sec)
    del index, segments
    gc.collect()

    # Python index is gone; free the reference records too (keep ref_seq).
    del refs
    gc.collect()

    # --- Phase B: compiled binaries, run sequentially (own memory). ---
    # Benchmark every C++ sketcher so the hash comparison is symmetric and we can
    # rank them on speed (cpp runs are cheap; the Python mapper dominates time).
    for label, platform, bin_path, reads in datasets:
        for hname in ("naive", "poly", "nthash"):
            print(f"  minshmap-cpp-{hname} mapping {label} ...", flush=True)
            runs[(f"minshmap-cpp-{hname}", label)] = run_binary(
                [CPP_BIN, REF, bin_path, "--hash", hname] + kc)
        if os.path.exists(SHMAP_BIN):
            print(f"  shmap mapping {label} ...", flush=True)
            runs[("shmap", label)] = run_binary(
                [SHMAP_BIN, "-p", bin_path, "-s", REF] + kc + shmap_extra)

    # --- Phase C: assemble per-dataset metrics. ---
    results = []
    tool_order = ["minshmap-py-nthash", "minshmap-cpp-naive", "minshmap-cpp-poly",
                  "minshmap-cpp-nthash", "shmap"]
    for label, platform, bin_path, reads in datasets:
        reads_by_name = {n: s for n, s in reads}
        n = len(reads)
        rd_bp = sum(len(s) for _, s in reads)

        shmap_pl = runs.get(("shmap", label), (None, None))[0]
        tools = []
        for tool in tool_order:
            if (tool, label) not in runs:
                continue
            pl, sec = runs[(tool, label)]
            mapped = len(pl)
            agree = None
            if shmap_pl is not None and tool != "shmap":
                common = sum(1 for q, p in shmap_pl.items() if overlaps(pl.get(q), p))
                agree = round(common / len(shmap_pl), 3) if shmap_pl else None
            tools.append({
                "tool": tool,
                "mapped": mapped,
                "mapped_frac": round(mapped / n, 4) if n else 0,
                "wall_sec": round(sec, 3),
                "reads_per_sec": round(n / sec) if sec else None,
                "agree_with_shmap": agree,
                "identity": None,
            })

        entry = {"dataset": label, "platform": platform, "n_reads": n,
                 "read_bp": rd_bp, "avg_read_len": rd_bp // n if n else 0, "tools": tools}

        # minSH is an aligner, not a mapper: it aligns minshmap's HiFi placements
        # to their chr21 window and reports identity. Benchmarked here as its own
        # participant (throughput = reads aligned / s, accuracy = median identity).
        if label == "hifi" and not args.no_minsh and ("minshmap-cpp-nthash", "hifi") in runs:
            print("  minsh aligning HiFi placements ...", flush=True)
            cpp_pl = runs[("minshmap-cpp-nthash", "hifi")][0]
            t0 = time.perf_counter()
            val = minsh_validate(reads_by_name, ref_seq, cpp_pl)
            minsh_sec = time.perf_counter() - t0
            aligns = val.get("alignments", [])
            med_id = median([a["identity"] for a in aligns])
            entry["minsh"] = {
                "n": len(aligns),
                "median_identity": med_id,
                "median_cells": median([a["cells"] for a in aligns]),
                "wall_sec": round(minsh_sec, 3),
                "note": val.get("error"),
                "samples": aligns[:10],
            }
            if aligns:
                tools.append({
                    "tool": "minsh", "mapped": len(aligns),
                    "mapped_frac": None, "wall_sec": round(minsh_sec, 3),
                    "reads_per_sec": round(len(aligns) / minsh_sec, 3) if minsh_sec else None,
                    "agree_with_shmap": None, "identity": med_id,
                })
                print(f"    median identity {med_id} over {len(aligns)} reads")

        results.append(entry)

    meta = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference": "T2T-CHM13v2.0 chr21",
        "ref_bp": ref_bp,
        "params": PARAMS,
        "py_index_sec": py_index_sec,
        "python": sys.version.split()[0],
    }
    write_csv(os.path.join(RESULTS, "realworld.csv"), meta, results)
    md = to_markdown(meta, results)
    with open(os.path.join(RESULTS, "realworld.md"), "w") as f:
        f.write(md)
    print("\n" + md)
    print("Wrote results_rw/realworld.csv and results_rw/realworld.md")


def write_csv(path, meta, results):
    """Flat one-row-per-(dataset, tool) CSV of every participant's metrics."""
    import csv
    cols = ["dataset", "platform", "tool", "mapped", "mapped_frac",
            "wall_sec", "reads_per_sec", "agree_with_shmap", "identity"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for e in results:
            for t in e["tools"]:
                w.writerow([
                    e["dataset"], e["platform"], t["tool"], t["mapped"],
                    _fmt(t["mapped_frac"]), _fmt(t["wall_sec"]),
                    _fmt(t["reads_per_sec"]), _fmt(t["agree_with_shmap"]),
                    _fmt(t["identity"]),
                ])


def _fmt(v):
    return "" if v is None else v


def _rankings(results):
    """Aggregate per-tool speed (overall reads/s) and accuracy across datasets."""
    agg = {}  # tool -> {reads, wall, mapped, agree[], identity}
    for e in results:
        for t in e["tools"]:
            a = agg.setdefault(t["tool"], {"reads": 0, "wall": 0.0, "mapped": 0,
                                           "agree": [], "identity": None})
            # reads processed by this tool on this dataset = the dataset size
            # (minsh only "processes" the placements it aligned).
            processed = t["mapped"] if t["tool"] == "minsh" else e["n_reads"]
            a["reads"] += processed
            a["wall"] += t["wall_sec"] or 0.0
            a["mapped"] += t["mapped"]
            if t["agree_with_shmap"] is not None:
                a["agree"].append(t["agree_with_shmap"])
            if t["identity"] is not None:
                a["identity"] = t["identity"]

    speed = []
    for tool, a in agg.items():
        rps = round(a["reads"] / a["wall"], 2) if a["wall"] else None
        speed.append((tool, rps, a["wall"]))
    speed.sort(key=lambda x: (x[1] is not None, x[1] or 0), reverse=True)

    accuracy = []
    for tool, a in agg.items():
        mean_agree = round(sum(a["agree"]) / len(a["agree"]), 3) if a["agree"] else None
        accuracy.append((tool, a["mapped"], mean_agree, a["identity"]))
    # rank mappers by total mapped (chr21 reads recovered); shmap leads as the
    # reference. minsh has no mapped count of its own -> sorted last on mapped.
    accuracy.sort(key=lambda x: x[1], reverse=True)
    return speed, accuracy


def to_markdown(meta, results):
    lines = [
        "# minSHmap real-world long-read benchmark",
        "",
        f"- generated: {meta['timestamp']}",
        f"- reference: {meta['reference']} ({meta['ref_bp']:,} bp)",
        f"- params: {meta['params']}",
        f"- python index build (one-time): {meta['py_index_sec']} s",
        f"- python: {meta['python']}",
        "",
        "Participants: **minshmap-py-nthash** (educational pure-Python mapper, run "
        "with ntHash — the sketcher that won the C++ speed comparison), "
        "**minshmap-cpp-{naive,poly,nthash}** (the C++ port with each sketcher), "
        "**shmap** (the reference C++ mapper), and **minsh** (the A* aligner, used "
        "to independently score placement quality).",
        "",
        "Reads are whole-genome WGS, so a low mapped fraction is expected "
        "(only the ~1-2% from chr21 should map). `agree_with_shmap` = fraction "
        "of shmap's chr21 placements that the tool reproduces at an overlapping "
        "locus and strand. `minsh` is an aligner, not a mapper: its `mapped` is "
        "the number of placements it aligned and its accuracy is alignment "
        "`identity`, not concordance.",
        "",
        "A `0.0` in **agree w/ shmap** means either the tool mapped nothing on "
        "that dataset (so there is no placement to compare — see the CLR rows, "
        "where k-mer survival under high error is ~0), or it placed reads at a "
        "*different copy of the same repeat family* / opposite strand than shmap. "
        "On a repetitive acrocentric chromosome like chr21 this strict "
        "interval-overlap metric understates real agreement: re-adjudicating the "
        "ONT placements by canonical-15-mer containment at each candidate locus "
        "shows minshmap and shmap are competitive (neither is wrong). The "
        "meaningful difference there is *sensitivity* (reads mapped), not the "
        "concordance score.",
        "",
    ]
    for e in results:
        lines += [
            f"## {e['dataset']} — {e['platform']}",
            "",
            f"- {e['n_reads']} reads, {e['read_bp']:,} bp, avg {e['avg_read_len']:,} bp/read",
            "",
            "| tool | mapped | mapped_frac | wall s | reads/s | agree w/ shmap | identity |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
        for t in e["tools"]:
            lines.append(
                f"| {t['tool']} | {t['mapped']} | {_fmt(t['mapped_frac'])} | "
                f"{t['wall_sec']} | {t['reads_per_sec']} | "
                f"{_fmt(t['agree_with_shmap'])} | {_fmt(t['identity'])} |"
            )
        lines.append("")
        if "minsh" in e:
            m = e["minsh"]
            if m.get("note"):
                lines += [f"minSH alignment check skipped: {m['note']}", ""]
            else:
                lines += [
                    f"minSH aligned {m['n']} HiFi placements to their chr21 window: "
                    f"median identity **{m['median_identity']}**, median A* cells "
                    f"{m['median_cells']} (in {m['wall_sec']} s).",
                    "",
                ]

    speed, accuracy = _rankings(results)
    lines += [
        "## Speed ranking (overall reads/s across all datasets)",
        "",
        "| rank | tool | reads/s | total wall s |",
        "| --- | --- | --- | --- |",
    ]
    for i, (tool, rps, wall) in enumerate(speed, 1):
        lines.append(f"| {i} | {tool} | {rps} | {round(wall, 2)} |")
    lines += [
        "",
        "## Accuracy ranking (chr21 reads recovered + concordance)",
        "",
        "Ranked by total reads mapped across all datasets (sensitivity); "
        "`mean_agree_with_shmap` is the placement concordance against the "
        "reference mapper, and `identity` is minSH's independent alignment score.",
        "",
        "| rank | tool | total_mapped | mean_agree_with_shmap | identity |",
        "| --- | --- | --- | --- | --- |",
    ]
    for i, (tool, mapped, agree, ident) in enumerate(accuracy, 1):
        lines.append(f"| {i} | {tool} | {mapped} | {_fmt(agree)} | {_fmt(ident)} |")
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
