"""10_bench_wgs.py - align a real WGS read set to the WHOLE human genome with minSHmap.

Downloads the T2T-CHM13v2.0 assembly (UCSC "hs1", one gzipped FASTA of all 24
chromosomes + chrM, ~1 GB) and runs the C++ binary (minshmap.exe / minshmap_linux)
on it for each WGS read set in data_rw/ (hifi, ont, clr - the same reads we used
against chr21, except now mapped against the FULL genome, so reads from any
chromosome can place). Times index build + mapping, counts mapped reads and mean
mapq, writes a CSV, prints a table.

WHY C++ ONLY: the pure-Python tool cannot index the whole genome. minSHmap's index
is a memory-heavy CSR (during build it holds raw (hash,hit) pairs AND the final hit
vector, ~40 bytes per minimizer at peak). The whole genome has ~3.1e9 * 2/(w+1)
minimizers; even in C++ that is several GB, and in CPython it is ~10x worse. So this
benchmark drives the C++ binary, and picks `w` large enough that the index fits in RAM.

MEMORY KNOB: minimizer density is ~2/(w+1), so a larger (odd) `w` => fewer index
entries => less RAM, at some cost in sensitivity (minimap2 does the same for big
genomes). The script estimates peak index RAM from the genome size and your chosen
`w` and refuses to start if it would not fit, telling you the smallest safe `w`.

Usage (Windows):
    C:\\Python314\\python.exe 10_bench_wgs.py --datasets hifi ont clr -w 31
    C:\\Python314\\python.exe 10_bench_wgs.py --max-reads 2000 -w 51   # quick, low-RAM
"""
import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data_rw")
RESULTS = os.path.join(HERE, "results_rw")

# UCSC "hs1" = T2T-CHM13v2.0, the whole human genome as a single gzipped FASTA.
GENOME_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
GENOME_GZ = os.path.join(DATA, "hs1.fa.gz")
GENOME_FA = os.path.join(DATA, "hs1.fa")

# Per-platform homology threshold (read error rate differs): same presets as the
# chr21 benchmark. Tune if a larger w shifts the score distribution.
THETA = {"hifi": 0.20, "ont": 0.15, "clr": 0.18}


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024


def binary():
    """The C++ minSHmap binary for this OS (built next to minshmap.cpp)."""
    root = os.path.dirname(HERE)
    exe = os.path.join(root, "minshmap.exe" if os.name == "nt" else "minshmap_linux")
    if not os.path.exists(exe):
        sys.exit(f"C++ binary not found: {exe}\n  build it first (see minshmap.cpp header).")
    return exe


def download_genome():
    """Fetch hs1.fa.gz (resumable-free, simple stream) and gunzip it to hs1.fa, once."""
    if os.path.exists(GENOME_FA):
        print(f"genome present: {GENOME_FA} ({human(os.path.getsize(GENOME_FA))})")
        return
    os.makedirs(DATA, exist_ok=True)
    if not os.path.exists(GENOME_GZ):
        print(f"downloading {GENOME_URL}\n  -> {GENOME_GZ}  (~1 GB, this takes a while)...")
        t0 = time.time()
        with urllib.request.urlopen(GENOME_URL) as r, open(GENOME_GZ + ".part", "wb") as out:
            total = int(r.headers.get("Content-Length", 0))
            done = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
                done += len(chunk)
                if total:
                    print(f"\r  {human(done)} / {human(total)} ({100*done/total:.0f}%)", end="")
        os.replace(GENOME_GZ + ".part", GENOME_GZ)
        print(f"\n  downloaded in {time.time()-t0:.0f}s")
    print(f"decompressing -> {GENOME_FA} ...")
    with gzip.open(GENOME_GZ, "rb") as gz, open(GENOME_FA + ".part", "wb") as out:
        shutil.copyfileobj(gz, out, 1 << 20)
    os.replace(GENOME_FA + ".part", GENOME_FA)
    print(f"  genome ready: {human(os.path.getsize(GENOME_FA))}")


def check_memory(genome_bp, w):
    """Estimate peak index-build RAM (~40 bytes per minimizer, density ~2/(w+1)) and
    bail with the smallest safe odd w if it would not fit in physical memory."""
    minimizers = genome_bp * 2 / (w + 1)
    peak_gb = minimizers * 40 / 1e9
    avail_gb = None
    try:
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        if os.name == "nt":
            ms = MS(); ms.dwLength = ctypes.sizeof(MS)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
            total_gb = ms.ullTotalPhys / 1e9
        else:
            total_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
        avail_gb = total_gb
    except Exception:
        total_gb = None
    print(f"index estimate (w={w}): ~{minimizers/1e6:.0f}M minimizers, "
          f"peak build RAM ~{peak_gb:.1f} GB"
          + (f"; machine has {total_gb:.1f} GB physical" if total_gb else ""))
    if total_gb and peak_gb > 0.8 * total_gb:
        # smallest odd w whose peak fits under 80% of RAM
        need = genome_bp * 2 * 40 / 1e9 / (0.8 * total_gb)  # = w+1
        safe = int(need) + 1
        safe += (safe % 2 == 0)  # make odd
        sys.exit(f"\nABORT: estimated index RAM (~{peak_gb:.1f} GB) exceeds 80% of "
                 f"physical RAM ({total_gb:.1f} GB).\n  Re-run with a larger odd window, "
                 f"e.g. -w {safe} (sparser minimizers), or use --max-reads on a smaller "
                 f"machine.\n  (mapping cost is unaffected by genome size; this is purely "
                 f"the index.)")


def count_genome_bp(path):
    """Total non-header bytes (approx bp) of the reference FASTA."""
    bp = 0
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b">"):
                bp += len(line.strip())
    return bp


def subset_reads(src, max_reads):
    """First max_reads records of a FASTA -> a temp file; returns its path (or src)."""
    if not max_reads:
        return src, None
    dst = src + f".first{max_reads}.fa"
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                n += 1
                if n > max_reads:
                    break
            fo.write(line)
    return dst, dst


def run_one(exe, ref, reads, k, w, theta, threads):
    """Run the C++ mapper, return (wall_s, n_reads, mapped, mean_mapq)."""
    n_reads = sum(1 for line in open(reads) if line.startswith(">"))
    t0 = time.time()
    proc = subprocess.run([exe, ref, reads, "-k", str(k), "-w", str(w),
                           "-t", str(theta), "-j", str(threads)],
                          capture_output=True, text=True)
    wall = time.time() - t0
    if proc.returncode != 0:
        print(f"  ! mapper failed (exit {proc.returncode}): {proc.stderr.strip()[:200]}")
        return wall, n_reads, 0, 0.0
    mapped = mapq_sum = 0
    for line in proc.stdout.splitlines():
        if not line:
            continue
        mapped += 1
        mapq_sum += int(line.rsplit("\t", 1)[-1])
    return wall, n_reads, mapped, (mapq_sum / mapped if mapped else 0.0)


def main():
    p = argparse.ArgumentParser(description="minSHmap whole-human-genome WGS benchmark (C++)")
    p.add_argument("--datasets", nargs="+", default=["hifi", "ont", "clr"])
    p.add_argument("-k", type=int, default=15)
    p.add_argument("-w", type=int, default=31, help="odd minimizer window; larger = less index RAM")
    p.add_argument("--max-reads", type=int, default=0, help="map only the first N reads per set")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    p.add_argument("--ref", default=GENOME_FA, help="reference FASTA (default: downloaded hs1.fa)")
    p.add_argument("--skip-download", action="store_true")
    a = p.parse_args()
    if a.w % 2 == 0:
        p.error("-w must be odd (canonical minimizers)")

    if a.ref == GENOME_FA and not a.skip_download:
        download_genome()
    if not os.path.exists(a.ref):
        sys.exit(f"reference not found: {a.ref}")

    exe = binary()
    print("counting genome size...", end=" ", flush=True)
    genome_bp = count_genome_bp(a.ref)
    print(f"{genome_bp/1e9:.2f} Gbp")
    check_memory(genome_bp, a.w)

    os.makedirs(RESULTS, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    out_csv = os.path.join(RESULTS, f"wgs_cpp_{stamp}.csv")
    rows = []
    for ds in a.datasets:
        reads = os.path.join(DATA, f"{ds}.fa")
        if not os.path.exists(reads):
            print(f"skip {ds}: {reads} missing")
            continue
        reads, tmp = subset_reads(reads, a.max_reads)
        theta = THETA.get(ds, 0.2)
        print(f"\n=== {ds} (k={a.k} w={a.w} t={theta}) on whole genome ===")
        wall, n, mapped, mq = run_one(exe, a.ref, reads, a.k, a.w, theta, a.threads)
        rps = n / wall if wall else 0
        print(f"  reads {n}  mapped {mapped} ({100*mapped/n:.1f}%)  "
              f"{wall:.1f}s ({rps:.1f} r/s)  mean mapq {mq:.1f}")
        rows.append({"dataset": ds, "reads": n, "k": a.k, "w": a.w, "theta": theta,
                     "wall_s": round(wall, 3), "reads_per_s": round(rps, 2),
                     "mapped": mapped, "mapped_frac": round(mapped / n, 4) if n else 0,
                     "mean_mapq": round(mq, 2)})
        if tmp:
            os.remove(tmp)

    if rows:
        with open(out_csv, "w", newline="") as f:
            wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wtr.writeheader()
            wtr.writerows(rows)
        print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
