"""Optional benchmark instrumentation for minSHmap - kept OUT of the pedagogical mapper.

The mapper (minshmap.py) stays minimal: it just does `b = Bench()`, calls `b.mark("index")`
and `b.mark("map")` at the phase boundaries, then `b.report(reads, mapped)`. All the timing,
/proc RSS parsing and formatting live here so they don't clutter the ~150-line teaching file.

Everything is a cheap no-op unless the environment variable MINSHMAP_BENCH is set; when it is,
exactly ONE line goes to STDERR (stdout PAF is never touched, so py/cpp stay byte-identical):

    [bench] index_s=.. map_s=.. reads=.. mapped=.. index_rss_mb=.. peak_rss_mb=..

The real-world benchmark (realworld/11_bench_3way.py) sets MINSHMAP_BENCH and regex-parses
that line for the per-phase timing/memory split.
"""
import os
import sys
import time


def _rss_mb():
    """(current VmRSS, peak VmHWM) in MB from /proc/self/status (Linux); (0, 0) elsewhere."""
    cur = peak = 0.0
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    cur = int(line.split()[1]) / 1024.0
                elif line.startswith("VmHWM:"):
                    peak = int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return cur, peak


class Bench:
    """Phase timer + resident-memory probe. `mark(label)` closes a phase (recording its
    duration and the current RSS); `report(reads, mapped)` prints the summary line. All
    methods are no-ops unless MINSHMAP_BENCH is set, so the mapper pays nothing by default."""

    def __init__(self):
        self.on = bool(os.environ.get("MINSHMAP_BENCH"))
        self._lap = time.perf_counter()
        self.secs = {}                                   # label -> phase seconds
        self.rss = {}                                    # label -> current VmRSS (MB) at that mark

    def mark(self, label):
        if not self.on:
            return
        now = time.perf_counter()
        self.secs[label] = now - self._lap
        self.rss[label] = _rss_mb()[0]
        self._lap = now                                  # next phase starts here

    def report(self, reads, mapped):
        if not self.on:
            return
        peak = _rss_mb()[1]
        sys.stderr.write(
            f"[bench] index_s={self.secs.get('index', 0):.6f} "
            f"map_s={self.secs.get('map', 0):.6f} reads={reads} mapped={mapped} "
            f"index_rss_mb={self.rss.get('index', 0):.6f} peak_rss_mb={peak:.6f}\n")
