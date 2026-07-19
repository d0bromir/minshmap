#!/usr/bin/env bash
# Full health + progress report for the current benchmark run.
set -u
BASE=~/Pesho/minshmap/realworld/pesho_table1
RES="$BASE/results"
LOG=$(ls -t "$RES"/run_a2_*.log 2>/dev/null | head -1)
CSV=$(ls -t "$RES"/table1_*.csv 2>/dev/null | head -1)

echo "########## DRIVER ##########"
pgrep -af 'benchmark\.py --threads' | grep -v grep || echo "DRIVER NOT RUNNING"

echo
echo "########## ACTIVE MAPPER (cmd + threads + elapsed + RSS) ##########"
pgrep -af '/bin/minimap2|winnowmap|/meryl|/release/shmap|/blend|mapquik|minshmap_linux' | grep -v grep | grep -vE '/bin/sh -c|/usr/bin/time|bash -c' || echo "(no mapper process this instant)"
echo "-- resource (top single proc) --"
ps -eo pid,etimes,rss,pcpu,comm --sort=-rss | awk 'NR==1 || /minimap2|winnowmap|meryl|shmap|blend|mapquik/' | head -6

echo
echo "########## COMPLETED CELLS (CSV: ${CSV:-none}) ##########"
if [ -n "${CSV:-}" ]; then column -s, -t "$CSV"; else echo "(no CSV yet — first mapper still running)"; fi

echo
echo "########## LOG TAIL (${LOG:-none}) ##########"
[ -n "${LOG:-}" ] && tail -15 "$LOG"

echo
echo "########## LIVE map-shmap PAF (if shmap running) ##########"
P=~/_paper_work/_mapshmap.paf
[ -f "$P" ] && printf "  _mapshmap.paf lines=%s size=%s\n" "$(wc -l < "$P")" "$(du -h "$P" | cut -f1)" || echo "  (map-shmap not at work)"
