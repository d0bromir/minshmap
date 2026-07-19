#!/usr/bin/env bash
# One-shot health check of the fresh benchmark run.
set -u
LOG=$(ls -t ~/Pesho/minshmap/realworld/pesho_table1/results/run_a2_*.log 2>/dev/null | head -1)
echo "== driver process =="
pgrep -af 'benchmark\.py --threads 1' | grep -v grep || echo "DRIVER DEAD"
echo "== current mapper =="
pgrep -af 'minimap2|/release/shmap|winnowmap|meryl|blend|mapquik|minshmap_linux' | grep -v grep || echo "(no mapper active this instant)"
echo "== log tail ($LOG) =="
tail -10 "$LOG" 2>/dev/null || echo "(no log)"
