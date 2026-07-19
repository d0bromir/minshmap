#!/usr/bin/env bash
set -uo pipefail
. "$HOME/.cargo/env" 2>/dev/null || true
REPO=/home/mpiuser/Pesho/minshmap/realworld/pesho_table1
cd "$REPO"
echo "===== mapquik ONLY, all 4 datasets (working binary + oneline ref) ====="
python3 scripts/benchmark.py --only mapquik --datasets all 2>&1 | tail -40
echo; echo "===== resulting mapquik rows ====="
NEW=$(ls -t results/table1_*.csv | head -1)
echo "CSV: $NEW"
grep -H "dataset,\|mapquik" "$NEW"
