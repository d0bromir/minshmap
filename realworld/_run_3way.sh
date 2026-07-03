#!/usr/bin/env bash
# Launch the full 3-way overnight benchmark with a clean environment. All paths
# and env vars are resolved inside bash to avoid PowerShell/WSL quoting issues.
set -uo pipefail
cd /mnt/c/Users/dobro/OneDrive/Desktop/PhD/Pesho/minshmap
export PYTHONPATH="$HOME/pylib"
LOG=realworld/results_rw/bench_full.log
echo "start $(date)" > "$LOG"
python3 -u realworld/11_bench_3way.py --scope both --datasets hifi ont clr \
        --max-reads 6000 --threads 4 >> "$LOG" 2>&1
echo "done $(date) rc=$?" >> "$LOG"
