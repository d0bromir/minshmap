#!/usr/bin/env bash
# Full real-world benchmark run. Unified sample size for every tool so counts
# are directly comparable; large enough to include on-target chr21 reads, small
# enough that the pure-Python mappers finish in reasonable time.
cd "$(dirname "$0")"
timeout 9000 python3 05_run_benchmark.py --max-reads 2000 > run_full.log 2>&1
echo "EXIT=$?"
echo "==== TAIL run_full.log ===="
tail -70 run_full.log
