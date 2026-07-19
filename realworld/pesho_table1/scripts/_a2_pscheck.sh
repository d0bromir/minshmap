#!/usr/bin/env bash
# Check for any running benchmark.py / mapper processes on a2.
set -u
echo "== benchmark.py / mapper processes =="
ps -eo pid,ppid,etime,rss,comm,args --sort=start_time \
  | grep -E 'benchmark\.py|shmap|minimap2|winnowmap|meryl|blend|mapquik' \
  | grep -v grep || echo "(none running)"
echo
echo "== newest results CSV + log =="
ls -lt ~/Pesho/minshmap/realworld/pesho_table1/results/*.csv 2>/dev/null | head -3
ls -lt ~/Pesho/minshmap/realworld/pesho_table1/results/*.log 2>/dev/null | head -3
