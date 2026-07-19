#!/usr/bin/env bash
# Decisively stop the invalid a2 benchmark run (buggy shmap params) and verify.
set -u
echo "== before =="
pgrep -af 'benchmark\.py|/Pesho/shmap/release/shmap' | grep -v grep || echo "(none)"

# SIGKILL the python driver, the shmap process, and its sh/time/bash wrappers
# (wrappers all carry the shmap path in their args -> matched by the 2nd pattern).
pkill -KILL -f 'benchmark\.py --threads 1' 2>/dev/null
pkill -KILL -f '/Pesho/shmap/release/shmap' 2>/dev/null
pkill -KILL -f '_paper_work/_time.txt' 2>/dev/null

echo "== after =="
if pgrep -af 'benchmark\.py|/Pesho/shmap/release/shmap' | grep -v grep; then
  echo "!! still running"
else
  echo "STOPPED-CLEAN"
fi
