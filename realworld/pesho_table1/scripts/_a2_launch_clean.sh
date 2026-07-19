#!/usr/bin/env bash
# Archive the old (invalid, buggy-param) results, then launch a fresh CLEAN benchmark
# run on a2 with the fixed benchmark.py. Detached so it survives SSH disconnect.
set -eu
BASE=~/Pesho/minshmap/realworld/pesho_table1
SCR="$BASE/scripts"
RES="$BASE/results"
WORK=~/_paper_work

# 1) archive prior CSV/log (keep as record, out of the way of the clean run)
mkdir -p "$RES/old-incomplete"
shopt -s nullglob
for f in "$RES"/table1_*.csv "$RES"/run_a2_*.log; do
  mv -f "$f" "$RES/old-incomplete/" && echo "archived $(basename "$f")"
done

# 2) clear stale per-mapper scratch from the killed run
rm -f "$WORK/_mapshmap.paf" "$WORK/_time.txt" 2>/dev/null || true

# 3) launch fresh, fully detached
STAMP=$(date +%Y%m%d-%H%M%S)
LOG="$RES/run_a2_${STAMP}.log"
cd "$SCR"
setsid nohup python3 -u benchmark.py --threads 1 --datasets all > "$LOG" 2>&1 < /dev/null &
PID=$!
sleep 1
echo "LAUNCHED pid=$PID"
echo "log=$LOG"
echo "== process =="
ps -o pid,etime,args -p "$PID" || echo "(pid not visible yet)"
echo "== log head =="
head -20 "$LOG" 2>/dev/null || true
