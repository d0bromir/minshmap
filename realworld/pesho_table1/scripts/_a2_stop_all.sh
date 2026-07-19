#!/usr/bin/env bash
# Decisively stop the a2 benchmark run: the python driver AND any active mapper child.
set -u
PATTERN='benchmark\.py --threads|/release/shmap|/bin/minimap2|winnowmap|/meryl|/blend|mapquik|minshmap_linux|_paper_work/_time.txt'
echo "== before =="
pgrep -af "$PATTERN" | grep -v grep || echo "(none)"

for p in 'benchmark\.py --threads' '/release/shmap' '/bin/minimap2' 'winnowmap' '/meryl' '/blend' 'mapquik' 'minshmap_linux' '_paper_work/_time.txt'; do
  pkill -KILL -f "$p" 2>/dev/null
done

echo "== after =="
if pgrep -af "$PATTERN" | grep -v grep; then
  echo "!! still running"
else
  echo "STOPPED-CLEAN"
fi
