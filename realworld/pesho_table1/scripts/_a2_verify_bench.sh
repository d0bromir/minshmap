#!/usr/bin/env bash
# Verify the synced benchmark.py on a2 has the fixed shmap params.
set -u
F=~/Pesho/minshmap/realworld/pesho_table1/scripts/benchmark.py
echo "== md5 on a2 =="
md5sum "$F"
echo "== fixed-param markers =="
grep -nE 'shmap-k|shmap-r|shmap-t|shmap-d|shmap-o|def run_map_shmap|min_diff|max_overlap' "$F" | head -30
