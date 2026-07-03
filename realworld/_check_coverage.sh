#!/usr/bin/env bash
# Diagnose the low map-rate: are reads simply whole-genome (so chr21 can only
# capture its ~1.5% share), or is the mapper missing reads it should catch?
# Map the SAME read subset against chr21 vs the whole genome (hs1.fa) and compare.
set -uo pipefail
cd /mnt/c/Users/dobro/OneDrive/Desktop/PhD/Pesho/minshmap
CPP=./minshmap_linux
DATA=realworld/data_rw
K=15; W=31; N=2000

declare -A THETA=( [hifi]=0.20 [ont]=0.15 [clr]=0.18 )

count() {  # count unique mapped query names on stdin
  cut -f1 | sort -u | wc -l
}

for ds in hifi ont clr; do
  src="$DATA/$ds.fa"
  sub="/tmp/sub_$ds.fa"
  head -n $((N*2)) "$src" > "$sub"
  reads=$(grep -c '^>' "$sub")
  t=${THETA[$ds]}

  m_chr21=$("$CPP" "$DATA/chr21.fa" "$sub" -k $K -w $W -t "$t" -j 4 2>/dev/null | count)
  t0=$(date +%s)
  m_wg=$("$CPP" "$DATA/hs1.fa" "$sub" -k $K -w $W -t "$t" -j 4 2>/dev/null | count)
  t1=$(date +%s)

  pc_chr21=$(awk -v m="$m_chr21" -v r="$reads" 'BEGIN{printf "%.1f", 100*m/r}')
  pc_wg=$(awk -v m="$m_wg" -v r="$reads" 'BEGIN{printf "%.1f", 100*m/r}')
  echo "$ds: reads=$reads  chr21 mapped=$m_chr21 (${pc_chr21}%)  wholegenome mapped=$m_wg (${pc_wg}%)  [wg took $((t1-t0))s]"
done
