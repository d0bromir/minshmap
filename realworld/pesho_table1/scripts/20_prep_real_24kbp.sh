#!/usr/bin/env bash
# Type 4: All chromosomes, REAL HiFi reads from HG002 (no ground truth).
# Reference = whole CHM13v2.0. Reads = real HG002 HiFi subset (hifi_sample.fastq)
# converted to FASTA. Because there is no ground truth, the benchmark reports
# Wrong Q60 = n/a and Mapped Q60 = number of reads with mapq = 60.
#
#   N=<count>  -> use only the first N reads (default: all reads in the sample)
# To fetch a larger real set first, run realworld/pesho_table1/scripts/20a_fetch...
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SAMPLE="${HIFI_SAMPLE:-$SD/../../data_rw/hifi_sample.fastq}"
OUT="$SD/../data/allchr_real_24kbp"
mkdir -p "$OUT"
N="${N:-0}"
[ -f "$SAMPLE" ] || { echo "missing real HiFi sample: $SAMPLE" >&2; exit 1; }

echo "[type 4] all chromosomes | REAL HG002 HiFi reads (no truth) -> reads.fa"
if [ "$N" -gt 0 ]; then
  head -n $((N * 4)) "$SAMPLE"
else
  cat "$SAMPLE"
fi | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' > "$OUT/reads.fa"
echo "reads: $(grep -c '^>' "$OUT/reads.fa")  ->  $OUT/reads.fa"
