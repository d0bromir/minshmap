#!/usr/bin/env bash
# Fetch a SMALL real HiFi FASTQ (with quality strings) to serve as the PBSIM3
# --method sample profile. Streams only the first N reads from the GIAB HG002
# PacBio CCS (HiFi) run -- the same source data_rw/hifi.fa was built from
# (see repo memory round 32) -- so PBSIM3 learns a genuine HiFi length+error profile.
#
# Usage:  N=2000 bash realworld/20a_fetch_hifi_sample.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$HERE/data_rw"
N="${N:-2000}"                       # reads to sample (4 FASTQ lines each)
OUT="$DATA/hifi_sample.fastq"
URL="${URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/m54238_180901_011437.Q20.fastq}"

echo "[20a] streaming first $N reads from GIAB HiFi CCS -> $OUT"
mkdir -p "$DATA"
# Stream and stop after 4*N lines; the pipe SIGPIPEs curl once head is done.
curl -fsSL "$URL" | head -n $((4 * N)) > "$OUT" || true
r=$(($(wc -l < "$OUT") / 4))
echo "[20a] wrote $r reads to $OUT ($(du -h "$OUT" | cut -f1))"
head -1 "$OUT"
