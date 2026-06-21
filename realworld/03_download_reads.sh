#!/usr/bin/env bash
# Download small subsets of three real, recently-sequenced long-read datasets
# for HG002/NA24385, one per common platform (the error-rate spectrum):
#
#   hifi : PacBio HiFi  (Sequel II,  ~99.9% acc)   SRR17284858
#   ont  : Oxford Nanopore (GridION, ~95% acc)     SRR11032657
#   clr  : PacBio CLR  (RS II subreads, ~87% acc)  SRR2036394
#
# We do NOT download the whole multi-GB files: we stream the gzip and stop after
# the first N reads (head closes the pipe -> curl/zcat get SIGPIPE and stop), so
# only a few hundred MB at most is transferred. FASTQ is converted to FASTA.
#
# Usage: bash 03_download_reads.sh [N_READS]   (default 6000)
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data_rw"
mkdir -p "$DATA"
N="${1:-6000}"

# label  accession
DATASETS=(
  "hifi SRR17284858"
  "ont  SRR11032657"
  "clr  SRR2036394"
)

resolve_url () {  # accession -> first https fastq url
  local acc="$1"
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv" \
    | awk -F'\t' 'NR>1{for(i=1;i<=NF;i++) if($i ~ /ftp\.sra\.ebi/){split($i,a,";"); print a[1]; exit}}'
}

for entry in "${DATASETS[@]}"; do
  set -- $entry; label="$1"; acc="$2"
  out_fa="$DATA/${label}.fa"
  if [[ -s "$out_fa" ]]; then
    echo "[$label] already present ($(grep -c '^>' "$out_fa") reads)"; continue
  fi
  echo "[$label] resolving $acc ..."
  url="$(resolve_url "$acc")"
  if [[ -z "$url" ]]; then echo "[$label] could not resolve fastq url" >&2; continue; fi
  echo "[$label] streaming first $N reads from https://$url"
  # FASTQ -> FASTA, keeping the run accession in the read name. head stops early.
  curl -fsSL "https://$url" 2>/dev/null \
    | zcat 2>/dev/null \
    | head -n $((N*4)) \
    | awk -v acc="$acc" 'NR%4==1{printf(">%s_%d\n", acc, ++i)} NR%4==2{print}' \
    > "$out_fa.tmp" || true
  mv "$out_fa.tmp" "$out_fa"
  reads=$(grep -c '^>' "$out_fa"); reads=${reads:-0}
  bp=$(grep -v '^>' "$out_fa" | tr -d '\n' | wc -c)
  avg=$(awk -v b="$bp" -v r="$reads" 'BEGIN{print (r>0)?int(b/r):0}')
  echo "[$label] wrote $out_fa : $reads reads, $bp bp (avg $avg bp/read)"
done
