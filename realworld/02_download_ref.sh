#!/usr/bin/env bash
# Download the T2T-CHM13v2.0 chromosome 21 reference (~45 Mb) only.
#
# Downloading the full analysis-set FASTA (~1 GB) just to extract one chromosome
# is wasteful, so we look up chr21's GenBank accession from the assembly report
# and stream just that single sequence from NCBI E-utilities.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data_rw"
mkdir -p "$DATA"
OUT="$DATA/chr21.fa"

if [[ -s "$OUT" ]]; then
    echo "chr21 already present: $OUT ($(wc -c <"$OUT") bytes)"; exit 0
fi

REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt"
echo "fetching assembly report ..."
# Match the assembled chromosome 21 row; strip CRs (report uses CRLF).
ACC="$(curl -fsSL "$REPORT_URL" | tr -d '\r' \
      | awk -F'\t' '$2=="assembled-molecule" && $3=="21" {print $5; exit}')"
if [[ -z "${ACC:-}" ]]; then echo "could not find chr21 accession" >&2; exit 1; fi
echo "chr21 GenBank accession: $ACC"

echo "downloading chr21 FASTA ..."
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text"
# Normalise the header to a clean ">chr21".
curl -fsSL "$EFETCH" | awk 'NR==1{print ">chr21"; next} {print}' > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"

BP="$(grep -v '^>' "$OUT" | tr -d '\n' | wc -c)"
echo "wrote $OUT  (${BP} bp)"
