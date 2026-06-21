#!/usr/bin/env bash
# Probe the WSL environment for the real-world benchmark prerequisites.
set -u
ok()   { command -v "$1" >/dev/null 2>&1 && echo "$1: yes" || echo "$1: NO"; }

echo "=== toolchain ==="
g++ --version | head -1
make --version | head -1
git --version
echo "cores: $(nproc)"

echo "=== libs ==="
if dpkg -s zlib1g-dev >/dev/null 2>&1; then echo "zlib1g-dev: yes"; else echo "zlib1g-dev: MISSING"; fi

echo "=== tools ==="
for t in curl wget awk zcat python3 samtools minimap2 seqkit; do ok "$t"; done

echo "=== network ==="
curl -sI https://github.com            2>&1 | head -1
curl -sI https://hgdownload.soe.ucsc.edu 2>&1 | head -1
curl -sI https://ftp.sra.ebi.ac.uk     2>&1 | head -1
