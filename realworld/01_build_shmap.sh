#!/usr/bin/env bash
# Build Pesho's shmap (release) inside WSL/Linux.
#
# shmap needs C++2a, zlib, OpenMP and five git submodules. Three of them are
# registered with git@github.com SSH URLs in .gitmodules; we rewrite those to
# HTTPS so they clone without SSH keys. Only the submodules actually #included
# by the release build are fetched (pdqsort / tensor-sketching are skipped).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHMAP_DIR="$(cd "$HERE/../../shmap" && pwd)"
echo "shmap dir: $SHMAP_DIR"
cd "$SHMAP_DIR"

# Clone submodules over HTTPS instead of SSH (no keys needed), shallow.
git config url."https://github.com/".insteadOf "git@github.com:"
NEEDED="ext/cmd_line_parser ext/unordered_dense ext/gtl ext/klib ext/tracy"
echo "init submodules: $NEEDED"
git submodule update --init --depth 1 $NEEDED

# Disable the Tracy profiler in the release build. Tracy is compiled in via
# `-DTRACY_ENABLE`; with no profiler server draining its event queue it buffers
# every zone in RAM and grows without bound, OOM-killing shmap on long-read
# datasets. Commenting the flag turns the Tracy macros into no-ops. Idempotent.
if grep -qE '^\s*CFLAGS \+= -DTRACY_ENABLE' Makefile; then
    sed -i -E 's/^(\s*CFLAGS \+= -DTRACY_ENABLE)/#\1  # disabled for benchmarking (unbounded RAM)/' Makefile
    echo "disabled -DTRACY_ENABLE in Makefile"
fi

echo "building release (make -j$(nproc)) ..."
make release -j"$(nproc)"

echo "=== shmap built ==="
./release/shmap -h 2>&1 | head -30 || true
