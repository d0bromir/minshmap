#!/usr/bin/env bash
# Build the Linux minshmap binary (weighted C++ source) against a freshly built
# Linux static minimizer_ext lib. Uses a separate cargo target dir so the Windows
# artifacts under minimizer_ext/target/ are left untouched.
set -uo pipefail
M=/mnt/c/Users/dobro/OneDrive/Desktop/PhD/Pesho/minshmap
LT=$HOME/minext_linux_target

echo "=== cargo build (linux staticlib, --no-default-features) ==="
cd "$M/minimizer_ext"
CARGO_TARGET_DIR="$LT" cargo build --release --no-default-features 2>&1 | tail -6
echo "cargo exit=${PIPESTATUS[0]}"
ls -la "$LT/release/libminimizer_ext.a"

echo "=== link minshmap_linux ==="
cd "$M"
g++ -O3 -std=c++17 -march=native -pthread \
    -I ../shmap/ext/unordered_dense/include \
    -o minshmap_linux minshmap.cpp \
    -L "$LT/release" -l:libminimizer_ext.a \
    -lpthread -ldl -lm 2>&1 | head -20
echo "link exit=${PIPESTATUS[0]}"
ls -la minshmap_linux
echo "=== smoke test ==="
./minshmap_linux data/ref.fa data/reads.fa -k 15 -w 5 -t 0.1 -j 1 2>/tmp/mm.err | head -2
echo "run exit=${PIPESTATUS[0]}"
head -5 /tmp/mm.err
echo "=== DONE ==="
