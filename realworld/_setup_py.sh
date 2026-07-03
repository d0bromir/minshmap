#!/usr/bin/env bash
# Set up the WSL Python side for minshmap.py:
#   - biopython (Bio.SeqIO)
#   - minimizer_ext extension module (.so) built from the same Rust crate
# The .so is placed in ~/pylib so the benchmark can add it via PYTHONPATH without
# touching the Windows wheel install or the Windows cargo target.
set -uo pipefail
M=/mnt/c/Users/dobro/OneDrive/Desktop/PhD/Pesho/minshmap
PT=$HOME/minext_py
PYLIB=$HOME/pylib

echo "=== ensure pip ==="
if ! python3 -m pip --version >/dev/null 2>&1; then
    curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py
    python3 /tmp/get-pip.py --user --break-system-packages
fi
python3 -m pip --version

echo "=== install biopython ==="
python3 -m pip install --user --break-system-packages --quiet biopython
python3 -c "import Bio; print('biopython', Bio.__version__)"

echo "=== build minimizer_ext cdylib (.so, default/python features) ==="
cd "$M/minimizer_ext"
CARGO_TARGET_DIR="$PT" cargo build --release 2>&1 | tail -5
echo "cargo exit=${PIPESTATUS[0]}"
mkdir -p "$PYLIB"
cp -f "$PT/release/libminimizer_ext.so" "$PYLIB/minimizer_ext.so"
ls -la "$PYLIB/minimizer_ext.so"

echo "=== import check ==="
PYTHONPATH="$PYLIB" python3 -c "from minimizer_ext import canonical_minimizers; print('minimizer_ext OK', canonical_minimizers('ACGTACGTACGTACGTACGT',5,3)[:2])"

echo "=== minshmap.py smoke test ==="
cd "$M"
PYTHONPATH="$PYLIB" python3 minshmap.py data/ref.fa data/reads.fa -k 15 -w 5 -t 0.1 2>/tmp/py.err | head -2
echo "run exit=${PIPESTATUS[0]}"
head -5 /tmp/py.err
echo "=== DONE ==="
