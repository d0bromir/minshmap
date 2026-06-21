#!/usr/bin/env bash
# Prepare the WSL side for the real-world harness:
#   - build a Linux minshmap binary from the same C++ source
#   - create a venv with numpy + fenwick (needed by minSH's A* alignment step)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"          # minshmap/

echo "=== building Linux minshmap binary ==="
g++ -O3 -std=c++17 -march=native -o "$ROOT/minshmap_linux" "$ROOT/minshmap.cpp"
"$ROOT/minshmap_linux" --demo --hash nthash | tail -1

echo "=== installing numpy + fenwick (minSH deps) without sudo ==="
# This Ubuntu python has no pip and venv needs apt; bootstrap pip into ~/.local.
if ! python3 -m pip --version >/dev/null 2>&1; then
    echo "bootstrapping pip via get-pip.py ..."
    curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py
    python3 /tmp/get-pip.py --user --break-system-packages
fi
python3 -m pip install --user --break-system-packages --quiet numpy fenwick
python3 -c "import numpy, fenwick; print('numpy', numpy.__version__, '+ fenwick OK')"
echo "=== setup complete ==="
