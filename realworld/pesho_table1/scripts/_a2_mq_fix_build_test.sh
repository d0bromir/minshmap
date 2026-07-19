#!/usr/bin/env bash
set -uo pipefail
. "$HOME/.cargo/env" 2>/dev/null || true
export CARGO_TARGET_DIR="$HOME/tools/mq_build"
NIGHTLY=nightly-2026-02-08
F="$HOME/tools/rust-seq2kminmers"
M="$HOME/tools/mapquik"

echo "===== 1. patch compressstoreu base_addr *mut u8 -> *mut i32 ====="
for src in "$F/src/hpc.rs" "$F/src/nthash_avx512_32.rs" "$F/src/nthash2_avx512_32.rs"; do
  sed -i '/_mm512_mask_compressstoreu_epi32/ s/as \*mut u8/as *mut i32/g' "$src"
done
grep -rn "compressstoreu_epi32" "$F/src" | sed 's/^/  /'

echo; echo "===== 2. override mapquik git dep -> local patched fork ====="
if ! grep -q 'patch."https://github.com/rchikhi/rust-seq2kminmers"' "$M/Cargo.toml"; then
cat >> "$M/Cargo.toml" <<EOF

[patch."https://github.com/rchikhi/rust-seq2kminmers"]
rust-seq2kminmers = { path = "$F" }
EOF
fi
tail -4 "$M/Cargo.toml"

echo; echo "===== 3. build (fat LTO) ====="
cd "$M"
if cargo "+$NIGHTLY" build --release 2>build.log; then
  echo "BUILD_OK"
else
  echo "BUILD_FAIL -- last 30 lines:"; tail -30 build.log; exit 1
fi
BIN="$CARGO_TARGET_DIR/release/mapquik"
ls -la "$BIN"

echo; echo "===== 4. functional test: 100 near-perfect ecoli reads -> ecoli genome ====="
cd "$M/example"
READS=nearperfect-ecoli.100.fa
REF=ecoli.genome.fa
run() { # $1=label $2..=extra flags
  local label="$1"; shift
  "$BIN" "$READS" --reference "$REF" -k 8 -d 0.01 -l 16 -p "out_$label" -g 100 --threads 11 "$@" >/dev/null 2>"log_$label.txt" || echo "  ($label exited non-zero)"
  local n=0; [ -f "out_$label.paf" ] && n=$(wc -l < "out_$label.paf")
  echo "  MODE $label: mapped_lines=$n / 100 reads"
}
run "hpcsimd"
run "hpc"    --nosimd
run "regular" --nosimd --nohpc
echo; echo "=== sample of best PAF (hpcsimd) ==="
head -2 "out_hpcsimd.paf" 2>/dev/null || echo "(empty)"
