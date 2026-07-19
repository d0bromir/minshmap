#!/usr/bin/env bash
# =============================================================================
# build_mappers.sh  --  Reproducible build recipe for EVERY mapper in Table 1
#                       of Ivanov & Medvedev, "map-shmap: Practical long-read
#                       mapping with seed heuristic on sketches".
#
# The paper compares five read mappers.  This script documents and performs the
# exact build of each one, from source, into ~/bin (no sudo / no root needed):
#
#   1. minimap2     github.com/lh3/minimap2        (C,   make)
#   2. Winnowmap2   github.com/marbl/Winnowmap     (C++, make)  + meryl
#   3. BLEND        github.com/CMU-SAFARI/BLEND    (C,   make)  minimap2 fork
#   4. mapquik      github.com/ekimb/mapquik       (Rust,cargo)
#   5. map-shmap    ../../../../shmap  (the ORIGINAL shmap in THIS repo; C++ make)
#
#   (minshmap -- our educational re-implementation -- is NOT a paper mapper and
#    is built separately; it already ships as minshmap/minshmap_linux.)
#
# Everything is cloned under ~/tools (ext4, fast).  Idempotent: a mapper whose
# binary already exists is skipped.  Run again any time to (re)build.
#
# Prerequisite handled automatically: Winnowmap's bundled meryl links against
# OpenSSL, whose -dev headers are missing on this box and cannot be apt-installed
# (no sudo).  Section 0 fetches + wires up OpenSSL entirely in user space.
# =============================================================================
set -euo pipefail

BIN="$HOME/bin"; SRC="$HOME/tools"
mkdir -p "$BIN" "$SRC"
JOBS="$(nproc)"
# Repo root of THIS checkout (…/Pesho), derived from this script's location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"   # scripts -> pesho_table1 -> realworld -> minshmap -> Pesho

have() { [ -x "$BIN/$1" ]; }

# -----------------------------------------------------------------------------
# 0. Build-time deps for Winnowmap's bundled meryl / htslib, all in user space.
#    meryl's main.mk resolves openssl + libcurl + liblzma via pkg-config and also
#    links -lz -lbz2 directly.  On this box:
#      * openssl -dev headers are missing  -> extracted from libssl-dev (no root)
#      * libcurl / liblzma -dev headers    -> taken from the existing miniconda3
#      * runtime .so's for curl/lzma/bz2/z -> already present in the system
#    We synthesise the three .pc files pointing at those, so meryl links cleanly.
# -----------------------------------------------------------------------------
setup_meryl_deps() {
  local DEV="$HOME/localdev/root"
  local INC="$DEV/usr/include"
  local LIBDIR="$DEV/usr/lib/x86_64-linux-gnu"
  local PCDIR="$LIBDIR/pkgconfig"
  local SYSLIB="/usr/lib/x86_64-linux-gnu"
  local CONDA_INC="$HOME/miniconda3/include"   # provides curl/curl.h, lzma.h

  # -- OpenSSL headers (extract libssl-dev .deb without root) ------------------
  if [ ! -f "$INC/openssl/hmac.h" ]; then
    echo "[deps] fetching libssl-dev (no root) ..."
    mkdir -p "$HOME/localdev/debs" "$DEV"
    ( cd "$HOME/localdev/debs" && apt-get download libssl-dev >/dev/null 2>&1 || true )
    for d in "$HOME"/localdev/debs/*.deb; do [ -e "$d" ] && dpkg-deb -x "$d" "$DEV"; done
  fi
  [ -f "$INC/openssl/hmac.h" ] || { echo "[deps] ERROR: openssl headers unavailable"; return 1; }

  # -- dev symlinks so -lssl/-lcrypto/-lcurl/-llzma resolve to runtime .so's ----
  mkdir -p "$LIBDIR"
  ln -sf "$SYSLIB/libssl.so.3"       "$LIBDIR/libssl.so"
  ln -sf "$SYSLIB/libcrypto.so.3"    "$LIBDIR/libcrypto.so"
  ln -sf "$SYSLIB/libcurl.so.4"      "$LIBDIR/libcurl.so"
  ln -sf "$SYSLIB/liblzma.so.5"      "$LIBDIR/liblzma.so"

  # -- synthesise pkg-config files ---------------------------------------------
  mkdir -p "$PCDIR"
  # openssl / libssl / libcrypto (headers from extracted deb)
  for pc in openssl libssl libcrypto; do
    cat > "$PCDIR/$pc.pc" <<EOF
prefix=$DEV/usr
libdir=$LIBDIR
includedir=$INC
Name: $pc
Description: OpenSSL (user-space)
Version: 3.0.13
Libs: -L\${libdir} -lssl -lcrypto
Cflags: -I\${includedir}
EOF
  done
  # libcurl (headers from miniconda, lib from system runtime)
  cat > "$PCDIR/libcurl.pc" <<EOF
libdir=$LIBDIR
includedir=$CONDA_INC
Name: libcurl
Description: libcurl (user-space)
Version: 8.0.0
Libs: -L\${libdir} -lcurl
Cflags: -I\${includedir}
EOF
  # liblzma (headers from miniconda, lib from system runtime)
  cat > "$PCDIR/liblzma.pc" <<EOF
libdir=$LIBDIR
includedir=$CONDA_INC
Name: liblzma
Description: liblzma (user-space)
Version: 5.4.5
Libs: -L\${libdir} -llzma
Cflags: -I\${includedir}
EOF

  export PKG_CONFIG_PATH="$PCDIR:${PKG_CONFIG_PATH:-}"
  export CPATH="$INC:$INC/x86_64-linux-gnu:$CONDA_INC:${CPATH:-}"
  export LIBRARY_PATH="$LIBDIR:$SYSLIB:${LIBRARY_PATH:-}"
  echo "[deps] openssl=$(pkg-config --modversion openssl 2>/dev/null) curl=$(pkg-config --exists libcurl && echo ok) lzma=$(pkg-config --exists liblzma && echo ok)"
}

# -----------------------------------------------------------------------------
# 1. minimap2   -- reference long-read mapper (baseline in the paper)
#      source : https://github.com/lh3/minimap2
#      build  : make          -> ./minimap2
# -----------------------------------------------------------------------------
build_minimap2() {
  if have minimap2; then echo "[skip] minimap2 already in $BIN"; return; fi
  echo "[build] minimap2 ..."
  cd "$SRC"
  [ -d minimap2 ] || git clone --depth 1 https://github.com/lh3/minimap2.git
  cd minimap2
  make -j"$JOBS"
  cp -f minimap2 "$BIN/minimap2"
  echo "[ok] minimap2 -> $BIN"
}

# -----------------------------------------------------------------------------
# 2. Winnowmap2 (+ meryl)  -- weighted-minimizer long-read mapper
#      source : https://github.com/marbl/Winnowmap
#      build  : make          -> bin/winnowmap , bin/meryl
#      needs  : OpenSSL (section 0) for the bundled meryl/htslib
# -----------------------------------------------------------------------------
build_winnowmap() {
  if have winnowmap && have meryl; then echo "[skip] winnowmap + meryl already in $BIN"; return; fi
  echo "[build] Winnowmap (winnowmap + meryl) ..."
  setup_meryl_deps
  cd "$SRC"
  [ -d Winnowmap ] || git clone --depth 1 https://github.com/marbl/Winnowmap.git
  cd Winnowmap
  make -j"$JOBS"
  cp -f bin/winnowmap "$BIN/winnowmap"
  cp -f bin/meryl     "$BIN/meryl"
  echo "[ok] winnowmap + meryl -> $BIN"
}

# -----------------------------------------------------------------------------
# 3. BLEND  -- strobemer-style mapper (minimap2 fork)
#      source : https://github.com/CMU-SAFARI/BLEND  (submodules required)
#      build  : make          -> bin/blend
# -----------------------------------------------------------------------------
build_blend() {
  if have blend; then echo "[skip] blend already in $BIN"; return; fi
  echo "[build] BLEND ..."
  cd "$SRC"
  [ -d BLEND ] || git clone --depth 1 --recursive https://github.com/CMU-SAFARI/BLEND.git
  cd BLEND
  make -j"$JOBS"
  cp -f bin/blend "$BIN/blend"
  echo "[ok] blend -> $BIN"
}

# -----------------------------------------------------------------------------
# 4. mapquik  -- k-min-mer mapper (Rust)
#      source : https://github.com/ekimb/mapquik
#               dep https://github.com/rchikhi/rust-seq2kminmers (rev a409c28)
#      build  : cargo +nightly build --release  -> target/release/mapquik
#      PROBLEM: with a *current* toolchain the dep rust-seq2kminmers fails to
#               compile because the AVX512 intrinsic _mm512_mask_compressstoreu_
#               epi32 changed its pointer arg from `*mut u8` to `*mut i32`
#               (stdarch, ~Rust 1.75).  An *old* nightly with the u8* signature
#               can't be used either: mapquik's modern transitive deps use
#               edition 2024, which old cargo cannot parse, and mapquik ships no
#               committed Cargo.lock to pin them.
#      FIX    : build a local, patched fork of rust-seq2kminmers (u8* -> i32* on
#               the uncommented compressstoreu_epi32 calls), point mapquik at it
#               with a [patch] override, and build with the CURRENT nightly
#               (needed for mapquik's `#![feature(iter_advance_by)]`).
#      toolchain: rustup (user-space, no root) installs the nightly on demand.
# -----------------------------------------------------------------------------
SEQ2KMINMERS_REV="a409c289fa55f014e73d00733e5e1a69317d94e3"
build_mapquik() {
  if have mapquik; then echo "[skip] mapquik already in $BIN"; return; fi
  echo "[build] mapquik (cargo +nightly, patched rust-seq2kminmers) ..."
  rustup toolchain list 2>/dev/null | grep -q '^nightly-' \
    || rustup toolchain install nightly --profile minimal

  # 4a. Local patched fork of the AVX512 seed dependency.
  cd "$SRC"
  if [ ! -d rust-seq2kminmers ]; then
    git clone https://github.com/rchikhi/rust-seq2kminmers.git
    ( cd rust-seq2kminmers && git checkout "$SEQ2KMINMERS_REV"
      # AVX512 mask-compress store now takes *mut i32 (was *mut u8).
      sed -i '/compressstoreu_epi32/ s/as \*mut u8/as *mut i32/g' \
        src/nthash_avx512_32.rs src/nthash2_avx512_32.rs src/hpc.rs )
  fi

  # 4b. mapquik itself, patched to use the local fork.
  [ -d mapquik ] || git clone --depth 1 https://github.com/ekimb/mapquik.git
  cd mapquik
  if ! grep -q 'patch."https://github.com/rchikhi/rust-seq2kminmers"' Cargo.toml; then
    cat >> Cargo.toml <<EOF

[patch."https://github.com/rchikhi/rust-seq2kminmers"]
rust-seq2kminmers = { path = "$SRC/rust-seq2kminmers" }
EOF
  fi
  rm -f Cargo.lock
  CARGO_TARGET_DIR="$SRC/mapquik_target" cargo +nightly build --release
  cp -f "$SRC/mapquik_target/release/mapquik" "$BIN/mapquik"
  echo "[ok] mapquik -> $BIN"
}

# -----------------------------------------------------------------------------
# 5. map-shmap  -- the ORIGINAL shmap (the tool under evaluation in the paper).
#      source : this repo, shmap/  (vendored deps under shmap/ext)
#      build  : make          -> release/shmap
#      The pre-built binary already ships at shmap/release/shmap; we skip if
#      present, otherwise build it here so the recipe is complete.
# -----------------------------------------------------------------------------
build_map_shmap() {
  local SHMAP_DIR="$REPO_ROOT/shmap"
  local SHMAP_BIN="$SHMAP_DIR/release/shmap"
  if [ -x "$SHMAP_BIN" ]; then
    echo "[skip] map-shmap already built at $SHMAP_BIN"
    return
  fi
  echo "[build] map-shmap (original shmap) ..."
  ( cd "$SHMAP_DIR" && make -j"$JOBS" )
  [ -x "$SHMAP_BIN" ] && echo "[ok] map-shmap -> $SHMAP_BIN" || echo "[warn] map-shmap build did not produce $SHMAP_BIN"
}

# ============================== run all ======================================
build_minimap2
build_winnowmap
build_blend
build_mapquik
build_map_shmap

echo
echo "================================ final check ================================"
for t in minimap2 winnowmap meryl blend mapquik; do
  printf '  %-12s %s\n' "$t" "$(command -v "$t" || echo MISSING)"
done
SHMAP_BIN="$REPO_ROOT/shmap/release/shmap"
printf '  %-12s %s\n' "map-shmap" "$([ -x "$SHMAP_BIN" ] && echo "$SHMAP_BIN" || echo MISSING)"
echo "============================================================================"
