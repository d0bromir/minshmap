# minSHmap

A minimalistic, **educational** sketch-based long-read mapper — a small, readable
re-implementation of [`shmap`](../shmap) in the spirit of
[`minSH`](../minSH). The goal is to be as simple to *understand* as possible,
not as fast as possible. The same algorithm is provided twice: a clear pure
[Python implementation](minshmap.py) and a faithful
[C++ port](minshmap.cpp) for the speed comparison.

A read **mapper** answers *"where in the reference does this read come from?"* —
it finds a **location**, not a base-by-base alignment. That is the one job
minSHmap does.

---

## 1. The minSHmap algorithm

minSHmap maps reads in three steps: **sketch → index → map**.

### Step 1 — Sketch (FracMinHash)

Each sequence is reduced to a small, strand-canonical set of k-mers:

```
sketch(seq, k, hfrac) -> [(pos, hash, strand), ...]
```

- Every k-mer gets a 64-bit hash. We keep a k-mer only if its hash falls in the
  smallest fraction `hfrac` of the hash space (**FracMinHash**). This keeps a
  roughly uniform ~`hfrac` sample of k-mers and lets two sequences be compared
  by the k-mers they share.
- The hash is **canonical**: a k-mer and its reverse complement hash to the same
  value (`h = min(h_fw, h_rc)`), so a read still matches when it was sequenced
  from the opposite strand. `strand` records which orientation won.

Three interchangeable sketchers are provided to study the speed/clarity
trade-off (selectable with `--hash`):

| sketcher | hashing | cost |
| --- | --- | --- |
| `naive` | re-hash every k-mer string (blake2b) | O(len · k) |
| `poly` | Rabin–Karp polynomial rolling hash | O(len) |
| `nthash` | rotation-based rolling hash (mirrors shmap) | O(len) |

### Step 2 — Index the reference

Sketch every reference segment and remember where each k-mer occurs:

```
index:  hash -> [(segment_id, position, strand), ...]
```

### Step 3 — Map a read

1. **Sketch the read** into `m` k-mers (its "seeds").
2. **Rank seeds rarest-first** by how many times each occurs in the reference
   index. Rare k-mers are the strongest, cheapest evidence of a true location.
3. **Seed candidate blocks.** Scatter the hits of the rarest seeds into
   read-length blocks over the reference. Consecutive blocks
   overlap (each hit is added to block `b` and `b-1`) so any homologous region
   lands fully inside at least one block.
4. **Refine each block with the seed heuristic (the "SH" in shmap).** Add the
   read's seeds one by one, rarest first, and track

   $$ sh = 1 - \frac{seeds\_used - matches}{m} $$

   `sh` is an **upper bound** on the containment the block can still reach. The
   moment `sh` drops below the homology threshold `theta`, the block can never
   be good enough, so it is **pruned immediately** without scanning the rest.
5. **Score and report.** A surviving block's score is its containment
   `matches / m`. Keep the best block (and a second-best from a different
   region to derive a minimap2-style `mapq`), then emit a minimal **PAF** line.

The key idea — and the thread linking this project to `minSH` — is the **seed
heuristic**: counting seeds that *cannot* be matched gives a cheap lower bound on
the work remaining. `minSH` uses it to bound remaining **edit distance** and
guide A\* alignment; minSHmap re-uses the very same counting to bound achievable
**containment** and prune candidate **locations** during mapping.

---

## 2. The shmap algorithm

[`shmap`](../shmap) is the production C++ mapper minSHmap is distilled from. It
follows the same sketch → index → map skeleton, but each stage is hardened for
real genomes and speed:

- **Sketch** — FracMinHash with a true ntHash rolling hash (`LUT_fw`/`LUT_rc`
  with `rotl`/`rotr`), producing `Kmer{r, h, strand}`.
- **Index** (`SketchIndex`) — splits hashes into `h2single` (unique) and
  `h2multi` (repeated, sorted by position for binary search), with optional
  `erase_frequent_kmers` to drop satellite-repeat k-mers.
- **Map** (`shmap.h::map_read`):
  - group read k-mers, record per-seed occurrences and hit counts, sort seeds
    rarest-first;
  - use `S = (1 - theta2)·m + 1` seeds to scatter hits into **overlapping
    buckets**;
  - `match_rest` applies the **seed-heuristic prune** lazily per bucket, with a
    threshold that *ratchets upward* toward the best score seen so far;
  - choose a similarity **metric** — default **Containment** (also Jaccard,
    bucket-SH, and LCS via a longest-increasing-subsequence of matched
    positions);
  - keep the best **and** second-best window and compute a minimap2-like
    **mapq** from their ratio;
  - emit a full PAF line (query/target coordinates, strand, matches, block
    length, mapq).

shmap additionally carries production machinery that minSHmap deliberately
omits: templated pruning modes, co-linear chaining safeguards, profiling
(Tracy), counters/timers, and multiple output metrics.

---

## 3. Differences between minSHmap and shmap

| aspect | minSHmap | shmap |
| --- | --- | --- |
| **purpose** | teaching: minimal, readable | production mapper |
| **size** | one ~400-line file (×2 for the C++ port) | many C++ headers/translation units |
| **sketch** | 3 swappable sketchers (`naive`/`poly`/`nthash`) | one tuned ntHash FracMinHash |
| **index** | single `hash -> list` dict | `h2single` + `h2multi`, binary search, frequent-k-mer erase |
| **similarity metric** | Containment only | Containment / Jaccard / bucket-SH / LCS |
| **prune threshold** | fixed `theta` | `theta` that **ratchets up** toward the best score |
| **acceptance rule** | exact containment `≥ theta` | mapq-filtered, with co-linear **chaining** |
| **mapq** | best vs second-best ratio (simple) | minimap2-style, chaining-aware |
| **engineering** | none (stdlib only in Python) | Tracy profiling, templated modes, counters/timers |
| **output** | minimal PAF | full PAF |

The most consequential difference for results is the **acceptance rule**:
minSHmap accepts a location purely on a single containment threshold, whereas
shmap accepts on a mapq computed from chained, second-best-aware evidence. That
extra machinery is exactly what gives shmap its higher **sensitivity** (it maps
more reads) — and reproducing it is what would push minSHmap past the simplicity
budget that is the whole point of the project (see the honest note below).

---

## 4. Benchmark — latest results and analysis

> **Implementation note.** minSHmap has since been unified on a *single*
> strand-canonical `(w, k)`-minimizer sketcher from the
> [`minimizer-iter`](minimizer_ext/) library — it replaces the earlier
> `naive`/`poly`/`nthash` FracMinHash sketchers described in §1. Both the
> [Python](minshmap.py) and the [C++](minshmap.cpp) port call that same library and
> share the same `bisect`/binary-search hit lookup, so they are **algorithmically
> equivalent** and emit **byte-identical** PAF. Mapping quality now uses the
> parameter-free ("φ-free") rule — see
> [NOTE_phi_elimination.md](NOTE_phi_elimination.md). The benchmark therefore has one
> minSHmap participant, run as both Python and C++, with `shmap` as the production
> reference.

### Real-world long-read benchmark — T2T-CHM13v2.0 chr21 (45,090,682 bp)

2,000 reads per dataset, `k=15, w=11`. Reads are whole-genome WGS, so a **low mapped
fraction is expected** — only the ~1–2% originating from chr21 should map. `shmap`
figures are its published 2026-06-21 baseline (its own settings), shown for reference.
Full data: [realworld.md](realworld/results_rw/realworld.md) and the latest CSV in
[results_rw/](realworld/results_rw/).

| dataset | py reads/s | cpp reads/s | shmap reads/s | py mapped | cpp mapped | shmap mapped |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hifi | 7.0 | 395.6 | 271.0 | 15 | 15 | 17 |
| ont | 1.7 | 154.3 | 69.0 | 32 | 32 | 32 |
| clr | 12.6 | 430.4 | 1142.0 | 2 | 2 | 2 |

### Synthetic short-read benchmark

8,000 reads (~300 bp, 5 % error), `k=15, w=11, theta=0.5`:

| implementation | reads/s | mapped | placement precision |
| --- | ---: | ---: | ---: |
| minSHmap-py | 14,499 | 3385 | 100 % |
| minSHmap-cpp | 50,084 | 3385 | — |

### Analysis

- **Python and C++ are algorithmically identical.** They map the *same* read count on
  every dataset (synthetic 3385 / 3385; chr21 15/15, 32/32, 2/2) and emit byte-identical
  PAF. The pure-Python mapper is for reading, not for production speed.
- **C++ now beats shmap on throughput** for the long-read sets — HiFi 395.6 vs 271 r/s
  and ONT 154.3 vs 69 r/s — trailing only on the short CLR reads (430 vs 1142 r/s), where
  shmap's lower per-read constant wins.
- **Sensitivity matches shmap** on ONT (32 = 32) and CLR (2 = 2), and is close on HiFi
  (15 vs 17). On the repetitive acrocentric chr21 the educational mapper is competitive
  with the production tool.
- **Why pure Python is slow on ONT (1.7 r/s).** Not an algorithm bug — the placements
  equal C++. ONT reads are ~35 kb, so each carries thousands of minimizers; the cost is the
  per-read minimizer work plus the 43 MB chr21 index build, times the ~100× CPython
  interpreter constant. The hit-list lookup itself is now `O(log n)` binary search in *both*
  implementations (see [NOTE_phi_elimination.md](NOTE_phi_elimination.md) §6, which also
  documents how this removed an earlier `O(hits²)` blow-up on repetitive references).

### Honest note on the remaining sensitivity gap

minSHmap still maps slightly fewer HiFi reads than shmap (15 vs 17). The missed reads sit
right at the containment threshold — 16 kb partial/junction reads only ~30 % chr21, so
whole-read containment is inherently diluted. Closing the gap *properly* needs shmap's
mapq filtering and co-linear chaining, which exceed the "keep it simple, keep it readable"
budget that is the entire point of this project. So the gap is left **visible** rather than
papered over with a precision-losing hack. minSHmap stays simple and honest; shmap stays
the sensitive production tool.

---

## Running it

```bash
# Python — map reads, print PAF
python minshmap.py reference.fa reads.fa -k 15 -w 11 -t 0.9

# C++ port (faster) — byte-identical output
g++ -O3 -std=c++17 -march=native -pthread \
    -I ../shmap/ext/unordered_dense/include \
    -o minshmap minshmap.cpp -L minimizer_ext/target/release \
    -l:libminimizer_ext.a -lws2_32 -luserenv -lbcrypt -lntdll
./minshmap reference.fa reads.fa -k 15 -w 11 -t 0.9
```

`-w` must be odd (canonical minimizers); the Windows link libraries
(`-lws2_32 -luserenv -lbcrypt -lntdll`) are not needed on Linux. The real-world benchmark
harness lives in [realworld/](realworld); `python realworld/09_bench_py_vs_cpp.py chr21`
races the Python and C++ implementations and writes a timestamped CSV (with `shmap`
reference columns) into [results_rw/](realworld/results_rw/).
