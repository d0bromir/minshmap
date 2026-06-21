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
3. **Seed candidate windows.** Scatter the hits of the rarest seeds into
   read-length windows ("buckets") over the reference. Consecutive windows
   overlap (each hit is added to bucket `b` and `b-1`) so any homologous region
   lands fully inside at least one window.
4. **Refine each window with the seed heuristic (the "SH" in shmap).** Add the
   read's seeds one by one, rarest first, and track

   $$ sh = 1 - \frac{seeds\_used - matches}{m} $$

   `sh` is an **upper bound** on the containment the window can still reach. The
   moment `sh` drops below the homology threshold `theta`, the window can never
   be good enough, so it is **pruned immediately** without scanning the rest.
5. **Score and report.** A surviving window's score is its containment
   `matches / m`. Keep the best window (and a second-best from a different
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

## 4. Benchmark — final results and analysis

Real-world long-read benchmark against **T2T-CHM13v2.0 chr21** (45,090,682 bp),
2,000 reads per dataset, unified parameters
`k=15, hfrac=0.05, theta=0.3, min_diff=0.02, max_matches=1000`. Reads are
whole-genome WGS, so a **low mapped fraction is expected** — only the ~1–2% of
reads originating from chr21 should map at all. Full data:
[realworld.md](realworld/results_rw/realworld.md) and
[realworld.csv](realworld/results_rw/realworld.csv).

### Speed ranking (overall reads/s across all datasets)

| rank | tool | reads/s |
| --- | --- | --- |
| 1 | `minshmap-cpp-nthash` | 196.8 |
| 2 | `minshmap-cpp-poly` | 171.7 |
| 3 | `shmap` | 157.6 |
| 4 | `minshmap-cpp-naive` | 146.1 |
| 5 | `minshmap-py-nthash` | 2.5 |
| 6 | `minsh` | 0.09 |

### Accuracy ranking (chr21 reads recovered)

| rank | tool | total mapped | mean agree w/ shmap | identity |
| --- | --- | --- | --- | --- |
| 1 | `shmap` | 51 | — (reference) | — |
| 2 | `minshmap-cpp-naive` | 22 | 0.197 | — |
| 3 | `minshmap-cpp-poly` | 20 | 0.187 | — |
| 4 | `minshmap-py-nthash` | 17 | 0.157 | — |
| 4 | `minshmap-cpp-nthash` | 17 | 0.157 | — |
| 6 | `minsh` | 1 | — | 0.869 |

Per-dataset mapped counts: **hifi** shmap 17 / cpp-naive 11 / py-nthash 9;
**ont** shmap 32 / cpp-naive 11 / py-nthash 8; **clr** shmap 2 / all minshmap 0.

### Brief analysis

- **C++ is ~70–190× faster than pure Python**, as expected for an educational
  implementation (~2 reads/s vs ~150–390 reads/s). The Python mapper is for
  reading, not for production.
- **ntHash is the fastest C++ sketcher on every dataset** — its rotation-based
  rolling hash maps directly onto native 64-bit rotate instructions. This is
  precisely shmap's design rationale, now confirmed on real data, which is why
  the Python participant is run with ntHash.
- **Speed-optimal ≠ sensitivity-optimal.** `naive` and `poly` map *marginally
  more* reads than `nthash` (22/20 vs 17 total) because each hash function
  samples a slightly different subset of k-mers. The differences are small and
  within sketching noise.
- **shmap is the sensitivity leader** (51 reads vs minSHmap's 17–22) — see the
  note below for why.

### Detailed analysis

- **Why minSHmap maps fewer reads than shmap.** The gap is an *acceptance-rule*
  problem, not a candidate-generation one. On **HiFi**, the missed reads sit at
  containment ≈ 0.289–0.304, right at the `theta=0.3` cliff: these 16 kb reads
  are only ~30% chr21 (partial/junction reads), so whole-read containment is
  inherently diluted. On **ONT**, true-locus containment is only ≈ 0.06 because
  k-mer survival under ~13% error is `(1-e)^k`, so *no* `theta=0.3` containment
  rule could ever map them. shmap clears both because its mapq + chaining accept
  on locally-supported, chained evidence rather than a single global threshold.
- **Concordance is low and that is partly expected.** On the repetitive
  acrocentric chr21, `agree_with_shmap` (interval overlap + strand) understates
  agreement: minSHmap and shmap often place the same read at *different copies*
  of the same repeat. On ONT, manual k-mer-containment adjudication of the
  shared reads gave minSHmap 5 wins vs shmap 3 — neither is "wrong", they pick
  different repeat instances.
- **CLR is too short/noisy for these parameters.** At avg 1.8 kb with ~10–15%
  error, k-mer survival at `k=15` is near zero; minSHmap maps 0 and even shmap
  only 2. This is a parameter regime, not a bug.
- **minSH as an independent ground-truth check.** minSH is an A\* *aligner*, not
  a mapper. It aligned one HiFi placement to its chr21 window at **identity
  0.869** (≫ 0.25 random), confirming minSHmap's placements are real. It can
  only validate a single read here because A\* on 16 kb reads at ~13% divergence
  degenerates toward Dijkstra (~3 M cells/read) — it is a spot-check, not a
  throughput metric, and is listed last in the speed ranking for exactly that
  reason.

### Honest note on closing the sensitivity gap

Matching shmap's mapped count (51 vs minSHmap's ~17–22) was **investigated and
deliberately reverted**. A simple relaxation — gating/scoring on shmap's
seed-heuristic lower bound instead of exact containment — lifted HiFi from 9 to
15 of shmap's 17, but introduced **~2% false positives** (13 spurious placements
on 600 reads). Closing the gap *properly* requires shmap's mapq filtering and
co-linear chaining safeguards, which exceed the "keep it simple, keep it
readable" constraint that is the entire purpose of this project. So the gap is
left **visible** in the accuracy ranking rather than papered over with a
precision-losing hack. minSHmap stays simple and honest; shmap stays the
sensitive production tool.

---

## Running it

```bash
# map reads, print PAF
python minshmap.py reference.fa reads.fa

# generate data and self-test
python minshmap.py --demo

# choose a sketcher
python minshmap.py reference.fa reads.fa --hash nthash

# C++ port (faster)
g++ -O3 -std=c++17 -march=native -o minshmap_linux minshmap.cpp
./minshmap_linux reference.fa reads.fa --hash nthash
```

The real-world benchmark harness lives in [realworld/](realworld); run
`realworld/05_run_benchmark.py` (it builds the C++ binary, downloads chr21 and
read sets, and writes `realworld.csv` / `realworld.md`).
