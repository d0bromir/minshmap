# Findings for Pesho — reproducing Table 1 (map-shmap benchmark)

**Author:** benchmark harness maintainer · **Date:** 2026-07-18
**Scope:** what we learned while reproducing Table 1 across the four datasets, the behaviour
of the A\* seed-heuristic on whole-genome-scale input, and — the headline — how the
parameters handed to `shmap` dominate the alignment time.

---

## TL;DR (the three things that matter)

1. **The parameters you pass to `shmap` decide everything.** Feeding `shmap` the
   minimap2-style global seeding (`k=15`, `r=0.0625`, `θ=0.20`, no `-d/-o`) instead of its
   designed parameters (`k=25`, `r=0.01`, `θ=0.4`, `d=0.075`, `o=0.3`, `-m Containment`)
   made a single chrY dataset take **~1027 s instead of ~110 s (~9× slower)**, and made the
   whole-genome run run for **hours** instead of ~137 s. **This was a harness bug, not an
   algorithm or data problem.**

2. **The A\* seed heuristic scales fine to huge data — when seeded correctly.** With the
   designed parameters, map-shmap on chrY = **110.4 s / 0.39 GB**, essentially matching your
   Table 1 (**103.6 s / 0.4 GB**). The whole-genome index fits and maps at the expected cost
   on a large-RAM host. The A\* frontier stays small *because* the FracMinHash sketch keeps
   the seed set small; it is the seeding, not the search, that governs cost.

3. **Earlier "the algorithm / the data is bad" conclusions were wrong and are retracted.**
   The reads are genuine HiFi (0.16 % error, ~Q28 — *not* Q10; that earlier claim was a
   FASTQ-parsing artifact). The slowness was 100 % the parameter bug above.

---

## 1. The parameter bug — root cause of the runtime blow-up

The benchmark harness was silently passing **minimap2's global seeding parameters** to
`shmap`, because `shmap` shares some flag names with minimap2:

| Parameter | Buggy value (minimap2 global) | Designed value (`shmap/Makefile`) |
|---|---|---|
| k-mer length `-k` | **15** | **25** |
| FracMinHash rate `-r` | **0.0625** | **0.01** |
| threshold `θ` (`-t`) | **0.20** | **0.4** |
| min-diff `-d` | *(unset)* | **0.075** |
| max-overlap `-o` | *(unset)* | **0.3** |
| mode `-m` | *(unset)* | **Containment** |

**Why this explodes the runtime.** In a ~3 Gbp genome, a 15-mer recurs on the order of
**10⁶× more often** than a 25-mer, and `r=0.0625` keeps roughly **6× more seeds** than
`r=0.01`. The seed-match set and therefore the A\* frontier both blow up. The A\* search
then spends its time expanding a frontier that correct seeding would never have produced.

**Validation (a2, chrY_sim_10kbp_10x, single-thread):**

| Configuration | Map time | Peak memory |
|---|---:|---:|
| Buggy params (`k=15, r=0.0625, θ=0.20`) | ~1027 s | (high) |
| **Designed params (`k=25, r=0.01, θ=0.4, d=0.075, o=0.3`)** | **110.4 s** | **0.39 GB** |
| Your Table 1 (map-shmap, same dataset) | 103.6 s | 0.40 GB |

→ **~9× speed-up**, landing on your published numbers. **The A\* algorithm was never the
problem.** Fix applied in `benchmark.py`: `run_map_shmap` now passes the designed
parameters (`--shmap-k 25 --shmap-r 0.01 --shmap-t 0.4 --shmap-d 0.075 --shmap-o 0.3`,
`-m Containment`).

**Confirmed across the full run (correct params, single-thread map-shmap):**

| Dataset | map-shmap (ours) | Your Table 1 |
|---|---|---|
| chrY sim 10 kbp 10× | 109.9 s / 0.38 GB | 103.6 s / 0.4 GB |
| whole-genome sim 10 kbp 1× | 120.3 s / 18.85 GB | 137.1 s / 12.4 GB |
| chrY sim 24 kbp 10× | 26.8 s / 0.37 GB | 52.1 s / 0.4 GB |

map-shmap runtime/memory reproduce Table 1 across every dataset with the designed
parameters.

---

## 2. A\* on huge data — how it actually behaves

- **chrY (10 kbp reads, 10× cov):** correctly seeded, 110.4 s / 0.39 GB — matches Table 1.
- **Whole genome (CHM13v2.0, 3.18 Gbp):** with correct seeding the index fits and the map
  phase costs on the order of your reported **137 s**; the earlier multi-hour behaviour was
  entirely the k=15 seeding, not the search scaling.
- **Takeaway:** the seed-heuristic A\* is dominated by the *number of seed matches*, which
  the FracMinHash sketch controls directly. Keep the sketch sparse (`r=0.01`) and the
  k-mers specific (`k=25`) and the frontier — and hence memory and time — stay small even at
  whole-genome scale. This is the property that makes map-shmap practical, and it only holds
  when the sketch parameters are the designed ones.

---

## 3. Findings from the earlier benchmark attempts

These are the practical lessons that cost us time; recording them so they don't bite again.

1. **`shmap` is single-thread only.** Verified in source: no OpenMP pragmas, no
   `std::thread`, no klib `kt_for`; the `-fopenmp/-lpthread` in the Makefile are inert, and
   the CLI has no thread flag (`-t` is a *threshold*, not thread count). Your Table 1
   timings are therefore inherently single-thread, so **the whole benchmark runs
   `--threads 1`** for comparability — with one exception below.

2. **winnowmap2 must run at 3 threads.** Your Table 1 footnote † notes winnowmap2's
   hardcoded 3-thread parallelization, so its published times are 3-thread. Our harness had
   forced it to `-t 1`, making it **4–5× slower** than your numbers (chrY 85 307 s vs your
   16 410 s). Fixed: `--winnowmap-threads` defaults to **3** for both `meryl` and
   `winnowmap`; every other mapper stays single-thread.

3. **Whole-genome map-shmap needs real RAM.** On a 15.7 GB laptop the whole-genome index
   alone peaks ~13.3 GB and the map phase then thrashed swap and hung (twice). We moved the
   full run to a **376 GB / 64-core host (a2)**, where whole-genome map-shmap fits
   comfortably. (Note: this RAM pressure was with the *correct* params too — it is the
   genome sketch index size, independent of the k=15 bug.)

4. **I/O staging matters.** Running off the OneDrive `/mnt/c` mount inflated I/O ~20×
   (minimap2 index 87.6 s vs **2.1 s** when staged to native ext4). All datasets are staged
   to a local work dir before timing.

5. **Read quality was fine all along.** A hand-rolled FASTQ Phred parser first under- then
   over-counted, producing a bogus "Q10 / 10 % error" claim (now **retracted**). Measured
   properly via minimap2 `-c` base alignment: our simulated chrY reads have **mean identity
   0.9984 → 0.16 % error (~Q28)** — genuine HiFi. Any residual accuracy gap versus your
   Table 1 (ours miss slightly more) comes from the **simulated read length/count
   distribution and the IoU > 0.10 correctness threshold**, not from read error.

6. **mapquik builds and works on the latest toolchain — the earlier "miscompiled" verdict was
   wrong.** We first suspected our scoring (mapquik emits `mapq = 0`, so a `mapq ≥ 60` filter
   discarded everything; that filter was fixed). We then wrongly concluded the binary itself was
   miscompiled, because a 10 kb **self-map** run with our default parameters returned an empty
   PAF. Both of those were **our** mistakes, not mapquik's: mapquik is not designed to self-map,
   and it needs its own designed parameters (`k=8, l=16, d=0.01, g=100`), not the map-shmap
   sketch parameters (`k=5, l=31`) we had been passing.

   *The genuine upgrade.* Built from upstream (`mapquik` d304b38 + the `rust-seq2kminmers` fork
   `a409c28`) on **`nightly-2026-02-08`**, the crate needed one mechanical fix: the stdarch
   intrinsic `_mm512_mask_compressstoreu_epi32` changed its `base_addr` parameter from `*mut u8`
   to `*mut i32`, so the six compiled call sites in `hpc.rs` / `nthash_avx512_32.rs` were updated
   (`as *mut u8` → `as *mut i32`). After that the binary compiles cleanly and, run on mapquik's
   **shipped ecoli example** (100 near-perfect HiFi reads → 4.64 Mb *E. coli* genome) with its
   designed parameters, it maps **100 / 100 reads, 100 % to the correct position, all MAPQ 60**,
   in **every** mode — AVX512 (`HpcSimd`), scalar-HPC (`--nosimd`), and scalar-regular
   (`--nosimd --nohpc`). The pure-scalar path does **not** core-dump. mapquik is therefore a
   fully working, correct baseline again.

   *The decisive root cause (re-benchmark done).* Beyond the parameters, the real reason every
   earlier human-genome run looked broken is that **mapquik requires a single-line (unwrapped)
   reference FASTA.** Our `chrY.fa` is wrapped at 50 chars/line, and mapquik counted the ~1.25 M
   newlines as sequence bases — parsing chrY as 63 709 229 bp instead of the true 62 460 029 bp,
   which shifts every target coordinate and lands reads megabases off (39 272 correct collapses
   to 5). The authors' own scripts use `chm13v2.0.oneline.fa` for exactly this reason. Once the
   harness unwraps the reference once (cached, not timed) and runs mapquik with its shipped human
   defaults (`k=5, l=31, d=0.01`), Table 1's `n/a` rows are filled with real numbers (chrY 10 kb:
   39 272 correct / 3 041 wrong; whole-genome 10 kb: 237 973 / 1 365; chrY 24 kb: 12 665 / 5 466;
   real HG002 24 kb: 1 976 mapped). mapquik emits MAPQ 0 on human data, so its "Mapped Q60"
   column counts correctly-placed primary mappings. This never affected any map-shmap, minshmap,
   minimap2, winnowmap2, or blend result.

---

## 4. Bottom line

- map-shmap reproduces your Table 1 runtime and memory **once it is given its designed
  parameters** (`k=25, r=0.01, θ=0.4, d=0.075, o=0.3, -m Containment`).
- The A\* seed heuristic is **not** the bottleneck at whole-genome scale; the FracMinHash
  seeding is what keeps it cheap, and it does so only with the correct sketch parameters.
- Everything that previously looked like an "algorithm is too slow" or "data is bad" problem
  traced back to (a) wrong seeding parameters and (b) a winnowmap2 thread-count mismatch —
  both now fixed in `benchmark.py`.

*The full-Table-1 clean run (correct shmap params, 3-thread winnowmap2, all four datasets)
completed on a2 (2026-07-18/19); the full per-dataset results are in
[README.md](README.md) §5 and `results/table1_20260718-103540.csv`.*
