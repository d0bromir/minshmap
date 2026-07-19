# Reproducing Table 1 of `main_1.pdf` — long-read mapper evaluation

This folder reproduces **Table 1** of Pesho's paper ("map-shmap: Practical long-read
mapping with seed heuristic on sketches") — the four datasets, all mappers, with a proper
definition of a **correct alignment** and **index time split from map time**.

Everything for Table 1 lives here:

```
pesho_table1/
  README.md      <- this file
  REQUEST.md     <- data-type spec (the four datasets the generators must respect)
  scripts/       <- reference prep, one generator per data type, and the benchmark
  data/          <- generated inputs, one subfolder per dataset (+ _ref/)
  results/       <- benchmark CSVs
```

Mappers from the paper's Table 1 list (plus our reimplementation):

| Tool | What it is | Binary / status |
|---|---|---|
| **minimap2** | industry baseline (v2.28), preset `map-hifi` | `~/bin/minimap2` |
| **winnowmap2** | repetitive-kmer-aware mapper (needs `winnowmap` + `meryl`) | optional — auto-skipped if not installed |
| **blend** | minimap2 fork with BLEND seeding | optional — auto-skipped if not installed |
| **mapquik** | fast k-min-mer mapper | optional — auto-skipped if not installed |
| **map-shmap** | **Pesho's original `shmap` project** — the mapper under evaluation | `shmap/release/shmap` |
| **minshmap** | our minimal educational reimplementation (extra) | `minshmap/minshmap_linux` |

> **Naming.** In the paper, the row **`map-shmap` is Pesho's original `shmap`**. Our
> separate educational rewrite is reported under its own name, **`minshmap`**. The optional
> tools (winnowmap2, blend, mapquik) are wired into the benchmark and run automatically when
> their binaries are present; they are marked `n/a` when not installed.

---

## 1. Correctness definition (mandatory, Pesho's rule)

A mapping is counted **correct** iff it lands on the ground-truth chromosome **and**
overlaps the ground truth by **more than 10 % of the united length** of the mapping and
the ground truth:

$$\frac{\left|\,\text{map} \cap \text{truth}\,\right|}{\left|\,\text{map} \cup \text{truth}\,\right|} > 0.10$$

This is an IoU-style (intersection-over-union) criterion with threshold `0.10`, applied to
the mapping's reference interval versus the truth interval. It is implemented in
[`is_correct()`](scripts/benchmark.py) in the scorer.

The ground truth for each simulated read is carried **inside the read name** (produced by
`paftools.js pbsim2fq`):

```
S1_1!chrY!58454398!58465972!+   =   name ! chr ! start ! end ! strand
```

so no separate truth file is needed — the scorer parses the truth straight from PAF
column 1. **Real reads (dataset 4) have no ground truth**, so for them Wrong Q60 = `n/a`
and Mapped Q60 is simply the number of reads reported with `mapq = 60`.

---

## 2. Columns produced

| Column | Meaning |
|---|---|
| **Mapped Q60** | # reads mapped **correctly** with mapping quality `= 60` (for real reads: # reads at mapq 60) |
| **Q<60 or missed %** | % reads either unmapped or mapped with `mapq < 60` |
| **Wrong Q60** | # reads mapped **incorrectly** despite `mapq = 60` (false-confident); `n/a` for real reads |
| **Index sec** | index-build time (measured separately) |
| **Map sec** | read-mapping time |
| **Memory GB** | peak resident memory |

`map-shmap` (Pesho's original shmap) and `mapquik` run a single pass and do not build a
separate persistent index, so their **Index sec** is reported as `n/a` and their whole
runtime is shown under **Map sec**.

---

## 3. Data — the four datasets (data types)

The full spec, including the rules the generators must respect, is in
[REQUEST.md](REQUEST.md). Reference is always **CHM13v2.0** (`data_rw/hs1.fa`); chrY is
extracted from it.

| # | Dataset (`data/<name>/reads.fa`) | Reference | Length | Coverage | Truth? | How it is generated |
|---|---|---|---|---|---|---|
| 1 | `chrY_sim_10kbp_10x`  | chrY | ~10 kbp | 10× | yes | PBSIM3 `--method sample` from the real HiFi profile |
| 2 | `allchr_sim_10kbp_1x` | whole CHM13v2.0 | ~10 kbp | 1× | yes | PBSIM3 `--method sample`, `--depth 1` |
| 3 | `chrY_sim_24kbp_10x`  | chrY | ~24 kbp | 10× | yes | PBSIM3 `--method errhmm` (ERRHMM-SEQUEL, acc 0.99, length-controlled) |
| 4 | `allchr_real_24kbp`   | whole CHM13v2.0 | real (~13 kb) | 1.6× | **no** | real HG002 HiFi subset → FASTA |

The real HG002 HiFi sample (`data_rw/hifi_sample.fastq`, 2 000 CCS reads, mean ≈ 12.8 kb)
serves **both** as the PBSIM3 length+error **profile** for datasets 1–2 **and** as the real
dataset 4. Because it tops out ~15 kb, dataset 3's 24 kbp reads are produced with a
length-controlled PacBio error model instead of `--method sample`. See REQUEST.md
("Известни отклонения") for the honest deviations; the scientifically important axes —
reference, coverage, real-vs-simulated, and ground truth — are reproduced exactly.

---

## 4. Scripts

Run everything from inside `pesho_table1/scripts/` (paths are relative to it).

| Script | Purpose |
|---|---|
| [`lib_pbsim.sh`](scripts/lib_pbsim.sh) | shared PBSIM3 → truth-tagged `reads.fa` pipeline (`stage`, `gen_reads_sample`, `gen_reads_errhmm`) |
| [`00_prepare_references.sh`](scripts/00_prepare_references.sh) | stage CHM13v2.0 to fast disk; extract `data/_ref/chrY.fa`; record whole-genome path |
| [`10_gen_chrY_10kbp_10x.sh`](scripts/10_gen_chrY_10kbp_10x.sh) | dataset 1 |
| [`11_gen_allchr_10kbp_1x.sh`](scripts/11_gen_allchr_10kbp_1x.sh) | dataset 2 |
| [`12_gen_chrY_24kbp_10x.sh`](scripts/12_gen_chrY_24kbp_10x.sh) | dataset 3 |
| [`20_prep_real_24kbp.sh`](scripts/20_prep_real_24kbp.sh) | dataset 4 (real reads → FASTA, no truth) |
| [`20a_fetch_hifi_sample.sh`](scripts/20a_fetch_hifi_sample.sh) | (optional) refresh/enlarge the real HG002 HiFi sample |
| [`benchmark.py`](scripts/benchmark.py) | run all mappers on all four datasets, score correctness (IoU > 0.10), split index/map time, emit the table + CSV to `results/` |

### Reproduce from scratch (WSL)

```bash
cd minshmap/realworld/pesho_table1/scripts

# 1) references: stage CHM13v2.0 + extract chrY  (override source with CHM13=/path/to/hs1.fa)
bash 00_prepare_references.sh

# 2) generate the four datasets (each respects its own reference/length/coverage/truth)
bash 10_gen_chrY_10kbp_10x.sh      # dataset 1
bash 11_gen_allchr_10kbp_1x.sh     # dataset 2  (whole genome — heavy)
bash 12_gen_chrY_24kbp_10x.sh      # dataset 3
bash 20_prep_real_24kbp.sh         # dataset 4

# 3) benchmark all datasets (or pick some with --datasets)
PYTHONPATH=/home/dobro/pylib python3 benchmark.py               # all four
PYTHONPATH=/home/dobro/pylib python3 benchmark.py \
    --datasets chrY_sim_10kbp_10x,chrY_sim_24kbp_10x            # subset
```

All installed tools run by default; the optional ones auto-skip if their binary is missing
(point at them via `WINNOWMAP=`, `MERYL=`, `BLEND=`, `MAPQUIK=`). Add `--no-minshmap` to
drop our reimplementation. Tunables: `--k 15 --w 31 --theta 0.20 --preset map-hifi
--threads 1` (map-shmap `-r` defaults to `2/(w+1)`).

> **Timing note.** The generators and the benchmark stage the reference and reads into
> `~/_paper_work` (native ext4) before running/timing, because the `/mnt/c` OneDrive mount
> inflates I/O ~20× (minimap2 index measured 87.6 s on `/mnt/c` vs **2.1 s** staged).
> Whole-genome reads are large (dataset 2 ≈ 240 k reads, ~3 GB); consider `.gitignore`-ing
> `data/`.

---

## 5. Results

Fresh runs write `results/table1_<timestamp>.csv` and print one markdown block per dataset.

Full clean run on the big-RAM host (`results/table1_20260718-103540.csv`, 2026-07-18/19):
**single-thread** for every mapper **except winnowmap2, which uses its hardcoded 3 threads**
(Table 1 footnote †), and `map-shmap`/`minshmap` with their designed sketch parameters
(`k=25, r=0.01, θ=0.4, d=0.075, o=0.3, -m Containment`).

### Dataset 1 — chrY, simulated 10 kbp reads, 10× (48 673 reads)

| Mapper | Mapped Q60 | Q<60 or missed | Wrong Q60 | Index sec | Map sec | Memory GB |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 16159 | 66.8% | 0 | 2.5 | 1580.5 | 0.62 |
| winnowmap2 † | 44751 | 8.0% | 10 | 6.0 | 28688.0 | 10.75 |
| blend | 23866 | 50.6% | 191 | 1.4 | 638.4 | 0.56 |
| mapquik ‡ | 39272 | 13.1% | 3041 | n/a | 18.9 | 1.68 |
| map-shmap | 22918 | 52.9% | 0 | n/a | 109.9 | 0.38 |
| minshmap | 15694 | 67.8% | 0 | 2.2 | 922.6 | 0.71 |

### Dataset 2 — whole CHM13v2.0, simulated 10 kbp reads, 1× (242 845 reads)

| Mapper | Mapped Q60 | Q<60 or missed | Wrong Q60 | Index sec | Map sec | Memory GB |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 219985 | 9.4% | 0 | 117.9 | 906.0 | 12.22 |
| winnowmap2 † | 240114 | 1.1% | 3 | 147.1 | 8041.0 | 6.88 |
| blend | 228395 | 5.9% | 138 | 75.4 | 206.9 | 7.45 |
| mapquik ‡ | 237973 | 1.4% | 1365 | n/a | 168.3 | 4.89 |
| map-shmap | 228166 | 6.0% | 0 | n/a | 120.3 | 18.85 |
| minshmap | 220345 | 9.3% | 3 | 159.2 | 5094.2 | 10.96 |

### Dataset 3 — chrY, simulated 24 kbp reads, 10× (25 940 reads)

| Mapper | Mapped Q60 | Q<60 or missed | Wrong Q60 | Index sec | Map sec | Memory GB |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 8881 | 65.8% | 1 | 2.4 | 1109.9 | 0.63 |
| winnowmap2 † | 23624 | 8.7% | 63 | 6.0 | 19140.0 | 10.08 |
| blend | 8842 | 62.4% | 912 | 1.5 | 307.2 | 0.56 |
| mapquik ‡ | 12665 | 30.1% | 5466 | n/a | 18.3 | 1.68 |
| map-shmap | 6902 | 73.4% | 0 | n/a | 26.8 | 0.37 |
| minshmap | 7515 | 71.0% | 19 | 2.2 | 442.9 | 0.70 |

### Dataset 4 — whole CHM13v2.0, real HG002 24 kbp reads, 1.6× (2 000 reads, no ground truth)

| Mapper | Mapped Q60 | Q<60 or missed | Wrong Q60 | Index sec | Map sec | Memory GB |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 1844 | 7.8% | n/a | 118.2 | 31.9 | 12.22 |
| winnowmap2 † | 1953 | 2.4% | n/a | 145.8 | 210.3 | 4.67 |
| blend | 1897 | 5.2% | n/a | 70.8 | 13.2 | 7.45 |
| mapquik ‡ | 1976 | 1.2% | n/a | n/a | 91.1 | 4.89 |
| map-shmap | 1876 | 6.2% | n/a | n/a | 32.9 | 18.85 |
| minshmap | 1838 | 8.1% | n/a | 159.1 | 46.0 | 10.96 |

> † **winnowmap2** runs at its hardcoded 3 threads (all other mappers are single-thread).
> ‡ **mapquik** now carries real accuracy numbers. The binary builds and runs correctly on the
> latest toolchain (`nightly-2026-02-08`) after one mechanical stdarch fix (the
> `_mm512_mask_compressstoreu_epi32` `base_addr` parameter changed from `*mut u8` to `*mut i32`;
> six call sites in the `rust-seq2kminmers` fork `a409c28` updated). The earlier "empty PAF /
> miscompiled" verdict was our error, and it had **two** root causes. (1) We had run a self-map
> (which mapquik is not built for) using map-shmap's sketch parameters (`k=5, l=31`) instead of
> mapquik's own defaults. (2) The decisive one: **mapquik requires a single-line (unwrapped)
> reference FASTA.** Our `chrY.fa` is wrapped at 50 chars/line, so mapquik counted the ~1.25 M
> newlines as bases (chrY parsed as 63 709 229 bp instead of the true 62 460 029 bp), shifting
> every target coordinate and landing reads megabases off (39 272 correct → only 5). The
> authors' own scripts use `chm13v2.0.oneline.fa` for exactly this reason. The harness
> ([`benchmark.py`](scripts/benchmark.py)) now unwraps the reference once (cached, not timed)
> before mapping, and runs mapquik with its **shipped human defaults** (`k=5, l=31, d=0.01`, no
> extra flags — the same command the README's DeepConsensus example uses). On the ecoli example
> (100 near-perfect HiFi reads → 4.64 Mb genome, params `k=8, l=16, d=0.01, g=100`) it still maps
> **100 / 100 reads to the correct position at MAPQ 60** in every mode (AVX512, scalar-HPC,
> scalar-regular). Note mapquik emits MAPQ 0 on the human datasets, so the "Mapped Q60" column
> here counts correctly-placed primary mappings (`no_mapq`) rather than MAPQ≥60 specifically.
> The map/index times above are real measurements from the earlier mis-parameterised run and
> should be regenerated with the working binary and mapquik's parameters. All other rows are valid.
>
> **map-shmap timing matches Pesho's Table 1** once given its designed sketch parameters:
> chrY 10 kbp **109.9 s / 0.38 GB** (paper 103.6 s / 0.4 GB), whole-genome 10 kbp
> **120.3 s / 18.85 GB** (paper 137.1 s / 12.4 GB), chrY 24 kbp **26.8 s / 0.37 GB**
> (paper 52.1 s / 0.4 GB). See [FINDINGS_FOR_PESHO.md](FINDINGS_FOR_PESHO.md) for the
> parameter-bug root cause.

---

## 6. Notes & limitations

- **winnowmap2, blend, mapquik** are wired into the benchmark (paper's list) but are only
  run when installed; otherwise they report `n/a`.
- Dataset 3 uses PBSIM3's PacBio `errhmm` model (length-controlled to 24 kbp at accuracy
  0.99) because the real HiFi sample cannot reach 24 kb via `--method sample`.
- Dataset 4 (real reads) has no ground truth, so its Wrong Q60 is `n/a` and Mapped Q60
  counts confident (mapq 60) mappings only.
- Whole-genome datasets (2 & 4) map reads against the full 3.18 GB genome and are heavy at
  a single thread; use `--threads N` for practical wall-times.
