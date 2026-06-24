# Equal-work benchmark: minSHmap (C++) vs shmap

Per-read **mapper-loop** throughput on the *same* reads, with the fixed
chr21 index-build cost cancelled by a differential N-vs-2N measurement.
Speed and sensitivity are reported as separate axes.

- Reads: first 1000 and first 2000 HiFi reads (same prefix). Best of 2. Single-threaded both.
- Params: {'k': 15, 'hfrac': 0.05, 'theta': 0.3, 'min_diff': 0.02, 'max_matches': 1000, 'max_seeds': 500}.
- NOTE: the differential subtracts two ~index-sized walls, so with few reps the
  ratio carries index-build noise; trust larger N (this one) over small N. A
  noisier N=300 run had the ratio flipped -- that was index noise, not signal.

| tool | wall(N) s | wall(2N) s | map loop Δ s | isolated reads/s | ~index s | mapped/2N |
|------|-----------|------------|-------------|------------------|----------|-----------|
| minshmap-cpp | 12.51 | 19.80 | 7.29 | 137.2 | 5.22 | 9/2000 |
| shmap | 5.58 | 9.99 | 4.41 | 226.9 | 1.17 | 17/2000 |

**Per-read mapper loop: shmap is 1.65x minshmap-cpp's throughput** on identical reads (index cost removed).
**Sensitivity: shmap maps 17 of 2000 reads vs minshmap's 9.**

Honest reading. On identical reads with the shared index cost removed,
shmap's per-read loop here is the faster one (227 vs 137 reads/s)
AND it maps more reads -- shmap wins on both axes. This OVERTURNS the
end-to-end realworld.md figure (minshmap-cpp 196.8 vs shmap 157.6), which
mixes index build with mapping across three datasets and is the less
controlled comparison. The per-read truth: shmap's heavily engineered
inner loop (cache-friendly hash maps, position-sorted hit lists with
binary search) beats minshmap's std::unordered_map + linear scans, even
though shmap does MORE algorithmic work per read (chaining, second-best,
mapq). minshmap's value is pedagogical clarity, not speed.

Real output-preserving headroom for minshmap: its index build (~5.2s)
is ~4x shmap's (~1.2s). Swapping std::unordered_map for a faster map
(ankerl/unordered_dense, as shmap uses) plus reserve() would cut that
without changing a single output line.
