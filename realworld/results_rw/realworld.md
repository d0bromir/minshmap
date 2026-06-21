# minSHmap real-world long-read benchmark

- generated: 2026-06-21T18:20:15+00:00
- reference: T2T-CHM13v2.0 chr21 (45,090,682 bp)
- params: {'k': 15, 'hfrac': 0.05, 'theta': 0.3, 'min_diff': 0.02, 'max_matches': 1000, 'max_seeds': 500}
- python index build (one-time): {'nthash': 63.08} s
- python: 3.12.3

Participants: **minshmap-py-nthash** (educational pure-Python mapper, run with ntHash — the sketcher that won the C++ speed comparison), **minshmap-cpp-{naive,poly,nthash}** (the C++ port with each sketcher), **shmap** (the reference C++ mapper), and **minsh** (the A* aligner, used to independently score placement quality).

Reads are whole-genome WGS, so a low mapped fraction is expected (only the ~1-2% from chr21 should map). `agree_with_shmap` = fraction of shmap's chr21 placements that the tool reproduces at an overlapping locus and strand. `minsh` is an aligner, not a mapper: its `mapped` is the number of placements it aligned and its accuracy is alignment `identity`, not concordance.

## hifi — PacBio HiFi

- 2000 reads, 32,519,382 bp, avg 16,259 bp/read

| tool | mapped | mapped_frac | wall s | reads/s | agree w/ shmap | identity |
| --- | --- | --- | --- | --- | --- | --- |
| minshmap-py-nthash | 9 | 0.0045 | 1163.772 | 2 | 0.471 |  |
| minshmap-cpp-naive | 11 | 0.0055 | 14.412 | 139 | 0.529 |  |
| minshmap-cpp-poly | 10 | 0.005 | 13.396 | 149 | 0.529 |  |
| minshmap-cpp-nthash | 9 | 0.0045 | 11.639 | 172 | 0.471 |  |
| shmap | 17 | 0.0085 | 7.383 | 271 |  |  |
| minsh | 1 |  | 11.693 | 0.086 |  | 0.8691 |

minSH aligned 1 HiFi placements to their chr21 window: median identity **0.8691**, median A* cells 3027482 (in 11.693 s).

## ont — Oxford Nanopore

- 2000 reads, 70,714,827 bp, avg 35,357 bp/read

| tool | mapped | mapped_frac | wall s | reads/s | agree w/ shmap | identity |
| --- | --- | --- | --- | --- | --- | --- |
| minshmap-py-nthash | 8 | 0.004 | 1022.893 | 2 | 0.0 |  |
| minshmap-cpp-naive | 11 | 0.0055 | 19.568 | 102 | 0.062 |  |
| minshmap-cpp-poly | 10 | 0.005 | 15.415 | 130 | 0.031 |  |
| minshmap-cpp-nthash | 8 | 0.004 | 13.691 | 146 | 0.0 |  |
| shmap | 32 | 0.016 | 28.932 | 69 |  |  |

## clr — PacBio CLR

- 2000 reads, 3,671,813 bp, avg 1,835 bp/read

| tool | mapped | mapped_frac | wall s | reads/s | agree w/ shmap | identity |
| --- | --- | --- | --- | --- | --- | --- |
| minshmap-py-nthash | 0 | 0.0 | 176.842 | 11 | 0.0 |  |
| minshmap-cpp-naive | 0 | 0.0 | 7.1 | 282 | 0.0 |  |
| minshmap-cpp-poly | 0 | 0.0 | 6.135 | 326 | 0.0 |  |
| minshmap-cpp-nthash | 0 | 0.0 | 5.151 | 388 | 0.0 |  |
| shmap | 2 | 0.001 | 1.751 | 1142 |  |  |

## Speed ranking (overall reads/s across all datasets)

| rank | tool | reads/s | total wall s |
| --- | --- | --- | --- |
| 1 | minshmap-cpp-nthash | 196.84 | 30.48 |
| 2 | minshmap-cpp-poly | 171.69 | 34.95 |
| 3 | shmap | 157.62 | 38.07 |
| 4 | minshmap-cpp-naive | 146.06 | 41.08 |
| 5 | minshmap-py-nthash | 2.54 | 2363.51 |
| 6 | minsh | 0.09 | 11.69 |

## Accuracy ranking (chr21 reads recovered + concordance)

Ranked by total reads mapped across all datasets (sensitivity); `mean_agree_with_shmap` is the placement concordance against the reference mapper, and `identity` is minSH's independent alignment score.

| rank | tool | total_mapped | mean_agree_with_shmap | identity |
| --- | --- | --- | --- | --- |
| 1 | shmap | 51 |  |  |
| 2 | minshmap-cpp-naive | 22 | 0.197 |  |
| 3 | minshmap-cpp-poly | 20 | 0.187 |  |
| 4 | minshmap-py-nthash | 17 | 0.157 |  |
| 5 | minshmap-cpp-nthash | 17 | 0.157 |  |
| 6 | minsh | 1 |  | 0.8691 |
