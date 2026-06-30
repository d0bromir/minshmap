# minSHmap cross-language benchmark (Python vs. C++)

> **Historical (obsolete).** This report predates the switch to the `minimizer-iter`
> library and still lists the old `naive`/`poly`/`nthash` sketchers. The current
> single-implementation synthetic benchmark (8,000 reads, `k=15 w=11 t=0.5`) gives
> **Python 14,499 reads/s** and **C++ 50,084 reads/s** (~3.5×), both mapping 3385
> reads at 100 % placement precision; on real chr21 long reads the two are
> algorithmically identical and emit byte-identical PAF. See
> [../realworld/results_rw/realworld.md](../realworld/results_rw/realworld.md).

- generated: 2026-06-21T19:19:58+00:00
- reference: 2 segment(s), 100,000 bp
- reads: 500 x ~305 bp
- params: {'k': 11, 'hfrac': 0.1, 'theta': 0.4, 'min_diff': 0.02}
- python: 3.14.0
- compiler: g++ (x86_64-posix-seh-rev0, Built by MinGW-Builds project) 15.2.0

| sketcher | lang | index s | ref bp/s | map s | reads/s | mapped | accuracy |
| --- | --- | --- | --- | --- | --- | --- | --- |
| naive | python | 0.3119 | 320658 | 0.4668 | 1071 | 0.984 | 0.984 |
| naive | cpp | 0.0065 | 15333895 | 0.0142 | 35285 | 0.978 | 0.978 |
| poly | python | 0.147 | 680489 | 0.3016 | 1658 | 0.982 | 0.982 |
| poly | cpp | 0.005 | 20179191 | 0.0165 | 30344 | 0.982 | 0.982 |
| nthash | python | 0.2044 | 489269 | 0.4013 | 1246 | 0.986 | 0.986 |
| nthash | cpp | 0.0033 | 30176836 | 0.0117 | 42894 | 0.986 | 0.986 |

C++ mapping speedup over Python (reads/s):
  - naive: 32.9x
  - poly: 18.3x
  - nthash: 34.4x
