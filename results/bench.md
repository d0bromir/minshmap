# minSHmap sketcher benchmark

> **Historical (obsolete).** This report compares the old `naive`/`poly`/`nthash`
> FracMinHash sketchers. minSHmap has since been unified on a single canonical
> `(w, k)`-minimizer from the `minimizer-iter` library, so there is no longer a
> sketcher choice. For the current Python-vs-C++ numbers see
> [../realworld/results_rw/realworld.md](../realworld/results_rw/realworld.md)
> (latest synthetic run: py 14,499 reads/s, cpp 50,084 reads/s, both mapped 3385).

- generated: 2026-06-21T19:19:33+00:00
- reference: 2 segment(s), 100,000 bp
- reads: 500 x ~305 bp
- params: {'k': 11, 'hfrac': 0.1, 'theta': 0.4, 'min_diff': 0.02}

| sketcher | index s | ref bp/s | map s | reads/s | mapped | accuracy |
| --- | --- | --- | --- | --- | --- | --- |
| naive | 0.3344 | 299032 | 0.6174 | 810 | 0.984 | 0.984 |
| poly | 0.1663 | 601174 | 0.3938 | 1270 | 0.982 | 0.982 |
| nthash | 0.2012 | 496953 | 0.3828 | 1306 | 0.986 | 0.986 |
