# minSHmap cross-language benchmark (Python vs. C++)

- generated: 2026-06-20T18:47:33+00:00
- reference: 2 segment(s), 100,000 bp
- reads: 500 x ~305 bp
- params: {'k': 11, 'hfrac': 0.1, 'theta': 0.4, 'min_diff': 0.02}
- python: 3.14.0
- compiler: g++ (x86_64-posix-seh-rev0, Built by MinGW-Builds project) 15.2.0

| sketcher | lang | index s | ref bp/s | map s | reads/s | mapped | accuracy |
| --- | --- | --- | --- | --- | --- | --- | --- |
| naive | python | 0.4444 | 225012 | 0.5364 | 932 | 0.984 | 0.984 |
| naive | cpp | 0.009 | 11087703 | 0.0127 | 39487 | 0.978 | 0.978 |
| poly | python | 0.1168 | 856476 | 0.2386 | 2096 | 0.982 | 0.982 |
| poly | cpp | 0.0089 | 11297009 | 0.0151 | 33101 | 0.982 | 0.982 |
| nthash | python | 0.2261 | 442335 | 1.0182 | 491 | 0.98 | 0.98 |
| nthash | cpp | 0.0024 | 41385589 | 0.025 | 20013 | 0.98 | 0.98 |

C++ mapping speedup over Python (reads/s):
  - naive: 42.4x
  - poly: 15.8x
  - nthash: 40.8x
