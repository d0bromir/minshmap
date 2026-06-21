# minSHmap sketcher benchmark

- generated: 2026-06-20T09:42:09+00:00
- reference: 2 segment(s), 100,000 bp
- reads: 500 x ~305 bp
- params: {'k': 11, 'hfrac': 0.1, 'theta': 0.4, 'min_diff': 0.02}

| sketcher | index s | ref bp/s | map s | reads/s | mapped | accuracy |
| --- | --- | --- | --- | --- | --- | --- |
| naive | 0.2864 | 349209 | 0.5284 | 946 | 0.984 | 0.984 |
| poly | 0.1147 | 871770 | 0.2453 | 2038 | 0.982 | 0.982 |
| nthash | 0.2669 | 374607 | 0.7354 | 680 | 0.98 | 0.98 |
