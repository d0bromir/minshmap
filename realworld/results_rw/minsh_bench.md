# minSH real-world alignment benchmark: original vs patched (banded) A*

- reference: T2T-CHM13v2.0 chr21 (real sequence content)
- tasks: 54 (lengths [1000, 2000, 4000], 6/length/category)
- per-alignment timeout: 8.0s; patched budget: 30% edits
- 'original' = align(max_cost=None) = upstream minSH; 'patched' = align(max_cost=30%·L)

`verdict` = finished without timing out (OK aligned or REJECT pruned). `mean_s` is over tasks that returned (timeouts counted at the cap).

| category | n | orig verdicts | orig timeouts | orig total s | patched verdicts | patched timeouts | patched total s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hifi | 18 | 17 | 1 | 32.1 | 17 | 1 | 30.3 |
| ont | 18 | 12 | 6 | 60.1 | 13 | 5 | 58.3 |
| decoy | 18 | 5 | 13 | 125.4 | 12 | 6 | 73.0 |
| ALL | 54 | 34 | 20 | 217.6 | 42 | 12 | 161.5 |

## Verdict

- exactness: of 29 tasks both versions aligned, **29/29** have identical edit distance (no mismatches).
- timeouts: original **20**, patched **12**.
- total wall: original **217.6s**, patched **161.5s** (**1.35x** faster overall).
- improvement: YES — the budget keeps results exact while bounding the worst case.
