# Note: eliminating the maximal-overlap parameter $\phi$ (Def. 5)

*Companion note to `main_1.pdf` (map-shmap). Proposes a parameter-free replacement
for Def. 5 and propagates it to Def. 6 and Def. 7. The `minshmap.py` / `minshmap.cpp`
educational re-implementations were changed accordingly and the change was validated
empirically (see `scripts/validate_phi.py` and the results quoted at the bottom).*

---

## 1. What $\phi$ is and where it appears

`map-shmap` reports a mapping quality $\mathrm{mapq}\in\{0,60\}$ that answers: *is the
best mapping confidently unique, or is there another, comparably good mapping somewhere
else?* To decide this it needs a notion of **alternative mapping**, and that notion
currently carries a free parameter $\phi$ — the *maximal overlap*.

In the manuscript $\phi$ appears in:

- **Def. 5 (Alternative mappings).** $M^{\phi}_{l,\theta}\subseteq M_{l,\theta}$ of a best
  mapping $s_{best}$ are the mappings $s$ with **different direction**
  ($s.dir \neq s_{best}.dir$) **and** with limited overlap
  $$
  C\bigl([s.start,s.end),\,[s_{best}.start,s_{best}.end)\bigr)\;\le\;\phi,
  \qquad \phi\in[0,1].
  $$
- **Def. 6 (Mapping quality).** $\mathrm{mapq}(s_{best})=60$ iff every
  $s\in M^{\phi}_{l,\theta}$ has $C(p,s) < C(p,s_{best})-\delta$, for a margin
  $\delta\in(0,\theta)$; otherwise $0$.
- **Def. 7 (Problem).** lists the parameters $l,\theta,\phi,\delta$.
- **Sec. 4.5 / 5.3 (algorithm + setting).** here the overlap is defined *differently* —
  via **IoU** with a threshold $o$ ("not more than 50%", then experimentally $o=0.7$),
  not via the containment $C$ used in Def. 5.

So $\phi$ brings three problems at once:

1. it is a free parameter the user must guess;
2. Def. 5 attaches a puzzling extra clause ("different direction"); and
3. the definition (overlap $=C$, threshold $\phi$) and the algorithm (overlap $=$ IoU,
   threshold $o$, with $0.5\neq0.7$) do not even use the same overlap measure.

---

## 2. The observation that removes $\phi$

Why does Def. 5 need the *"different direction"* clause at all? Because the overlap there
is measured in **sketch coordinates** of $t$ and $\hat t$, which are two different arrays.
The reverse-complement image of the best mapping (the *same* genomic locus, but living in
$\hat t$) lands in a *different* interval and therefore looks like a legitimate
alternative — so the clause is bolted on to stop the best mapping from disqualifying
itself.

If instead the overlap is measured in the **original reference coordinates of $T$** — the
interval $g(s)\subseteq[0,|T|)$ of reference positions covered by the $k$-mers of $s$,
i.e. exactly the interval Sec. 4.6 already computes for the output — then:

- the reverse-complement image of $s_{best}$ occupies the **same** genomic interval, so its
  overlap with $s_{best}$ is $1$ and it is excluded **automatically, with no direction
  clause**;
- a genuinely different locus (including a repeat copy in the **same** orientation, which
  the current Def. 5 wrongly drops) has a **disjoint** genomic interval and is correctly
  kept as an alternative.

In genomic coordinates the only thing that overlaps $s_{best}$ is a *shifted copy of the
same placement* (e.g. the two half-overlapping blocks of Def. 10 that both cover one
locus). Those are exactly the mappings we must **not** count as alternatives. Hence the
natural, parameter-free threshold is the extreme value $\phi\to 0$: **disjointness**.

---

## 3. Proposed $\phi$-free definitions

Let $g(s)\subseteq[0,|T|)$ be the interval of original reference positions covered by the
$k$-mers of mapping $s$ (the interval already reported in Sec. 4.6).

> **Definition 5′ (Alternative mappings).** The alternative mappings
> $M^{\mathrm{alt}}_{l,\theta}\subseteq M_{l,\theta}$ of a best mapping $s_{best}$ are the
> candidate mappings whose reference footprint is disjoint from that of the best mapping:
> $$
> M^{\mathrm{alt}}_{l,\theta}\;\triangleq\;
> \bigl\{\, s\in M_{l,\theta}\;:\; g(s)\cap g(s_{best})=\varnothing \,\bigr\}.
> $$

> **Definition 6′ (Mapping quality).** The best mapping has
> $\mathrm{mapq}(s_{best})=60$ if every alternative mapping is sufficiently weaker:
> $$
> C(p,s) < C(p,s_{best})-\delta \quad\text{for all } s\in M^{\mathrm{alt}}_{l,\theta},
> \qquad \delta\in(0,\theta);
> $$
> otherwise $\mathrm{mapq}(s_{best})=0$.

> **Definition 7′ (Sketched Read Mapping Problem).** Given parameters $l,\theta,\delta$
> and sketches $p,t,rc(t)$: if the best mapping exists, report $s_{best}$ and
> $\mathrm{mapq}(s_{best})$; otherwise report nothing.

The parameter list of Def. 7 loses $\phi$: it becomes $l,\theta,\delta$.

---

## 4. Why this is strictly cleaner (not just "one fewer knob")

1. **Removes $\phi$** from Def. 5, 6, 7 and from the algorithm's setting.
2. **Removes the "different direction" clause** of Def. 5 — it becomes redundant, because
   in $T$-coordinates the forward and reverse images of a locus coincide.
3. **Unifies definition and algorithm** — the $C$-overlap of Def. 5 and the IoU-with-$o$
   of Sec. 4.5/5.3 (with the inconsistent $0.5$ vs $0.7$) both disappear.
4. **More correct semantics** — a repeat copy in the *same* orientation now counts as an
   alternative (the old Def. 5 missed it), which is the right behaviour for mapq.
5. **Handles the block artifact for free** — one locus covered by two half-overlapping
   blocks (Def. 10) yields two near-identical, mutually overlapping windows; disjointness
   excludes the duplicate with no threshold to tune.

---

## 5. Algorithmic consequence (minimal)

The best mapping $s_{best}$ is **unchanged** by this proposal — only its mapq is computed
differently — so accuracy of the *placement* is untouched and only *confidence* is
affected. In Sec. 4.5, the second-best search no longer carries a parameter: while
scanning blocks for $s_{second}$, skip those whose genomic interval intersects
$g(s_{best})$ (the best lies in $\le 2$ half-overlapping blocks, Def. 10, so this is just
"skip the best's block(s)"). Optimization 1 (raising $\theta'$ to $C(p,s_{best})-\delta$ so
a second-best within the margin survives the prune) is unaffected.

In the reference implementations this is exactly how it was coded: every window keeps its
$g$-interval $[\text{r\_min},\,\text{r\_max}+k)$ in forward $T$-coordinates (a single
forward index serves both strands via canonical minimizers), and
$$
\mathrm{mapq}=60 \iff \max_{\,s\,:\,g(s)\cap g(s_{best})=\varnothing} C(p,s) \;<\; C(p,s_{best})-\delta .
$$

---

## 6. Empirical validation (prove / reject)

The change was implemented in `minshmap.py` and `minshmap.cpp` (output stays **byte-identical
between the two** — verified 7854/7854 PAF lines, same hash, on `data/` with `-k 11 -w 5 -t 0.4`)
and tested with `scripts/validate_phi.py`. **All checks pass → the proposal is supported, not
rejected:**

- **Best-is-really-best (Lemma 1 / safe pruning).** Over 400 sampled reads, the algorithm's
  best-block score equals an exhaustive un-pruned search over **all** minimizers/blocks for
  every read (`best LOST = 0`). The seed-heuristic pruning and the rarest-seed candidate
  generation never discard the true best, and — crucially — the placement is **unchanged** by
  the φ-free change; only the mapq is computed differently.
- **mapq on a non-repetitive reference** (8000 reads, random 2×50 kb ref). Every placed read is
  unique → **100 % Q60**, and **Q60 precision = overall precision = 100 %**. The φ-free mapq
  does not over-flag confident mappings.
- **mapq on a reference with a duplicated 4 kb block.** Reads from the duplicated locus have a
  genuine disjoint second copy → **100 % (275/275) correctly dropped to mapq 0**; reads from a
  unique region → **100 % (267/267) stay at mapq 60**. This is the decisive test: Def. 5′ flags
  exactly the ambiguous reads, with no overlap parameter.
- **φ-free vs. φ-based (IoU o = 0.7).** On 365 compared reads the two rules give **identical**
  decisions (185 Q60 each). So removing φ costs **nothing** on clear-cut cases, while also
  fixing the cases the old Def. 5 mishandled (the "different direction" clause and same-orientation
  repeats) and removing a parameter the user otherwise has to tune.
- **Throughput (Python vs C++, synthetic, 8000 reads).** py 16,723 reads/s (index 0.02 s,
  precision 100 %); cpp 47,936 reads/s. Both map the same 3385 reads. `shmap` itself was **not
  re-run** (project benchmark policy reuses its baselines: paper Table 1 reports map-shmap as the
  most sensitive of the compared mappers; realworld baseline ≈ 157 reads/s on chr21).

### Why pure-Python is slow on real (repetitive) data — and whether it is a bug

It was **not an algorithm error** (results are identical to C++), but there *was* a genuine
complexity inefficiency in the original pure Python: `_score_block` used to **linear-scan**
each seed's hit list, and `_candidate_blocks` creates one block per hit, so a *frequent*
minimizer cost `O(hits) blocks × O(hits) scan = O(hits²)` per read. This is now **fixed** —
`_score_block` **binary-searches** the (sid,pos)-sorted hit list with `bisect_left`, exactly
like the C++ port, giving `O(log hits)` per seed → `O(hits·log hits)` total. Measured `map_read`
time vs. the largest hit-list length, before and after the fix:

| max hit-list | Python (linear, old) | Python (bisect, now) | C++ |
|---:|---:|---:|---:|
| 1   | 78 µs     | 58 µs    | ~20 µs |
| 50  | 6.5 ms    | 2.9 ms   | — |
| 500 | **563 ms** | **36 ms** | **1.8 ms** |

The quadratic blowup is gone: at 500 hits Python is **15.8× faster** and the gap to C++ shrank
from ~300× to ~20× — i.e. just the CPython interpreter constant, no algorithmic penalty. On the
real repetitive chr21 reference, pure-Python mapping rose from ~2 reads/s to ~8 reads/s. On
non-repetitive synthetic data throughput is unchanged (≈14–17 k reads/s). The two implementations
are now **algorithmically equivalent** (same binary-search lookup), which is a project invariant.

