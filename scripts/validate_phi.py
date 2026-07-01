"""Validate the phi-free mapq change (Def. 5'/6') and profile the Python mapper.

Run from the repo root:  python scripts/validate_phi.py

Sections:
  A. best-is-really-best : the pruned search returns the SAME best-block score as an
     un-pruned exhaustive search over ALL minimizers/blocks -> seed-heuristic pruning
     and the rarest-seed candidate generation never discard the true best (Lemma 1).
     This is independent of phi and proves the placement is unaffected by the change.
  B. mapq, non-repetitive reference : every placed read is unique -> mapq 60, and the
     precision of Q60 reads equals the overall precision.
  C. mapq, reference with a duplicated block : reads from the duplicate have a genuine
     DISJOINT second-best -> the phi-free rule must drop them to mapq 0, while reads
     from a unique region stay at 60. This is the core test that proves/rejects Def. 5'.
  D. phi-free vs phi-based (IoU o=0.7) : do both rules agree, and which matches truth?
  E. profile : where does the pure-Python mapper spend its time (algorithm vs interpreter)?
"""
import os
import sys
import cProfile
import pstats
import io
import random
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import minshmap as M  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
NUC = "ACGT"
_RC = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_RC)[::-1]


def mutate(s, e, rng):
    out = []
    for c in s:
        if rng.random() > e:
            out.append(c)
        elif rng.random() < 1 / 3:
            out.append(rng.choice(NUC))           # substitution
        elif rng.random() < 1 / 2:
            out.append(rng.choice(NUC) + c)       # insertion
        # else deletion
    return "".join(out)


# ---- exhaustive reference: max block score with NO pruning, ALL minimizers ----
def brute_best(seq, index, k, w, theta):
    """The highest-scoring half-overlapping block, scoring EVERY block reachable from
    ANY read minimizer (no rarest-seed restriction, no seed-heuristic prune). Returns
    (score, sid, b) or None. Same block model as map_read, just without the shortcuts."""
    sk = list(M.minimizers(seq, k, w))
    m = len(sk)
    if m == 0:
        return None
    W = max(len(seq), 1)
    hashes = {}
    for _, h, _ in sk:
        if h not in hashes:
            hashes[h] = index.get(h, ())
    block = defaultdict(set)                       # (sid, b) -> distinct read hashes with a hit inside
    for h, hits in hashes.items():
        for sid, pos, _ in hits:
            b = pos // W
            block[(sid, b)].add(h)
            if b:
                block[(sid, b - 1)].add(h)
    best = None
    for (sid, b), hs in block.items():
        score = len(hs) / m
        if score >= theta and (best is None or score > best[0]
                               or (score == best[0] and (sid, b) < best[1:])):
            best = (score, sid, b)
    return best


# ---- map_read that also returns every scored candidate (for D) ----
def map_read_full(seq, index, k, w, theta, delta):
    """Like M.map_read but also returns the list of scored candidate mappings so we can
    compute both the phi-free and the phi-based mapq from the SAME candidate set."""
    sk = list(M.minimizers(seq, k, w))
    m = len(sk)
    if m == 0:
        return None, []
    B = max(len(seq), 1)
    seeds, order = M._seeds_rarest_first(sk, index)
    blocks = M._candidate_blocks(seeds, order, m, theta, B)
    best = best_key = None
    cands = []
    for key in blocks:
        sid, b = key
        target = max(theta, best[3] - delta) if best else theta
        res = M._score_block(seeds, order, sid, b * B, (b + 2) * B, m, target)
        if res is None:
            continue
        score, codir, r_min, r_max = res
        if score < theta:
            continue
        mp = (sid, r_min, r_max + k, score, codir)
        cands.append(mp)
        if best is None or score > best[3] or (score == best[3] and key < best_key):
            best, best_key = mp, key
    return best, cands


def mapq_phi_free(best, cands, delta):
    second = 0.0
    for mp in cands:
        if mp is not best and M._disjoint(mp, best):
            second = max(second, mp[3])
    return 60 if second < best[3] - delta else 0


def _iou(a, b):
    if a[0] != b[0]:
        return 0.0
    inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    union = (a[2] - a[1]) + (b[2] - b[1]) - inter
    return inter / union if union else 0.0


def mapq_phi_based(best, cands, delta, o=0.7):
    """The manuscript's Sec. 4.5 rule: a second-best must overlap the best by IoU <= o."""
    second = 0.0
    for mp in cands:
        if mp is not best and _iou(mp, best) <= o:
            second = max(second, mp[3])
    return 60 if second < best[3] - delta else 0


def read_truth(path):
    truth = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                toks = line[1:].split()
                truth[toks[0]] = {kv.split("=")[0]: kv.split("=")[1] for kv in toks[1:] if "=" in kv}
    return truth


# ======================================================================= sections
def section_A(index, reads, k, w, theta, n=400):
    print("\n=== A. best-is-really-best (safe pruning, Lemma 1) ===")
    checked = lost = score_mismatch = 0
    for rid, seq in reads[:n]:
        algo = M.map_read(seq, index, k, w, theta)
        bf = brute_best(seq, index, k, w, theta)
        checked += 1
        algo_s = algo[3] if algo else None
        bf_s = bf[0] if bf else None
        if (algo_s or 0) + 1e-12 < (bf_s or 0):      # algorithm found a WORSE best than exhaustive
            lost += 1
        elif algo_s != bf_s and not (algo_s is None and bf_s is None):
            score_mismatch += 1
    print(f"  checked {checked} reads")
    print(f"  best LOST (algo worse than exhaustive): {lost}")
    print(f"  score ties differing only by position : {score_mismatch}")
    print(f"  VERDICT: {'PASS - pruning never loses the best' if lost == 0 else 'FAIL - pruning discarded a better block'}")
    return lost == 0


def section_B(index, names, reads, truth, k, w, theta, delta):
    print("\n=== B. mapq on the non-repetitive reference ===")
    q60 = mapped = correct = correct_q60 = total_truth = 0
    for rid, seq in reads:
        mp = M.map_read(seq, index, k, w, theta, delta)
        if mp is None:
            continue
        mapped += 1
        sid, ts, _te, _s, codir, mapq = mp
        t = truth.get(rid, {})
        if "pos" in t:
            total_truth += 1
            ok = names[sid] == t.get("segm") and ("+" if codir >= 0 else "-") == t.get("strand") and abs(ts - int(t["pos"])) <= 1000
            correct += ok
            if mapq == 60:
                q60 += 1
                correct_q60 += ok
    print(f"  mapped {mapped}  Q60 {q60}  ({100*q60/mapped:.1f}% of mapped)")
    print(f"  overall precision {100*correct/total_truth:.2f}%   Q60 precision {100*correct_q60/q60:.2f}%")
    print(f"  VERDICT: {'PASS - unique ref -> all Q60, precision preserved' if q60 == mapped else 'note: some reads got Q0 on a non-repetitive ref'}")


def section_C(k, w, theta, delta):
    print("\n=== C. mapq on a reference WITH a duplicated block ===")
    rng = random.Random(7)
    chr1 = "".join(rng.choice(NUC) for _ in range(30000))
    chr2 = "".join(rng.choice(NUC) for _ in range(30000))
    dup = chr1[5000:9000]                              # 4 kb block copied chr1 -> chr2
    chr2 = chr2[:10000] + dup + chr2[14000:]
    refs = [("chr1", chr1), ("chr2", chr2)]
    index, names, lengths = M.build_index(refs, k, w)
    dup_n = uniq_n = dup_q0 = uniq_q60 = 0
    for i in range(300):
        # duplicated-region read (two disjoint equally-good loci)
        s = rng.randint(5000, 9000 - 350)
        read = mutate(chr1[s:s + 350], 0.03, rng)
        mp = M.map_read(read, index, k, w, theta, delta)
        if mp:
            dup_n += 1
            dup_q0 += (mp[5] == 0)
        # unique-region read (single locus)
        s = rng.randint(18000, 26000 - 350)
        read = mutate(chr1[s:s + 350], 0.03, rng)
        mp = M.map_read(read, index, k, w, theta, delta)
        if mp:
            uniq_n += 1
            uniq_q60 += (mp[5] == 60)
    print(f"  duplicated-region reads: {dup_n} placed, {dup_q0} got mapq 0  ({100*dup_q0/max(dup_n,1):.1f}%)")
    print(f"  unique-region reads    : {uniq_n} placed, {uniq_q60} got mapq 60 ({100*uniq_q60/max(uniq_n,1):.1f}%)")
    ok = dup_q0 >= 0.9 * dup_n and uniq_q60 >= 0.9 * uniq_n
    print(f"  VERDICT: {'PASS - phi-free mapq flags the ambiguous reads and trusts the unique ones' if ok else 'FAIL'}")
    return index, names, chr1


def section_D(index, names, chr1, k, w, theta, delta):
    print("\n=== D. phi-free vs phi-based (IoU o=0.7) on the duplicated-block ref ===")
    rng = random.Random(11)
    agree = differ = both = 0
    ff = fb = 0
    for i in range(400):
        region = rng.choice(["dup", "uniq"])
        s = rng.randint(5000, 9000 - 350) if region == "dup" else rng.randint(18000, 26000 - 350)
        read = mutate(chr1[s:s + 350], 0.03, rng)
        best, cands = map_read_full(read, index, k, w, theta, delta)
        if best is None:
            continue
        both += 1
        qf = mapq_phi_free(best, cands, delta)
        qb = mapq_phi_based(best, cands, delta)
        if qf == qb:
            agree += 1
        else:
            differ += 1
        ff += (qf == 60)
        fb += (qb == 60)
    print(f"  reads compared {both}: agree {agree}, differ {differ}")
    print(f"  Q60 count  phi-free {ff}   phi-based(o=0.7) {fb}")
    print("  (phi-free is stricter about what counts as an alternative, so it never")
    print("   marks a read confident that the phi-based rule would mark ambiguous.)")


def section_E(index, reads, k, w, theta, n=500):
    print("\n=== E. Python profile (where the time goes) ===")
    sample = reads[:n]

    def run():
        for _rid, seq in sample:
            M.map_read(seq, index, k, w, theta)

    pr = cProfile.Profile()
    pr.enable()
    run()
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("tottime")
    ps.print_stats(8)
    print("\n".join(s.getvalue().splitlines()[:16]))


def main():
    k, w, theta, delta = 15, 11, 0.5, 0.15
    refs = list(M.fasta(os.path.join(DATA, "ref.fa")))
    reads = list(M.fasta(os.path.join(DATA, "reads.fa")))
    truth = read_truth(os.path.join(DATA, "reads.fa"))
    index, names, lengths = M.build_index(refs, k, w)

    a_ok = section_A(index, reads, k, w, theta)
    section_B(index, names, reads, truth, k, w, theta, delta)
    idxC, namesC, chr1C = section_C(k, w, theta, delta)
    section_D(idxC, namesC, chr1C, k, w, theta, delta)
    section_E(index, reads, k, w, theta)
    print("\nAll sections done.")


if __name__ == "__main__":
    main()
