"""q07: look-elsewhere check of the post-hoc S39+S66 vs S103+S104 depth contrast (analyst c07 per-sector depths).
Null: six sector depths ~ N(mu, sigma_s) with the analyst's errors. Families: (a) 5 contiguous early/late
bipartitions; (b) all 31 bipartitions; (c) all pairs of disjoint groups with >= 2 sectors each.
Statistic: max |z| of the weighted group-mean difference over the family. p = P(max|z|_null >= z_obs).
Output: out/q07_lee.json"""
import json, itertools, numpy as np
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
S = [12, 39, 66, 93, 103, 104]
data = {"INV": ([19.7, 13.1, 11.9, 15.8, 13.7, 19.9], [6.3, 1.6, 2.5, 1.9, 2.1, 2.6]),
        "REF": ([12.2, 14.6, 14.5, 14.7, 17.5, 18.7], [3.8, 2.0, 2.7, 2.3, 1.9, 1.9]),
        "REF_noE58": ([12.2, 14.6, 14.5, 14.7, 16.0, 18.7], [3.8, 2.0, 2.7, 2.3, 2.0, 1.9])}
idx = range(6)
fam_a = [(tuple(range(0, k)), tuple(range(k, 6))) for k in range(1, 6)]
fam_b = []
for r in range(1, 6):
    for A in itertools.combinations(idx, r):
        B = tuple(i for i in idx if i not in A)
        if (B, A) not in fam_b: fam_b.append((A, B))
fam_c = []
for ra in range(2, 5):
    for A in itertools.combinations(idx, ra):
        rest = [i for i in idx if i not in A]
        for rb in range(2, len(rest) + 1):
            for B in itertools.combinations(rest, rb):
                if (B, A) not in fam_c: fam_c.append((A, B))
rng = np.random.default_rng(11)
out = {"n_family": dict(a=len(fam_a), b=len(fam_b), c=len(fam_c))}
def zstat(d, e, A, B):
    wA = 1 / e[..., list(A)] ** 2; wB = 1 / e[..., list(B)] ** 2
    mA = (d[..., list(A)] * wA).sum(-1) / wA.sum(-1); mB = (d[..., list(B)] * wB).sum(-1) / wB.sum(-1)
    return (mB - mA) / np.sqrt(1 / wA.sum(-1) + 1 / wB.sum(-1))
for name, (d, e) in data.items():
    d = np.array(d); e = np.array(e)
    zobs = float(zstat(d, e, (1, 2), (4, 5)))
    N = 200000
    sim = rng.normal(0, 1, (N, 6)) * e
    r = {"z_obs": zobs}
    for fn, fam in (("a", fam_a), ("b", fam_b), ("c", fam_c)):
        zmax = np.max(np.abs(np.column_stack([zstat(sim, np.broadcast_to(e, sim.shape), A, B) for A, B in fam])), axis=1)
        r[f"p_{fn}"] = float(np.mean(zmax >= abs(zobs)))
        r[f"maxz_obs_{fn}"] = float(np.max(np.abs([zstat(d, e, A, B) for A, B in fam])))
        r[f"p_maxobs_{fn}"] = float(np.mean(zmax >= r[f"maxz_obs_{fn}"]))
    w = 1 / e ** 2; mu = (d * w).sum() / w.sum(); r["chi2"] = float((((d - mu) / e) ** 2).sum())
    out[name] = r
    print(name, {k: round(v, 3) for k, v in r.items()})
json.dump(out, open(f"{R3}/out/q07_lee.json", "w"), indent=1)
