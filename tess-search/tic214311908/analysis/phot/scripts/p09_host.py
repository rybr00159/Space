"""Which star dims? Forward-model test using the deblended light curves and the injection-calibrated response.

For each sector the dips measured at the transit ephemeris on six light curves (B:M3V, B:Mira, Bp:pair, A:M3V,
A:pair, A:Mira) form a vector m. A transit of total flux loss delta on host h predicts delta * k_h, with k_h the
response vector from the pixel injections (p07; per sector, mean over 7 fake ephemerides). The covariance C of m is
estimated empirically from 300 random fake ephemerides (same P, random phase, real transit windows excluded)
measured simultaneously on all six light curves (captures red noise and correlations between light curves),
plus the injection scatter of k_h (times the measured dip) added in quadrature to the diagonal.
chi2_h = min over delta of (m - delta k_h)^T C^-1 (m - delta k_h), per sector and summed (delta common to all sectors
and per-sector free variants). Hosts: M3V, M3V with a mismatched PSF (template error), pair, Mira.
Output: out/host_test.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from tlib import ephem_frac
from scipy.stats import chi2 as chi2dist

KEYS = ["B:M3V", "B:Mira", "Bp:pair", "A:M3V", "A:pair", "A:Mira"]
LCK = {"B:M3V": "M3V_B", "B:Mira": "Mira_B", "Bp:pair": "pair_Bp", "A:M3V": "M3V_A", "A:pair": "pair_A", "A:Mira": "Mira_A"}
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"]
INJ = json.load(open(f"{ROOT}/out/inject_summary.json"))
RNG = np.random.default_rng(7)


def dips(L, fref, P_, T0_, real_mask):
    t = L["t_tdb"]; ok = L["quality_ok"]; ex = float(L["exptime_s"]) / 86400
    y = ephem_frac(np.nan_to_num(t), P_, T0_, T14, FL, ex)
    out = []
    for k in KEYS:
        d = L[f"det_{LCK[k]}"]; m = ok & np.isfinite(d) & ~real_mask
        if np.sum(m & (y > 0.5)) < 2: return None
        A = np.column_stack([np.ones(m.sum()), -y[m]]); c, *_ = np.linalg.lstsq(A, d[m], rcond=None)
        out.append(c[1] * fref[k.split(":")[1]])
    return np.array(out)


res = {}
tot = {h: [] for h in ("M3V", "M3V_mismatch", "pair", "Mira")}
Mvec, Cmat, Kmat, Kstd = {}, {}, {}, {}
st = json.load(open(f"{ROOT}/out/lc_stats.json"))
for s in SECTORS:
    L = dict(np.load(f"{ROOT}/data/lc_s{s}.npz"))
    fref = st[str(s)]["fref"]
    t = L["t_tdb"]
    real = np.abs((((np.nan_to_num(t) - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1
    m = dips(L, fref, P, T0, np.zeros(len(t), bool))
    nulls = []
    for k in range(300):
        v = dips(L, fref, P, T0 + RNG.uniform(0.06, 0.94) * P, real)
        if v is not None: nulls.append(v)
    nulls = np.array(nulls)
    C = np.cov(nulls.T)
    Mvec[s], Cmat[s] = m, C
    row = dict(measured=dict(zip(KEYS, m.tolist())), null_std=dict(zip(KEYS, np.sqrt(np.diag(C)).tolist())),
               null_corr=(C / np.outer(np.sqrt(np.diag(C)), np.sqrt(np.diag(C)))).round(3).tolist(), hosts={})
    for h in tot:
        k_ = np.array([INJ["per_sector_mean"][f"{h}|{s}"][kk] for kk in KEYS])
        ks = np.array([INJ["per_sector_std"][f"{h}|{s}"][kk] for kk in KEYS])
        Kmat[(s, h)], Kstd[(s, h)] = k_, ks
        Ci = np.linalg.inv(C)
        dlt = float(k_ @ Ci @ m / (k_ @ Ci @ k_))
        C2 = C + np.diag((ks * max(dlt, 0)) ** 2)
        Ci2 = np.linalg.inv(C2)
        dlt = float(k_ @ Ci2 @ m / (k_ @ Ci2 @ k_)); r = m - dlt * k_
        c2 = float(r @ Ci2 @ r)
        row["hosts"][h] = dict(delta=dlt, delta_err=float(1 / np.sqrt(k_ @ Ci2 @ k_)), chi2=c2, dof=len(KEYS) - 1, p=float(chi2dist.sf(c2, len(KEYS) - 1)))
        tot[h].append(c2)
    res[s] = row
    print(f"S{s}: measured " + ", ".join(f"{k} {v:.1f}" for k, v in zip(KEYS, m)) + " | null sd " + ", ".join(f"{v:.1f}" for v in np.sqrt(np.diag(C))))
    print("      " + "; ".join(f"{h}: delta {v['delta']:.1f}+-{v['delta_err']:.1f} chi2 {v['chi2']:.1f}/5 (p {v['p']:.2g})" for h, v in row["hosts"].items()))
# combined with a common delta (same flux loss in every sector)
summ = {}
for h in tot:
    num = 0.0; den = 0.0
    for s in SECTORS:
        k_ = Kmat[(s, h)]; C = Cmat[s] + np.diag((Kstd[(s, h)] * 15.0) ** 2); Ci = np.linalg.inv(C)
        num += k_ @ Ci @ Mvec[s]; den += k_ @ Ci @ k_
    dlt = num / den; c2 = 0.0
    for s in SECTORS:
        k_ = Kmat[(s, h)]; C = Cmat[s] + np.diag((Kstd[(s, h)] * max(dlt, 0)) ** 2); r = Mvec[s] - dlt * k_
        c2 += float(r @ np.linalg.inv(C) @ r)
    dof = len(SECTORS) * len(KEYS) - 1
    summ[h] = dict(delta_common=float(dlt), delta_err=float(1 / np.sqrt(den)), chi2_common=c2, dof_common=dof, p_common=float(chi2dist.sf(c2, dof)),
                   chi2_persector_sum=float(sum(tot[h])), dof_persector=len(SECTORS) * (len(KEYS) - 1))
    print(f"ALL host {h:13s}: common delta {dlt:.1f}+-{1/np.sqrt(den):.1f} e/s, chi2 {c2:.1f}/{dof} (p {chi2dist.sf(c2, dof):.2g}); "
          f"sum of per-sector chi2 {sum(tot[h]):.1f}/{len(SECTORS)*(len(KEYS)-1)}")
best = min(summ, key=lambda h: summ[h]["chi2_common"])
for h in summ: summ[h]["dchi2_vs_best"] = summ[h]["chi2_common"] - summ[best]["chi2_common"]
save_json(dict(ephemeris=dict(P=P, T0=T0, T14=T14, flat=FL), per_sector={str(k): v for k, v in res.items()}, combined=summ), f"{ROOT}/out/host_test.json")
