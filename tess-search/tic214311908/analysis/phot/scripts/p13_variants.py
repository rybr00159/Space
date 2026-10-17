"""Systematic uncertainty of the deblended depth from photometric-model choices.

Variants of p05 (all else identical): default; Mira template = SPOC PRF at the b-image position (instead of the
empirical b image); fit radius 2.0 and 4.0 px; other stars (T < 15.5 within the fit region) free per cadence;
M-dwarf templates broadened by a 0.15 px Gaussian; M-dwarf templates displaced by (+0.1, -0.1) and (-0.1, +0.1) px.
For each variant and sector: dip (e/s) at the global ephemeris on B:M3V, A:M3V, A:pair, A:Mira and the odd-even
depth difference of B:M3V (white errors; same detrending and quality cuts as p06).
Output: out/variants.json, out/variants.csv
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p05_phot import setup, photometry
from p06_lc import detrend, FREF, BKG_EXC_MAX
from tlib import ephem_frac, lin_depth

TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"]
VARIANTS = {"default": {}, "mira_prf": dict(mira_template="prf"), "rfit2": dict(rfit=2.0), "rfit4": dict(rfit=4.0),
            "free_others": dict(free_others=True), "blur015": dict(prf_perturb=dict(blur=0.15)),
            "shift+": dict(prf_perturb=dict(dx=0.1, dy=-0.1)), "shift-": dict(prf_perturb=dict(dx=-0.1, dy=0.1))}
KEYS = (("B", "M3V"), ("A", "M3V"), ("A", "pair"), ("A", "Mira"))


def measure_real(R, fref):
    t = R["t"]; chiB = R["FIT"]["B"]["chi2"]; bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
    ok = R["qok"] & (bexc < BKG_EXC_MAX) & R["jit_ok"] & np.isfinite(chiB) & (chiB < 3 * np.nanmedian(chiB[R["q0"]]))
    mask = np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1
    y = ephem_frac(np.nan_to_num(t), P, T0, T14, FL, R["exptime"])
    ep = np.round((t - T0) / P)
    out = {}
    for mk, nm in KEYS:
        names = R["FIT"][mk]["names"]; dF = R["FIT"][mk]["F"][:, names.index(nm)]
        raw = 1 + dF / fref[nm]; tr = detrend(t, raw, ok, mask); d = raw / tr
        m = ok & np.isfinite(d)
        sd = robust_std(d[m & ~mask]); e = np.full(m.sum(), sd)
        dep, de, _, _ = lin_depth(d[m], e, y[m])
        out[f"{mk}:{nm}"] = (dep * fref[nm], de * fref[nm])
        if mk == "B":
            oe = []
            for par in (0, 1):
                mm = m & ((ep % 2 == par) | (y == 0))
                oe.append(lin_depth(d[mm], np.full(mm.sum(), sd), y[mm])[:2])
            out["oddeven_B"] = ((oe[1][0] - oe[0][0]) * fref[nm], np.hypot(oe[1][1], oe[0][1]) * fref[nm])
    return out


rows = []
for s in SECTORS:
    S = setup(s)
    for vn, kw in VARIANTS.items():
        R = photometry(s, S=S, verbose=False, **kw)
        fref = dict(FREF, Mira=float(R["mu"]["F"] - 0.9 * (FREF["M3V"] + FREF["pair"])))
        r = measure_real(R, fref)
        row = dict(sector=s, variant=vn)
        for k, v in r.items(): row[k] = v[0]; row[k + "_err"] = v[1]
        rows.append(row)
        print(f"S{s} {vn:12s}: " + ", ".join(f"{k} {v[0]:.1f}+-{v[1]:.1f}" for k, v in r.items()), flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{ROOT}/out/variants.csv", index=False, float_format="%.3f")
summ = {}
for vn in VARIANTS:
    d = D[D.variant == vn]; s_ = {}
    for k in [f"{a}:{b}" for a, b in KEYS] + ["oddeven_B"]:
        w = 1 / d[k + "_err"] ** 2; mu = float(np.sum(w * d[k]) / w.sum()); s_[k] = (mu, float(1 / np.sqrt(w.sum())))
    summ[vn] = s_
    print(f"ALL {vn:12s}: " + ", ".join(f"{k} {v[0]:.1f}+-{v[1]:.1f}" for k, v in s_.items()))
save_json(summ, f"{ROOT}/out/variants.json")
