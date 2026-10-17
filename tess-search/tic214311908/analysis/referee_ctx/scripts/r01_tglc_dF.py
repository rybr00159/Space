"""R01: re-derive the TGLC Sector 39 (and 12) flux lost per star (e-/s) with an independent estimator.

Estimator: joint fit of a common trapezoid amplitude with a separate quadratic baseline in each +-0.30 d window around
each predicted transit (no global detrending). Noise: the same estimator applied to fake triplets of windows
(the 3 real windows shifted by a common offset dt, |dt| >= 0.6 d from any real or phase-0.5 window), which keeps the
cross-star correlations. Offsets are stepped by 0.05 d but only offsets >= 0.6 d apart are independent, so we also
report the number of independent samples.

Output: out/r01_tglc_dF.json
"""
import glob, json, sys, os
import numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(__file__))
import rlib as R

STARS = {"M3V": 5925373272850784768, "pairA": 5925373272850782208, "pairB": 5925373272827253632,
         "Mira": 5925373268528684288, "ref1": 5925373204131295872, "ref2": 5925373204131294336}
LC = "tic214311908/analysis/ctx/lc/tglc"


def load(sid, sector, col):
    f = glob.glob(f"{LC}/hlsp_tglc_tess_ffi_gaiaid-{sid}-s{sector:04d}-*.fits")[0]
    d = fits.getdata(f, 1)
    ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0) & np.isfinite(d[col])
    return np.array(d["time"][ok], float), np.array(d[col][ok], float)


def run(sector, col, half=0.30):
    ex = R.exptime(sector)
    res = {}
    t_all = {}
    for k, sid in STARS.items():
        t, f = load(sid, sector, col); t_all[k] = (t, f)
    t0 = t_all["M3V"][0]
    tcs, ks = R.centres(t0.min(), t0.max())
    tcs = tcs[(tcs > t0.min() - 0.1) & (tcs < t0.max() + 0.1)]
    t5, _ = R.centres(t0.min(), t0.max(), 0.5)
    avoid = np.concatenate([tcs, t5])
    # real
    for k, (t, f) in t_all.items():
        d, e, n, per = R.joint_fit(t, f, tcs, ex, half)
        res[k] = dict(dF=float(d), err_formal=float(e), n_win=n, per=[(float(a), float(b), float(c)) for a, b, c in per])
    # fake triplets with common offset
    dts = np.arange(-R.P / 2, R.P / 2, 0.05)
    good = []
    for dt in dts:
        sh = tcs + dt
        if np.min(np.abs(sh[:, None] - avoid[None, :])) < 2 * half:
            continue
        good.append(dt)
    good = np.array(good)
    F = np.full((len(good), len(STARS)), np.nan)
    for i, dt in enumerate(good):
        for j, (k, (t, f)) in enumerate(t_all.items()):
            d, e, n, _ = R.joint_fit(t, f, tcs + dt, ex, half)
            if n == len(tcs) or n >= 2:
                F[i, j] = d
    ok = np.all(np.isfinite(F), axis=1); F = F[ok]; good = good[ok]
    # independent subset: offsets at least 2*half apart
    ind = [0]
    for i in range(1, len(good)):
        if good[i] - good[ind[-1]] >= 2 * half:
            ind.append(i)
    C = np.cov(F.T); Ci = np.cov(F[ind].T)
    out = dict(sector=sector, column=col, half_window_d=half, n_transit_windows=len(tcs),
               dF={k: round(v["dF"], 2) for k, v in res.items()},
               err_formal={k: round(v["err_formal"], 2) for k, v in res.items()},
               err_fake_all={k: round(float(np.sqrt(C[j, j])), 2) for j, k in enumerate(STARS)},
               err_fake_indep={k: round(float(np.sqrt(Ci[j, j])), 2) for j, k in enumerate(STARS)},
               fake_mean={k: round(float(F[:, j].mean()), 2) for j, k in enumerate(STARS)},
               corr_fake_all=np.round(np.corrcoef(F.T), 3).tolist(),
               n_fake_all=int(len(good)), n_fake_indep=len(ind),
               per_window={k: v["per"] for k, v in res.items()})
    np.save(f"out/r01_fake_{col}_s{sector}.npy", F)
    np.save(f"out/r01_fake_dt_{col}_s{sector}.npy", good)
    return out


if __name__ == "__main__":
    allres = {}
    for sector in (39, 12):
        for col in ("psf_flux", "aperture_flux"):
            for half in ((0.30, 0.20, 0.45) if sector == 39 else (0.40,)):
                o = run(sector, col, half)
                allres[f"S{sector}_{col}_h{half}"] = o
                print(f"S{sector} {col} half={half}: dF", o["dF"])
                print("    err_formal", o["err_formal"])
                print("    err_fake  ", o["err_fake_all"], "n", o["n_fake_all"], "indep", o["n_fake_indep"])
                print("    fake mean ", o["fake_mean"])
    json.dump(allres, open("out/r01_tglc_dF.json", "w"), indent=1)
