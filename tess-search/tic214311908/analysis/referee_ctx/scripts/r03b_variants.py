"""R03b: robustness of the direct S39 pixel fit (r03) to WCS calibration, fit radius, background handling and
baseline window. Reports chi2 differences relative to the M3V and the free-fit offset (arcsec E,N from the M3V).
Output: out/r03b_variants.json
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
import leak as L, rlib as R

def run(use_cal=True, radius=3.5, bkg=True, half=0.30, deg=2):
    F = L.Field(39, use_cal=use_cal)
    with fits.open(F.file) as h:
        d = h[1].data; t = np.array(d["TIME"], float); fl = np.array(d["FLUX"], float); q = d["QUALITY"]
    ok = np.isfinite(t) & ((q == 0) | (q == 2048)) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    t, fl = t[ok], fl[ok]; yy, xx = np.mgrid[0:21, 0:21]
    if bkg:
        med = np.median(fl, axis=0); faint = med <= np.percentile(med, 30)
        Ab = np.column_stack([np.ones(faint.sum()), xx[faint], yy[faint]])
        c = np.linalg.lstsq(Ab, fl[:, faint].T, rcond=None)[0]
        fl = fl - (c[0][:, None, None] + c[1][:, None, None] * xx + c[2][:, None, None] * yy)
    ex = R.exptime(39)
    tcs, _ = R.centres(t.min(), t.max()); tcs = tcs[(tcs > t.min() + 0.3) & (tcs < t.max() - 0.3)]
    t5, _ = R.centres(t.min(), t.max(), 0.5); avoid = np.concatenate([tcs, t5])
    def amp(centres):
        rows, ys = [], []; nb = len(centres); npar = 1 + (deg + 1) * nb
        for b, tc in enumerate(centres):
            m = np.abs(t - tc) < half; s = R.shape(t[m], tc, ex); x = (t[m] - tc) / half
            A = np.zeros((m.sum(), npar)); A[:, 0] = -s
            for k in range(deg + 1): A[:, 1 + (deg + 1) * b + k] = x ** k
            rows.append(A); ys.append(fl[m].reshape(m.sum(), -1))
        return np.linalg.lstsq(np.vstack(rows), np.vstack(ys), rcond=None)[0][0].reshape(21, 21)
    Areal = amp(tcs)
    cs = np.arange(t.min() + half, t.max() - half, 0.02); cs = cs[np.min(np.abs(cs[:, None] - avoid[None, :]), axis=1) > 2 * half]
    sel = [cs[0]]
    for c in cs[1:]:
        if c - sel[-1] >= 2 * half: sel.append(c)
    fakes = np.array([amp([c]) for c in sel if (np.abs(t - c) < half).sum() > 40])
    var = np.var(fakes, axis=0, ddof=1) / len(tcs)
    x3, y3 = F.xy(L.STARS["M3V"]); fm = np.hypot(xx - x3, yy - y3) < radius; w = 1 / var[fm]
    def fit(src):
        P = F.stamp(*src)[fm]; y = Areal[fm]; A = np.column_stack([P, np.ones(len(P))]); Aw = A * w[:, None]
        c = np.linalg.solve(A.T @ Aw, Aw.T @ y); r = y - A @ c; return float(c[0]), float(np.sum(r * r * w))
    def sh(xy, a, pa=165.0):
        u = np.array([np.sin(np.radians(pa)), np.cos(np.radians(pa))]); dd = F.Ji @ (u * a); return (xy[0] + dd[0], xy[1] + dd[1])
    H = {"M3V": (x3, y3), "M3V_colour3.5": sh((x3, y3), 3.5), "G17.5_3.1": F.xy(L.EXTRA["G17.5_3.1"]),
         "G17.7_5.5": F.xy(L.EXTRA["G17.7_5.5"]), "pairA": F.xy(L.STARS["pairA"]), "Mira_gaia": F.xy(L.STARS["Mira"])}
    res = {h: fit(v) for h, v in H.items()}
    m3en = F.xy_to_en(x3, y3); best = (1e18, None)
    for dn in np.arange(-15, 15.1, 0.5):
        for de in np.arange(-35, 5.1, 0.5):
            c2 = fit(F.en_to_xy(de, dn))[1]
            if c2 < best[0]: best = (c2, (de, dn))
    out = dict(n_pix=int(fm.sum()), n_fake=len(fakes), amp_M3V=res["M3V"][0], chi2_M3V=res["M3V"][1],
               dchi2_vs_M3V={h: round(v[1] - res["M3V"][1], 2) for h, v in res.items()},
               free_chi2=best[0], free_offset_from_M3V_EN=[round(best[1][0] - m3en[0], 1), round(best[1][1] - m3en[1], 1)])
    return out

V = {}
for tag, kw in {"baseline": {}, "no_wcs_cal": dict(use_cal=False), "radius2.5": dict(radius=2.5), "radius5": dict(radius=5.0),
                "no_bkg_plane": dict(bkg=False), "half0.20": dict(half=0.20), "half0.45": dict(half=0.45), "linear_baseline": dict(deg=1)}.items():
    V[tag] = run(**kw); print(tag, V[tag], flush=True)
json.dump(V, open("out/r03b_variants.json", "w"), indent=1, default=float)
