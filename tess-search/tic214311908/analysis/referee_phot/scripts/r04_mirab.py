"""Referee r04: the Mira's empirical variable-light image b = d pixel / d F_Mira and its position.
Pixel time series (QUALITY==0, low background, 2-h bins) are regressed on [1, F_ap - <F_ap>, B, dx, dy] (variant 1) and
additionally on a linear time term (variant 2, more conservative: only the part of the Mira variation that is not a
linear trend defines b). F_ap = background-subtracted flux within 2.5 px of the Mira's Gaia position.
The b image is fitted within 3 px with the SPOC PRF at a free position (+ constant) and compared with the Mira's
Gaia position (field-star WCS from r01). Also: flux-weighted centroid of b vs that of the PRF at the Gaia position
(model-light check). Output: out/mirab_s{sector}.json, data/mirab_s{sector}.npz"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from scipy.optimize import least_squares


def binit(t, Y, w):
    k = np.floor((t - t[0]) / w).astype(int); u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    o = np.zeros((len(u),) + Y.shape[1:]); np.add.at(o, inv, Y); return o / n.reshape((-1,) + (1,) * (Y.ndim - 1)), n


for s in SECTORS:
    C = Cut(s); g = C.stars(); W = json.load(open(f"{R}/out/wcs_s{s}.json")); ox, oy = W["offset_px"]
    a = np.load(f"{R}/data/aux_s{s}.npz")
    mira = g[g.source_id == GID["Mira"]].iloc[0]; m3 = g[g.source_id == GID["M3V"]].iloc[0]
    mx, my = mira.x + ox, mira.y + oy
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    ok = C.fin & (C.q == 0) & (a["B"] < np.nanmedian(a["B"][C.q == 0]) + 60) & np.isfinite(a["dx"])
    fl = C.f.reshape(C.n, -1)
    ap = (np.hypot(xx - mx, yy - my) <= 2.5).ravel()
    Fap = fl[:, ap].sum(1) - a["B"] * ap.sum()
    res = {}
    for var in (1, 2):
        cols = [np.ones(ok.sum()), Fap[ok] - Fap[ok].mean(), a["B"][ok] - a["B"][ok].mean(), a["dx"][ok], a["dy"][ok]]
        if var == 2: cols.append(C.t[ok] - C.t[ok].mean())
        X = np.column_stack(cols)
        Xb, n = binit(C.t[ok], X, 2 / 24); Yb, _ = binit(C.t[ok], fl[ok], 2 / 24)
        kb = n >= 0.5 * np.median(n); Xb, Yb = Xb[kb], Yb[kb]
        c, *_ = np.linalg.lstsq(Xb, Yb, rcond=None)
        r = Yb - Xb @ c; cu = np.linalg.pinv(Xb.T @ Xb)
        b = c[1].reshape(C.ny, C.nx); be = np.sqrt((r ** 2).sum(0) / (len(Yb) - Xb.shape[1]) * cu[1, 1]).reshape(C.ny, C.nx)
        m = np.hypot(xx - mx, yy - my) <= 3.0
        w = 1 / be[m]
        def fit(p):
            A = np.column_stack([C.psf(p[0], p[1])[m], np.ones(m.sum())]); cc, *_ = np.linalg.lstsq(A * w[:, None], b[m] * w, rcond=None); return A, cc
        def rr(p):
            A, cc = fit(p); return (A @ cc - b[m]) * w
        best = min(((np.sum(rr([mx + dx_, my + dy_]) ** 2), [mx + dx_, my + dy_]) for dx_ in np.linspace(-1, 1, 9) for dy_ in np.linspace(-1, 1, 9)), key=lambda z: z[0])
        q = least_squares(rr, best[1], diff_step=1e-3, x_scale=0.05)
        cov = np.linalg.pinv(q.jac.T @ q.jac) * max(1, 2 * q.cost / (m.sum() - 4))
        off_en = C.en(mx, my, q.x[0] - mx, q.x[1] - my)
        err_en = np.abs(C.en(mx, my, np.sqrt(cov[0, 0]), 0))[0], np.abs(C.en(mx, my, 0, np.sqrt(cov[1, 1])))[1]
        # model-free check: flux-weighted centroid of b within 2 px vs centroid of the PRF at the Gaia position in the same pixels
        m2 = np.hypot(xx - mx, yy - my) <= 2.0
        Pg = C.psf(mx, my)
        cb = np.array([(b[m2] * xx[m2]).sum() / b[m2].sum(), (b[m2] * yy[m2]).sum() / b[m2].sum()])
        cp = np.array([(Pg[m2] * xx[m2]).sum() / Pg[m2].sum(), (Pg[m2] * yy[m2]).sum() / Pg[m2].sum()])
        cen_en = C.en(mx, my, *(cb - cp))
        m3_en = C.en(mx, my, m3.x + ox - mx, m3.y + oy - my)
        res[var] = dict(pos_px=q.x.tolist(), offset_EN_arcsec=off_en.tolist(), offset_err_arcsec=[float(e) for e in err_en],
                        chi2=float(2 * q.cost), npix=int(m.sum()), amp=float(fit(q.x)[1][0]),
                        centroid_offset_EN_arcsec=cen_en.tolist(), M3V_EN_from_Mira=m3_en.tolist(), frac_var=float(np.std(Fap[ok]) / np.mean(Fap[ok])))
        if var == 1: np.savez(f"{R}/data/mirab_s{s}.npz", b=b, be=be, mx=mx, my=my)
    save_json(res, f"{R}/out/mirab_s{s}.json")
    print(f"S{s}: Mira var {res[1]['frac_var']*100:.1f}%; b-image PRF offset from Gaia (E,N): v1 {np.round(res[1]['offset_EN_arcsec'],1).tolist()} +- {np.round(res[1]['offset_err_arcsec'],1).tolist()} chi2 {res[1]['chi2']:.0f}/{res[1]['npix']-4}; "
          f"v2 {np.round(res[2]['offset_EN_arcsec'],1).tolist()}; centroid check v1 {np.round(res[1]['centroid_offset_EN_arcsec'],1).tolist()} v2 {np.round(res[2]['centroid_offset_EN_arcsec'],1).tolist()}; M3V at {np.round(res[1]['M3V_EN_from_Mira'],1).tolist()}", flush=True)
