"""Referee r01: per-sector WCS offset from isolated normal-colour field stars, and the PRF flux scale
(PRF-fit flux / TIC catalogue flux) for field stars, incl. red (BP-RP > 2) ones, used to put the M3V's catalogue flux
on the same scale as PRF-fitted transit flux losses.
Method: median image of QUALITY==0 cadences with low background; every star with T < 14 lying 2.5..17.5 px inside the
cutout, > 4 px from the Mira, with no Gaia neighbour brighter than T+2.5 within 2 px, is fitted (free x, y, flux,
plane background) within r <= 2.5 px after subtracting all other stars (T<17) at catalogue flux x global scale.
Output: out/wcs_s{sector}.json; print summary."""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from scipy.optimize import least_squares

tic = pd.read_csv(f"{R}/data/tic330.csv", dtype={"GAIA": str, "ID": str})
tic = tic[tic.GAIA.notna()].copy(); tic["GAIA"] = tic.GAIA.astype(np.int64)
TT = dict(zip(tic.GAIA, tic.Tmag))


def run(s):
    C = Cut(s); g = C.stars(17.5)
    g["Ttic"] = g.source_id.map(TT)
    good = C.fin & (C.q == 0)
    fl = C.f.reshape(C.n, -1)
    bk = np.nanpercentile(fl, 10, axis=1)
    good &= bk < np.nanmedian(bk[good]) + 3 * rstd(bk[good])
    img = np.median(C.f[good], 0)
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mira = g[g.source_id == GID["Mira"]].iloc[0]
    ins = g[(g.x > -4) & (g.x < C.nx + 3) & (g.y > -4) & (g.y < C.ny + 3)].reset_index(drop=True)
    # first pass: global flux scale using the catalogue fluxes (no positional refinement)
    rows = []
    cand = ins[(ins["T"] < 14.5) & (ins.x > 2.5) & (ins.x < C.nx - 3.5) & (ins.y > 2.5) & (ins.y < C.ny - 3.5)]
    cand = cand[np.hypot(cand.x - mira.x, cand.y - mira.y) > 4]
    for r in cand.itertuples():
        d = np.hypot(ins.x - r.x, ins.y - r.y)
        if np.any((d > 0) & (d < 1.5) & (ins["T"] < r.T + 1.5)): continue
        m = np.hypot(xx - r.x, yy - r.y) <= 2.5
        oth = np.zeros((C.ny, C.nx))
        for o in ins.itertuples():
            if o.source_id == r.source_id or np.hypot(o.x - r.x, o.y - r.y) > 5: continue
            oth += tflux(o.T) * C.psf(o.x, o.y)
        dat = (img - oth)[m]; w = 1 / np.sqrt(np.abs(img[m]) + 100)
        u = (xx[m] - r.x); v = (yy[m] - r.y)
        def fit(p):
            A = np.column_stack([C.psf(p[0], p[1])[m], np.ones(m.sum()), u, v])
            c, *_ = np.linalg.lstsq(A * w[:, None], dat * w, rcond=None); return A, c
        def res(p):
            A, c = fit(p); return (A @ c - dat) * w
        q = least_squares(res, [r.x, r.y], diff_step=1e-3, x_scale=0.05)
        A, c = fit(q.x)
        cov = np.linalg.pinv(q.jac.T @ q.jac) * max(1, 2 * q.cost / (m.sum() - 6))
        rows.append(dict(source_id=int(r.source_id), T=float(r.T), Ttic=float(r.Ttic) if np.isfinite(r.Ttic) else np.nan, bp_rp=float(r.bp_rp),
                         ruwe=float(r.ruwe), x=float(r.x), y=float(r.y), ox=float(q.x[0] - r.x), oy=float(q.x[1] - r.y),
                         ex=float(np.sqrt(cov[0, 0])), ey=float(np.sqrt(cov[1, 1])), flux=float(c[0]),
                         rms_rel=float(np.std(A @ c - dat) / c[0]), dist_mira=float(np.hypot(r.x - mira.x, r.y - mira.y))))
    d = pd.DataFrame(rows)
    d["ratio_tic"] = d.flux / tflux(d.Ttic)
    nrm = d[(d.bp_rp < 2.0) & (d.ruwe < 1.4) & (d.rms_rel < 0.05)]
    off = (float(np.median(nrm.ox)), float(np.median(nrm.oy)))
    offe = (float(1.2533 * rstd(nrm.ox) / np.sqrt(len(nrm))), float(1.2533 * rstd(nrm.oy) / np.sqrt(len(nrm))))
    red = d[(d.bp_rp >= 2.0) & (d.ruwe < 1.4)]
    out = dict(sector=s, n_good=int(good.sum()), offset_px=off, offset_err_px=offe, n_normal=len(nrm),
               scatter_px=(float(rstd(nrm.ox)), float(rstd(nrm.oy))),
               flux_ratio_tic_normal=(float(np.median(nrm.ratio_tic)), float(rstd(nrm.ratio_tic)), int(nrm.ratio_tic.notna().sum())),
               flux_ratio_tic_red=(float(np.median(red.ratio_tic)) if len(red) else np.nan, float(rstd(red.ratio_tic)) if len(red) > 2 else np.nan, int(len(red))),
               red_offsets=red[["source_id", "bp_rp", "ox", "oy", "ex", "ey", "ratio_tic", "T"]].values.tolist(),
               stars=rows)
    save_json(out, f"{R}/out/wcs_s{s}.json")
    print(f"S{s}: offset ({off[0]:+.3f}, {off[1]:+.3f}) +- ({offe[0]:.3f}, {offe[1]:.3f}) px from {len(nrm)} normal stars, scatter "
          f"({out['scatter_px'][0]:.3f}, {out['scatter_px'][1]:.3f}); flux/TIC normal {out['flux_ratio_tic_normal'][0]:.3f} rstd {out['flux_ratio_tic_normal'][1]:.3f} (n {out['flux_ratio_tic_normal'][2]}); "
          f"red n={len(red)} ratio {out['flux_ratio_tic_red'][0]:.3f}; red offsets {np.round(red[['bp_rp','ox','oy']].values,2).tolist()}", flush=True)


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
