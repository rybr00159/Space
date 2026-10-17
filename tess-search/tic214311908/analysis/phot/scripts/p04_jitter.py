"""Per-cadence pointing offsets (dx, dy) and field-star flux scale from the field outside 3 px of the Mira.

For each cadence: image = s * M0 + dx * dM0/dx + dy * dM0/dy + quadratic background, with M0 the static star model of
p02 (Gaia positions + global shift, p02 fluxes). Linear least squares per cadence, pixel weights from the temporal
scatter of each pixel. Also records the median background level in the faintest 15% of pixels (scattered-light flag).
Output: data/jitter_s{sector}.npz (t, dx, dy, scale, bkg, chi2r for every cadence; NaN where not finite)
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *


def static_model(C, W, prf, shift=(0.0, 0.0)):
    g = C.gaia(18.5)
    g = g[(g.x > -4) & (g.x < C.nx + 3) & (g.y > -4) & (g.y < C.ny + 3)].reset_index(drop=True)
    dx, dy = W["global_shift"]; ff = W["free_fluxes"]
    M = np.zeros((C.ny, C.nx))
    for s_, T, x, y in zip(g.source_id, g["T"], g.x, g.y):
        if s_ == IDS["Mira"]: continue
        f_ = ff.get(str(int(s_)), tmag_flux(T) * W["faint_scale"])
        M += f_ * prf.image(x + dx + shift[0], y + dy + shift[1])
    return M, g


def run(sector):
    C = Cutout(sector); W = json.load(open(f"{ROOT}/out/wcs_s{sector}.json"))
    prf = PRFModel(C); prf.blur = W["blur"]
    M0, g = static_model(C, W, prf)
    h = 0.05
    Gx = (static_model(C, W, prf, (h, 0))[0] - static_model(C, W, prf, (-h, 0))[0]) / (2 * h)
    Gy = (static_model(C, W, prf, (0, h))[0] - static_model(C, W, prf, (0, -h))[0]) / (2 * h)
    mira = g[g.source_id == IDS["Mira"]].iloc[0]; dx0, dy0 = W["global_shift"]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mask = np.hypot(xx - mira.x - dx0, yy - mira.y - dy0) > 3.0
    u = (xx[mask] - 10) / 10; v = (yy[mask] - 10) / 10
    A = np.column_stack([M0[mask], Gx[mask], Gy[mask], np.ones(mask.sum()), u, v, u * u, v * v, u * v])
    cube = C.flux.reshape(len(C.t), -1)[:, mask.ravel()]
    fin = C.fin
    # per-pixel noise from point-to-point scatter of good cadences
    good = fin & (C.q == 0)
    sig = 1.4826 * np.median(np.abs(np.diff(cube[good], axis=0)), axis=0) / np.sqrt(2) + 1e-3
    Aw = A / sig[:, None]
    out = np.full((len(C.t), A.shape[1]), np.nan); chi2r = np.full(len(C.t), np.nan)
    Y = (cube[fin] / sig).T
    coef, *_ = np.linalg.lstsq(Aw, Y, rcond=None)
    out[fin] = coef.T
    resid = Aw @ coef - Y
    chi2r[fin] = np.sum(resid ** 2, axis=0) / (mask.sum() - A.shape[1])
    s_ = out[:, 0]
    dx = out[:, 1] / s_; dy = out[:, 2] / s_
    allpix = C.flux.reshape(len(C.t), -1)
    med = np.nanmedian(allpix[good], 0)
    faint = med < np.percentile(med, 15)
    bkg = np.nanmedian(allpix[:, faint], 1)
    np.savez(f"{ROOT}/data/jitter_s{sector}.npz", t=C.t, q=C.q, dx=dx, dy=dy, scale=s_, bkg=bkg, chi2r=chi2r, fin=fin)
    gd = good & np.isfinite(dx)
    print(f"S{sector}: jitter rms (q=0) dx {np.nanstd(dx[gd]):.4f} dy {np.nanstd(dy[gd]):.4f} px; robust {robust_std(dx[gd]):.4f} {robust_std(dy[gd]):.4f}; "
          f"p1-p99 dx {np.nanpercentile(dx[gd],1):+.3f}..{np.nanpercentile(dx[gd],99):+.3f}; scale median {np.nanmedian(s_[gd]):.3f} rstd {robust_std(s_[gd]):.4f}; "
          f"chi2r median {np.nanmedian(chi2r[gd]):.1f}; bkg median {np.nanmedian(bkg[gd]):.0f} p99 {np.nanpercentile(bkg[gd],99):.0f}", flush=True)


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
