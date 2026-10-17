"""Referee r03: per-cadence auxiliary regressors, model-free: background B(t) = median of the faintest 20% of pixels
farther than 6 px from the Mira; pointing (dx, dy) from a linearised fit of each background-subtracted frame to the
sector median frame and its finite-difference gradients, using pixels > 4 px from the Mira with median flux above
the 70th percentile (field stars). Output: data/aux_s{sector}.npz (t, q, fin, B, dx, dy, scale)."""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
for s in SECTORS:
    C = Cut(s); g = C.stars(); mira = g[g.source_id == GID["Mira"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    dm = np.hypot(xx - mira.x, yy - mira.y)
    fl = C.f.reshape(C.n, -1)
    q0 = C.fin & (C.q == 0)
    med = np.nanmedian(fl[q0], 0)
    far = (dm > 6).ravel()
    faint = far & (med < np.percentile(med[far], 20))
    B = np.full(C.n, np.nan); B[C.fin] = np.median(fl[C.fin][:, faint], 1)
    M = np.median(C.f[q0] - B[q0, None, None], 0)
    gx = np.zeros_like(M); gy = np.zeros_like(M)
    gx[:, 1:-1] = (M[:, 2:] - M[:, :-2]) / 2; gy[1:-1, :] = (M[2:, :] - M[:-2, :]) / 2
    sel = ((dm > 4) & (M > np.percentile(M[dm > 4], 70))).ravel()
    sel &= ~np.isin(np.arange(C.ny * C.nx), np.r_[np.arange(C.nx), np.arange(C.nx) + (C.ny - 1) * C.nx])
    sel &= (xx.ravel() > 0) & (xx.ravel() < C.nx - 1)
    A = np.column_stack([M.ravel()[sel], -gx.ravel()[sel], -gy.ravel()[sel]])
    Y = (fl[C.fin][:, sel] - B[C.fin, None]).T
    c, *_ = np.linalg.lstsq(A, Y, rcond=None)
    sc = np.full(C.n, np.nan); dx = np.full(C.n, np.nan); dy = np.full(C.n, np.nan)
    sc[C.fin] = c[0]; dx[C.fin] = c[1] / c[0]; dy[C.fin] = c[2] / c[0]
    np.savez(f"{R}/data/aux_s{s}.npz", t=C.t, q=C.q, fin=C.fin, B=B, dx=dx, dy=dy, scale=sc)
    print(f"S{s}: B median {np.nanmedian(B[q0]):.0f}; dx rstd {rstd(dx[q0]):.4f} dy rstd {rstd(dy[q0]):.4f} px; scale rstd {rstd(sc[q0]):.4f}; n sel px {sel.sum()}", flush=True)
