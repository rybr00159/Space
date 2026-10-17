"""refclean: the referee's ('REF') local-window photometry (../referee_phot/scripts/r03_aux.py + r05_local.py) with
the field variables excluded from its field-wide estimators.

r03 builds, per cadence, B(t) = median of the faintest 20% of the pixels > 6 px from the Mira, and the pointing (dx, dy)
and scale from a linearised fit of the field-star pixels (> 4 px from the Mira, median flux above the 70th percentile).
CLEAN: pixels in cleanlib.var_mask(sector, mode) are removed from both pixel sets. Everything else (window baselines
per pixel: 1, tau, tau^2, B, dx, dy; PRF fit of the residual frames within 2.5 px of the M3V) is the referee's code,
imported read-only from ../referee_phot/scripts/r05_local.py.
"""
import os, sys
for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[k] = "1"
sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
import r05_local as R5
from rlib import Cut, GID, rstd, SECTORS
from r05_local import trap_shape, half_for, MASK
import numpy as np, json
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import var_mask, CL, RF, jd

RFR = R5.R          # ../referee_phot (read-only)


def aux_clean(s, mode, C=None, save=True):
    C = C or Cut(s); g = C.stars(); mira = g[g.source_id == GID["Mira"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    dm = np.hypot(xx - mira.x, yy - mira.y)
    vm = var_mask(s, mode, C.ny, C.nx).ravel()
    fl = C.f.reshape(C.n, -1)
    q0 = C.fin & (C.q == 0)
    med = np.nanmedian(fl[q0], 0)
    far = (dm > 6).ravel() & ~vm                                                  # CLEAN
    faint = far & (med < np.percentile(med[far], 20))
    B = np.full(C.n, np.nan); B[C.fin] = np.median(fl[C.fin][:, faint], 1)
    M = np.median(C.f[q0] - B[q0, None, None], 0)
    gx = np.zeros_like(M); gy = np.zeros_like(M)
    gx[:, 1:-1] = (M[:, 2:] - M[:, :-2]) / 2; gy[1:-1, :] = (M[2:, :] - M[:-2, :]) / 2
    sel = ((dm > 4) & (M > np.percentile(M[dm > 4], 70))).ravel()
    sel &= ~np.isin(np.arange(C.ny * C.nx), np.r_[np.arange(C.nx), np.arange(C.nx) + (C.ny - 1) * C.nx])
    sel &= (xx.ravel() > 0) & (xx.ravel() < C.nx - 1)
    sel &= ~vm                                                                    # CLEAN
    A = np.column_stack([M.ravel()[sel], -gx.ravel()[sel], -gy.ravel()[sel]])
    Y = (fl[C.fin][:, sel] - B[C.fin, None]).T
    c, *_ = np.linalg.lstsq(A, Y, rcond=None)
    sc = np.full(C.n, np.nan); dx = np.full(C.n, np.nan); dy = np.full(C.n, np.nan)
    sc[C.fin] = c[0]; dx[C.fin] = c[1] / c[0]; dy[C.fin] = c[2] / c[0]
    out = dict(t=C.t, q=C.q, fin=C.fin, B=B, dx=dx, dy=dy, scale=sc, nsel=int(sel.sum()), nfaint=int(faint.sum()))
    if save: np.savez(f"{CL}/data/refaux_{mode}_s{s}.npz", **out)
    return out


class SectorC(R5.Sector):
    """r05_local.Sector with the auxiliary regressors taken from aux_clean(mode) ('orig' = the referee's own file)."""
    def __init__(self, s, mode="all", aux=None):
        self.s = s; C = self.C = Cut(s)
        g = C.stars(); W = json.load(open(f"{RFR}/out/wcs_s{s}.json")); ox, oy = W["offset_px"]
        g["x"] += ox; g["y"] += oy; self.g = g
        if aux is not None: self.a = aux
        elif mode == "orig": self.a = dict(np.load(f"{RFR}/data/aux_s{s}.npz"))
        else:
            fn = f"{CL}/data/refaux_{mode}_s{s}.npz"
            self.a = dict(np.load(fn)) if os.path.exists(fn) else aux_clean(s, mode, C)
        z = np.load(f"{RFR}/data/tcorr_s{s}.npy"); self.ttdb = C.t + np.interp(C.t, z[0], z[1]) / 86400
        pos = {k: g[g.source_id == GID[k]].iloc[0] for k in ("Mira", "M3V", "pA", "pB", "bg3")}
        self.pos = {k: (float(v.x), float(v.y)) for k, v in pos.items()}
        yy, xx = np.mgrid[0:C.ny, 0:C.nx]; self.xx, self.yy = xx, yy
        mx, my = self.pos["Mira"]; m3x, m3y = self.pos["M3V"]
        self.reg = ((np.hypot(xx - mx, yy - my) <= 3.5) | (np.hypot(xx - m3x, yy - m3y) <= 3.5)).ravel()
        self.regM = (np.hypot(xx - m3x, yy - m3y) <= 2.5).ravel()
        self.P_M3V = C.psf(m3x, m3y)
        self.P_pair = (C.psf(*self.pos["pA"]) + 0.367 * C.psf(*self.pos["pB"])) / 1.367
        b = np.load(f"{RFR}/data/mirab_s{s}.npz")["b"]
        bm = np.where(np.hypot(xx - mx, yy - my) <= 4.0, b, 0.0); self.b_Mira = bm / bm.sum()
        self.P_Mira = C.psf(mx, my)
        q = C.q; B = self.a["B"]
        self.qok = C.fin & ((q == 0) | (q == 2048)) & np.isfinite(self.a["dx"])
        q0 = C.fin & (q == 0)
        self.Bmed = float(np.nanmedian(B[q0]))
        jit = np.ones(C.n, bool)
        for k in ("dx", "dy"):
            v = self.a[k]; jit &= np.abs(v - np.nanmedian(v[q0])) < 6 * rstd(v[q0]) + 0.02
        self.qok &= jit
        self.flat = C.f.reshape(C.n, -1)
        self.exptime = C.dt


def tiled(S, step=0.30, keep=0.15, bmax=60.0, **kw):
    """Contiguous REF-method light curve: windows every `step` d (half-width per sector as r05), keeping |t - tc| < keep
    of each (the per-pixel baseline masks |t - tc| < 0.07 d as for a transit). Returns t, aM (e-/s), sigma, window id."""
    t = S.ttdb[S.qok]
    cs = np.arange(t.min() + 0.2, t.max() - 0.2, step)
    T, A, SG, W = [], [], [], []
    for i, tc in enumerate(cs):
        r = S.window(tc, bmax=bmax, **kw)
        if r is None: continue
        m = np.abs(r["t"] - tc) < keep
        if m.sum() < 3: continue
        sg = rstd(r["aM"][r["oot"]])
        T.append(r["t"][m]); A.append(r["aM"][m]); SG.append(np.full(m.sum(), sg)); W.append(np.full(m.sum(), i))
    return np.concatenate(T), np.concatenate(A), np.concatenate(SG), np.concatenate(W)
