"""Per-window per-pixel transit depth images from the TESScut cubes (real and fake epochs), own code.

For every window (+-0.30 d around a real or fake mid-time) and every pixel:
    F(t) = c0 + c1 dt + c2 dt^2 [+ kx dx(t) + ky dy(t)] - A * trapezoid(t)
A > 0 is a dimming (e-/s). dx, dy: field-star pointing (flux-weighted centroids of isolated T < 12.5 stars more than
5 px from the Mira), median-subtracted.  The same fit is applied to the SPOC background cube FLUX_BKG and to the
pointing series themselves (in-transit pointing excursion).
Outputs: out/s01_img_s{sec}.npz
"""
import sys, numpy as np
sys.path.insert(0, 'tic214311908/analysis/panel_systematics/scripts')
from ps_lib import *

HALF = 0.30
g = gaia()


def pointing(c, pos_all):
    xm, ym = c.xy(g[g.source_id == ID['Mira']])[0]
    sel = []
    for i, (x, y) in enumerate(pos_all):
        if not (2 <= x <= c.nx - 3 and 2 <= y <= c.ny - 3):
            continue
        if g['T'].values[i] > 12.5 or np.hypot(x - xm, y - ym) < 5:
            continue
        d = np.hypot(pos_all[:, 0] - x, pos_all[:, 1] - y)
        nb = (d < 2.0) & (d > 0) & (g['T'].values < g['T'].values[i] + 2.5)
        if nb.any():
            continue
        if g.phot_variable_flag.values[i] == 'VARIABLE':
            continue
        sel.append(i)
    F = c.flux
    DX, DY, used = [], [], []
    for i in sel:
        x, y = pos_all[i]
        ix, iy = int(round(x)), int(round(y))
        sub = F[:, iy - 1:iy + 2, ix - 1:ix + 2]
        s = np.nansum(sub, axis=(1, 2))
        xx, yy = np.meshgrid(np.arange(-1, 2), np.arange(-1, 2))
        cx = np.nansum(sub * xx, axis=(1, 2)) / s
        cy = np.nansum(sub * yy, axis=(1, 2)) / s
        if np.nanmedian(s) < 2000:
            continue
        DX.append(cx - np.nanmedian(cx[c.good])); DY.append(cy - np.nanmedian(cy[c.good])); used.append(int(g.source_id.values[i]))
    DX, DY = np.array(DX), np.array(DY)
    return np.nanmedian(DX, 0), np.nanmedian(DY, 0), used


def fitwin(c, tc, dx, dy, use_pt):
    m = c.good & (np.abs(c.t - tc) < HALF)
    t = c.t[m]
    tmpl = trap(t, tc, c.dt)
    nin_exp = T14 * 0.6 / c.dt
    nin = (tmpl > 0.5).sum()
    pre = ((t < tc - T14 / 2) & (t > tc - HALF)).sum(); post = ((t > tc + T14 / 2) & (t < tc + HALF)).sum()
    nside = (HALF - T14 / 2) / c.dt
    if nin < 0.5 * nin_exp or pre < 0.25 * nside or post < 0.25 * nside:
        return None
    d = t - tc
    cols = [np.ones_like(d), d, d ** 2]
    if use_pt:
        cols += [dx[m], dy[m]]
    cols.append(-tmpl)
    X = np.column_stack(cols)
    Y = c.flux[m].reshape(m.sum(), -1)
    B = c.bkg[m].reshape(m.sum(), -1)
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ X.T @ Y
    res = Y - X @ beta
    A = beta[-1].reshape(c.ny, c.nx)
    sig = np.sqrt(np.sum(res ** 2, 0) / (len(t) - X.shape[1]) * XtXi[-1, -1]).reshape(c.ny, c.nx)
    Ab = (XtXi @ X.T @ B)[-1].reshape(c.ny, c.nx)
    # pointing excursion (template amplitude in dx, dy with quadratic baseline only)
    X0 = np.column_stack([np.ones_like(d), d, d ** 2, -tmpl])
    bp = np.linalg.lstsq(X0, np.column_stack([dx[m], dy[m]]), rcond=None)[0][-1]
    # quality flags inside the full window (including flagged cadences)
    mw = np.isfinite(c.t) & (np.abs(c.t - tc) < HALF)
    mi = np.isfinite(c.t) & (np.abs(c.t - tc) < T14 / 2)
    return dict(A=A, sig=sig, Ab=Ab, ptx=bp[0], pty=bp[1], n=len(t), nin=int(nin),
                qbits_win=int(np.bitwise_or.reduce(c.q[mw])) if mw.any() else 0,
                qbits_in=int(np.bitwise_or.reduce(c.q[mi])) if mi.any() else 0,
                nflag_in=int((mi & ~c.good).sum()), nin_all=int(mi.sum()),
                bkg_in=float(np.nanmedian(c.bkg[m][:, 9, 9])), mean_img=np.nanmean(c.flux[m], 0))


def main(sec):
    c = Cube(sec)
    pos = c.xy(g)
    dx, dy, used = pointing(c, pos)
    tg = c.t[c.good]
    e0, e1 = int(np.floor((tg.min() - T0) / P)) - 1, int(np.ceil((tg.max() - T0) / P)) + 1
    md = c.momentum_dumps()
    out = {}
    for use_pt in (1, 0):
        R, Fk = [], []
        for E in range(e0, e1 + 1):
            tc = T0 + E * P
            r = fitwin(c, tc, dx, dy, use_pt)
            if r is not None:
                r.update(E=E, tc=tc, dump_dt=float(np.min(np.abs(md - tc))) if len(md) else np.nan)
                R.append(r)
            for ph in np.arange(0.12, 0.885, 0.055):
                tf = T0 + (E + ph) * P
                if np.min(np.abs(((tf - T0) / P + 0.5) % 1 - 0.5)) * P < 0.62:
                    continue
                rf = fitwin(c, tf, dx, dy, use_pt)
                if rf is not None:
                    rf.update(E=E, tc=tf, ph=ph, dump_dt=float(np.min(np.abs(md - tf))) if len(md) else np.nan)
                    Fk.append(rf)
        tag = 'pt' if use_pt else 'nopt'
        for nm, L in (('real', R), ('fake', Fk)):
            for k in ('A', 'sig', 'Ab', 'mean_img'):
                out[f'{tag}_{nm}_{k}'] = np.array([r[k] for r in L])
            for k in ('ptx', 'pty', 'n', 'nin', 'qbits_win', 'qbits_in', 'nflag_in', 'nin_all', 'bkg_in', 'E', 'tc', 'dump_dt'):
                out[f'{tag}_{nm}_{k}'] = np.array([r[k] for r in L])
    out['mean'] = np.nanmedian(c.flux[c.good], 0)
    out['pos'] = pos
    out['source_id'] = g.source_id.values
    out['dt'] = c.dt
    out['used_pointing_stars'] = np.array(used)
    out['dumps'] = md
    out['tgood'] = np.array([tg.min(), tg.max()])
    np.savez_compressed(f'{ROOT}/out/s01_img_s{sec}.npz', **out)
    print(sec, 'real', len(out['pt_real_A']), 'fake', len(out['pt_fake_A']), 'pointing stars', len(used), flush=True)


if __name__ == '__main__':
    for s in (map(int, sys.argv[1:]) if len(sys.argv) > 1 else SECTORS):
        main(s)
