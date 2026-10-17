"""Referee per-window transit-depth images (method 'LQ': local quadratic baseline fitted SIMULTANEOUSLY with the
transit template, per pixel, in a +-W window; NO time-domain background subtraction - background is absorbed later by
a plane in the image fit).  Linear in the pixel data, so an injected A*PRF transit with the template shape adds
exactly A*PRF to every window image.

For every sector: real windows (centred on T0 + nP) and fake windows (0.05-d grid, centres > W + T14/2 + 0.05 d from
any real transit so that no real transit enters a fake window).  Saved per window: D (depth image, e-/s, positive =
the pixel loses flux), V (formal variance), and the window centre.
Usage: rq03_img.py [W]  -> out/img_W{W}.npz  (default W = 0.30 d)
"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *

W = float(sys.argv[1]) if len(sys.argv) > 1 else 0.30
shape = dict(t14=T14, flat=FLAT)
if len(sys.argv) > 2:
    shape = dict(t14=float(sys.argv[2]) / 24, flat=float(sys.argv[3]))
tag = f'W{W:.2f}' + (f'_T{shape["t14"]*24:.2f}_f{shape["flat"]:.2f}' if len(sys.argv) > 2 else '')


def window(t, Y, tc, dt):
    m = np.abs(t - tc) <= W
    if m.sum() < 8: return None
    tt = t[m] - tc
    tau = trapezoid(t[m], tc, dt, **shape)
    exp_in = shape['t14'] * (1 + shape['flat']) / 2 / dt          # expected cadences-equivalent of transit
    if tau.sum() < 0.5 * exp_in: return None
    side = (W - shape['t14'] / 2) / dt
    if (tt < -shape['t14'] / 2).sum() < 0.4 * side or (tt > shape['t14'] / 2).sum() < 0.4 * side: return None
    X = np.column_stack([np.ones(m.sum()), tt, tt ** 2, -tau])
    XtXi = np.linalg.inv(X.T @ X)
    B = XtXi @ X.T @ Y[m]
    R = Y[m] - X @ B
    s2 = (R ** 2).sum(0) / (m.sum() - 4)
    return B[3], s2 * XtXi[3, 3], int(m.sum()), float(tau.sum())


out = {}
for sec in SECTORS:
    S = Sec(sec)
    gd = S.good
    t = S.t[gd]; Y = S.flux[gd].reshape(gd.sum(), -1)
    n = np.arange(np.ceil((t.min() - T0) / P) - 1, np.floor((t.max() - T0) / P) + 2)
    tr = T0 + n * P
    tr = tr[(tr > t.min() - 0.2) & (tr < t.max() + 0.2)]
    fk = np.arange(t.min() + W, t.max() - W, 0.05)
    fk = fk[np.all(np.abs(fk[:, None] - tr[None, :]) > W + shape['t14'] / 2 + 0.05, axis=1)]
    for kind, tcs in (('real', tr), ('fake', fk)):
        D, V, TC, NC = [], [], [], []
        for tc in tcs:
            r = window(t, Y, tc, S.dt)
            if r is None: continue
            D.append(r[0].reshape(S.ny, S.nx)); V.append(r[1].reshape(S.ny, S.nx)); TC.append(tc); NC.append(r[3])
        out[f's{sec}_{kind}_D'] = np.array(D, np.float32); out[f's{sec}_{kind}_V'] = np.array(V, np.float32)
        out[f's{sec}_{kind}_tc'] = np.array(TC); out[f's{sec}_{kind}_tau'] = np.array(NC)
    print(f'S{sec}: real windows {len(out[f"s{sec}_real_tc"])} of {len(tr)} ({np.round(out[f"s{sec}_real_tc"],3).tolist()}); fake {len(out[f"s{sec}_fake_tc"])}', flush=True)
np.savez_compressed(f'{ROOT}/out/img_{tag}.npz', **out)
