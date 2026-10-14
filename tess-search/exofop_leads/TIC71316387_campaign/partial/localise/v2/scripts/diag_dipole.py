"""Dipole content of stacked difference images: D = a*M + gx*dM/dx + gy*dM/dy + c at the direct-image target position.
gx, gy (e/s * px) measure a pointing-like shift; real eclipse stacks vs. the stack of all fake windows (null)."""
from ana import *
import json
RATIO = 0.49
res = {}
for sec in SECTORS:
    ctx = Ctx(sec)
    k = ctx.kind
    deep = np.where(k == 'deep')[0]; sh = np.where(k == 'shallow')[0]; fk = np.where(k == 'fake')[0]
    Ds, ss, Bs = ctx.stack(np.r_[deep, sh], np.r_[np.ones(len(deep)), RATIO * np.ones(len(sh))])
    rd = ctx.fit_direct(Bs, ss)
    x, y = rd['x'], rd['y']; h = 0.01
    M = ctx.psf(x, y); Mx = (ctx.psf(x + h, y) - ctx.psf(x - h, y)) / (2 * h); My = (ctx.psf(x, y + h) - ctx.psf(x, y - h)) / (2 * h)
    m = ctx.mask
    def dip(idx, an):
        D, s, B = ctx.stack(idx, an)
        w = 1 / s[m]
        A = np.column_stack([M[m], Mx[m], My[m], np.ones(m.sum())]) * w[:, None]
        c, *_ = np.linalg.lstsq(A, D[m] * w, rcond=None)
        cov = np.linalg.pinv(A.T @ A)
        return c, np.sqrt(np.diag(cov))
    out = {}
    for name, idx, an in [('deep', deep, np.ones(len(deep))), ('shallow', sh, np.ones(len(sh))),
                          ('fake_all_as_deep', fk, np.ones(len(fk)))]:
        c, e = dip(idx, an)
        # NB: for deep/shallow, units: D normalised by sbar -> peak amplitude units; gx/gy also divided by sbar
        # convert to raw e/s*px (multiply back by mean sbar)
        sb = ctx.z['sbar'][idx].mean()
        out[name] = dict(a=c[0], gx=c[1], gy=c[2], e_gx=e[1], e_gy=e[2], sbar=sb, gx_raw=c[1] * sb, gy_raw=c[2] * sb,
                         egx_raw=e[1] * sb, egy_raw=e[2] * sb, dx_equiv=-c[1] / c[0] if name != 'fake_all_as_deep' else np.nan,
                         dy_equiv=-c[2] / c[0] if name != 'fake_all_as_deep' else np.nan)
        print(sec, name, 'a=%.2f gx_raw=%.2f+-%.2f gy_raw=%.2f+-%.2f  (sbar %.2f)' % (c[0], c[1] * sb, e[1] * sb, c[2] * sb, e[2] * sb, sb))
    res[sec] = out
to_json(res, f'{R}/results/diag_dipole.json')
