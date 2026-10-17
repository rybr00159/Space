"""Raw-pixel aperture light curves at the M3V (own code): blind period search and transit shape.

Aperture: 3x3 box on the pixel nearest the M3V's apparent position (contains the Mira core) and r <= 2.5 px.
Background: per-cadence median of faint pixels (faintest 20% of pixels > 6 px from the Mira and > 3.5 px from any Gaia VARIABLE) x npix.
Systematics model per contiguous segment (gaps > 0.5 d split segments): cubic B-spline in time (knots every KN d)
+ field-star pointing dx, dy, dx^2, dy^2, dx*dy; 4-sigma iterative clipping.
  mode 'blind': transits NOT masked (for the period search);
  mode 'mask' : transits (+-1.2 T14) masked while fitting the systematics (for the shape fit).
Outputs: out/s03_lc_{mode}.npz (time, flux in ppt of 979 e-/s after EE correction, sector)
"""
import sys, numpy as np
sys.path.insert(0, 'tic214311908/analysis/panel_systematics/scripts')
from ps_lib import *
from scipy.interpolate import BSpline
from s01_depthimg import pointing

MODE = sys.argv[1] if len(sys.argv) > 1 else 'blind'
KN = float(sys.argv[2]) if len(sys.argv) > 2 else 0.75
g = gaia()


def bspl_basis(t, t0, t1, kn):
    nk = max(int(np.ceil((t1 - t0) / kn)), 1)
    inner = np.linspace(t0, t1, nk + 1)
    k = 3
    tk = np.r_[[t0] * k, inner, [t1] * k]
    nb = len(tk) - k - 1
    B = np.zeros((len(t), nb))
    for i in range(nb):
        c = np.zeros(nb); c[i] = 1
        B[:, i] = BSpline(tk, c, k, extrapolate=False)(t)
    return np.nan_to_num(B)


T, Fb, Fr, S = [], [], [], []
for sec in SECTORS:
    c = Cube(sec)
    pos = c.xy(g)
    sid = g.source_id.values
    im = int(np.where(sid == ID['M3V'])[0][0]); imi = int(np.where(sid == ID['Mira'])[0][0])
    x, y = pos[im]; xm, ym = pos[imi]
    dx, dy, _ = pointing(c, pos)
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    ix, iy = int(round(x)), int(round(y))
    box = np.zeros((c.ny, c.nx), bool); box[iy - 1:iy + 2, ix - 1:ix + 2] = True
    r25 = np.hypot(xx - x, yy - y) <= 2.5
    faint = np.hypot(xx - xm, yy - ym) > 6
    for i, (sx, sy) in enumerate(pos):
        if g.phot_variable_flag.values[i] == 'VARIABLE':
            faint &= np.hypot(xx - sx, yy - sy) > 3.5
    # faintest 20% of the remaining pixels (median image)
    medimg = np.median(c.flux[c.good], 0)
    faint &= medimg <= np.percentile(medimg[faint], 20)
    faint &= (xx > 0) & (yy > 0) & (xx < c.nx - 1) & (yy < c.ny - 1)
    m = c.good & np.isfinite(dx) & np.isfinite(dy)
    t = c.t[m]
    bk = np.median(c.flux[m][:, faint], 1)
    ee_b = c.prf(x, y)[box].sum(); ee_r = c.prf(x, y)[r25].sum()
    fb = c.flux[m][:, box].sum(1) - box.sum() * bk
    fr = c.flux[m][:, r25].sum(1) - r25.sum() * bk
    ddx, ddy = dx[m], dy[m]
    # segments
    br = np.where(np.diff(t) > 0.5)[0]
    segs = np.split(np.arange(len(t)), br + 1)
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    intr = np.abs(ph * P) < 1.2 * T14 / 2 if MODE == 'mask' else np.zeros(len(t), bool)
    outb, outr = np.full(len(t), np.nan), np.full(len(t), np.nan)
    for sg in segs:
        if len(sg) < 50:
            continue
        ts = t[sg]
        B = bspl_basis(ts, ts.min() - 1e-6, ts.max() + 1e-6, KN)
        pdx, pdy = ddx[sg] / 0.001, ddy[sg] / 0.001
        X = np.column_stack([B, pdx, pdy, pdx ** 2, pdy ** 2, pdx * pdy])
        for f, o in ((fb, outb), (fr, outr)):
            yv = f[sg]
            use = ~intr[sg]
            for it in range(5):
                bb, *_ = np.linalg.lstsq(X[use], yv[use], rcond=None)
                rr = yv - X @ bb
                s = 1.4826 * np.median(np.abs(rr[use] - np.median(rr[use])))
                new = (~intr[sg]) & (np.abs(rr) < 4 * s)
                if (new == use).all():
                    break
                use = new
            o[sg] = rr
    T.append(t); Fb.append(outb / ee_b / FM3V * 1e3); Fr.append(outr / ee_r / FM3V * 1e3); S.append(np.full(len(t), sec))
    print(sec, 'nbkg', int(faint.sum()), 'EE box %.3f r25 %.3f' % (ee_b, ee_r), 'npts', len(t), 'rms box ppt %.2f' % np.nanstd(outb / ee_b / FM3V * 1e3), flush=True)
T, Fb, Fr, S = map(np.concatenate, (T, Fb, Fr, S))
np.savez_compressed(f'{ROOT}/out/s03_lc_{MODE}_k{KN}.npz', t=T, box=Fb, r25=Fr, sec=S)
