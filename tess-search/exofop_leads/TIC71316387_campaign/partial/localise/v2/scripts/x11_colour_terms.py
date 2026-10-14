"""Fit the colour-dependent (lateral colour) astrometric term from x10 per-star tables.

Model per component (pixel frame):
  r = poly2(x, y) + m*(G-12) + k_rad*(C-C0)*u_rad + k_tan*(C-C0)*u_tan  [+ optional quadratic colour term]
Robust: iterative 4-sigma clipping with an error floor. Errors from a bootstrap over stars.
Also returns binned radial residuals vs BP-RP (after removing poly + magnitude terms) for plotting, and the
predicted radial shift (arcsec, + = towards boresight) of a source of colour C relative to the target (C0=0.786).
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from loclib import R, to_json

C0 = 0.786
FLOOR = float(os.environ.get('X11_FLOOR', 0.01))   # px error floor per star
QUAD = os.environ.get('X11_QUAD', '0') == '1'


def design(a, xb, yb, xt, yt):
    x = (a['xw'] - xt) / 300.0; y = (a['yw'] - yt) / 300.0
    ux = xb - a['xw']; uy = yb - a['yw']; n = np.hypot(ux, uy); ux /= n; uy /= n   # towards boresight
    tx, ty = uy, -ux
    dc = a['C'] - C0
    poly = [np.ones_like(x), x, y, x * x, x * y, y * y]
    dg = a['G'] - 12.0
    # stacked design for both components: unknowns = poly_x(6), poly_y(6), m_x, m_y, k_rad, k_tan [, q_rad]
    nP = len(poly)
    N = len(a)
    cols_x = []; cols_y = []
    for p in poly:
        cols_x.append(np.r_[p, np.zeros(N)])
    for p in poly:
        cols_y.append(np.r_[np.zeros(N), p])
    A = np.column_stack(cols_x + cols_y + [np.r_[dg, np.zeros(N)], np.r_[np.zeros(N), dg],
                                          np.r_[dc * ux, dc * uy], np.r_[dc * tx, dc * ty]] +
                        ([np.r_[dc ** 2 * ux, dc ** 2 * uy]] if QUAD else []))
    b = np.r_[a['xf'] - a['xw'], a['yf'] - a['yw']]
    e = np.r_[np.hypot(a['ex'], FLOOR), np.hypot(a['ey'], FLOOR)]
    return A, b, e, (ux, uy, tx, ty)


def robust_fit(A, b, e, nclip=6):
    ok = np.ones(len(b), bool)
    for _ in range(nclip):
        w = 1 / e[ok]
        p, *_ = np.linalg.lstsq(A[ok] * w[:, None], b[ok] * w, rcond=None)
        r = (b - A @ p) / e
        s = 1.4826 * np.median(np.abs(r[ok]))
        new = np.abs(r) < 4 * max(s, 1e-9)
        n = len(b) // 2
        new = new[:n] & new[n:]; new = np.r_[new, new]
        if np.array_equal(new, ok):
            break
        ok = new
    return p, ok, s


def analyse(files, label):
    As, bs, es, meta = [], [], [], []
    per = []
    for f in files:
        z = np.load(f)
        a = z['stars']
        good = (a['chi2r'] < 20) & (a['ruwe'] < 1.4) & (a['peak'] < 1.5e5)
        a = a[good]
        A, b, e, u = design(a, float(z['xb']), float(z['yb']), float(z['xt']), float(z['yt']))
        p, ok, s = robust_fit(A, b, e)
        # bootstrap over stars
        N = len(a); rng = np.random.default_rng(1)
        boots = []
        for _ in range(300):
            i = rng.integers(0, N, N); ii = np.r_[i, i + N]
            pb, okb, _ = robust_fit(A[ii], b[ii], e[ii], nclip=3)
            boots.append(pb[14:])
        boots = np.array(boots)
        k = p[14:]; ek = boots.std(0)
        # residual after removing poly + mag terms, projected radially
        rr = b - A[:, :14] @ p[:14]
        n = N
        rad = rr[:n] * u[0] + rr[n:] * u[1]; tan = rr[:n] * u[2] + rr[n:] * u[3]
        okk = ok[:n]
        bins = [(-0.2, 0.6), (0.6, 0.75), (0.75, 0.9), (0.9, 1.05), (1.05, 1.25), (1.25, 1.5), (1.5, 1.9), (1.9, 2.4), (2.4, 4.0)]
        bt = []
        for lo, hi in bins:
            s_ = okk & (a['C'] >= lo) & (a['C'] < hi)
            if s_.sum() >= 3:
                bt.append(dict(lo=lo, hi=hi, n=int(s_.sum()), C_med=float(np.median(a['C'][s_])),
                               rad=float(np.median(rad[s_])), rad_err=float(1.2533 * np.std(rad[s_]) / np.sqrt(s_.sum())),
                               tan=float(np.median(tan[s_])), tan_err=float(1.2533 * np.std(tan[s_]) / np.sqrt(s_.sum()))))
        nM = int(((a['C'] > 1.9) & okk).sum())
        res = dict(file=os.path.basename(f), sec=int(z['sec']), cam=int(z['cam']), ccd=int(z['ccd']), n=int(N), n_used=int(okk.sum()),
                   n_red_gt1p9=nM, rms_clip=float(s), k_rad=float(k[0]), k_rad_err=float(ek[0]), k_tan=float(k[1]),
                   k_tan_err=float(ek[1]), m_x=float(p[12]), m_y=float(p[13]),
                   q_rad=(float(k[2]) if QUAD else None), q_rad_err=(float(ek[2]) if QUAD else None), bins=bt)
        per.append(res)
        print(f"{res['file'][:40]} S{res['sec']} n={res['n_used']}/{N} (red>1.9: {nM}) k_rad={k[0]*1e3:+.1f}+-{ek[0]*1e3:.1f} mpx/mag "
              f"k_tan={k[1]*1e3:+.1f}+-{ek[1]*1e3:.1f}" + (f" q_rad={k[2]*1e3:+.1f}+-{ek[2]*1e3:.1f}" if QUAD else '') +
              f"  mag-term ({p[12]*1e3:+.1f},{p[13]*1e3:+.1f}) mpx/mag", flush=True)
        for bb in bt:
            print(f"     C {bb['lo']:.2f}-{bb['hi']:.2f} (n={bb['n']:3d}, med {bb['C_med']:.2f}): rad {bb['rad']*1e3:+6.1f}+-{bb['rad_err']*1e3:4.1f} mpx  tan {bb['tan']*1e3:+6.1f}+-{bb['tan_err']*1e3:4.1f}")
    to_json(per, f'{R}/results/x11_colour_terms{label}.json')
    return per


if __name__ == '__main__':
    files = sys.argv[1:] or sorted(glob.glob(f'{R}/results/x10_*.npz'))
    analyse(files, os.environ.get('X11_LABEL', ''))
