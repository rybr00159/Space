"""Red end of the TESS lateral-colour relation.

Per FFI (wide sample, x10w_*):  fit poly2(x,y) + magnitude term + linear colour term (radial & tangential, scaled by
field angle theta/theta_target) using only stars with BP-RP < 1.7 (robust).  Then for ALL stars compute the radial
residual after removing poly + magnitude terms, rescaled to the target's field angle (x theta_t/theta), and bin in
BP-RP.  Combine frames after normalising each frame by its own k_rad (to k_ref = mean), giving the colour curve
  s(C) [arcsec, radial, + = towards boresight] at the reference lateral-colour strength.
Also fits the field-angle dependence (colour term proportional to theta^a) to test the scaling.
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from loclib import R, to_json
C0 = 0.786
PXDEG = 21.0 / 3600.
FLOOR = 0.01


def frame(f):
    z = np.load(f); a = z['stars']
    a = a[(a['chi2r'] < 20) & (a['ruwe'] < 1.4) & (a['peak'] < 1.5e5) & np.isfinite(a['C'])]
    xt, yt, xb, yb = float(z['xt']), float(z['yt']), float(z['xb']), float(z['yb'])
    ux = xb - a['xw']; uy = yb - a['yw']; n = np.hypot(ux, uy)
    th = n * PXDEG; tht = np.hypot(xb - xt, yb - yt) * PXDEG
    ux /= n; uy /= n; tx, ty = uy, -ux
    sc = th / tht
    x = (a['xw'] - xt) / 400.; y = (a['yw'] - yt) / 400.
    dc = a['C'] - C0; dg = a['G'] - 12
    N = len(a)
    poly = [np.ones(N), x, y, x * x, x * y, y * y, x ** 3, y ** 3, x * x * y, x * y * y]
    Z = np.zeros(N)
    cols = [np.r_[p, Z] for p in poly] + [np.r_[Z, p] for p in poly] + [np.r_[dg, Z], np.r_[Z, dg],
            np.r_[dc * sc * ux, dc * sc * uy], np.r_[dc * sc * tx, dc * sc * ty]]
    A = np.column_stack(cols)
    b = np.r_[a['xf'] - a['xw'], a['yf'] - a['yw']]
    e = np.r_[np.hypot(a['ex'], FLOOR), np.hypot(a['ey'], FLOOR)]
    use = np.r_[a['C'] < 1.7, a['C'] < 1.7]
    ok = use.copy()
    for _ in range(8):
        w = 1 / e[ok]
        p, *_ = np.linalg.lstsq(A[ok] * w[:, None], b[ok] * w, rcond=None)
        r = (b - A @ p) / e
        s = 1.4826 * np.median(np.abs(r[ok]))
        good = np.abs(r) < 4 * s
        good = good[:N] & good[N:]
        new = use & np.r_[good, good]
        if np.array_equal(new, ok):
            break
        ok = new
    npoly = 2 * len(poly) + 2
    rr = b - A[:, :npoly] @ p[:npoly]
    rad = (rr[:N] * ux + rr[N:] * uy) / sc            # rescaled to target field angle
    tan = (rr[:N] * tx + rr[N:] * ty) / sc
    # outlier flag for all stars (incl. red): use robust scatter of blue stars around linear model
    lin = p[npoly] * dc
    res_rad = rad - lin
    sblue = 1.4826 * np.median(np.abs(res_rad[a['C'] < 1.7] - np.median(res_rad[a['C'] < 1.7])))
    keep = np.abs(res_rad - np.median(res_rad)) < 6 * max(sblue, 0.005) + np.abs(lin) * 0 + 0.2   # loose: keep red stars
    # field-angle exponent test: fit colour term ~ theta^a on blue stars
    return dict(file=os.path.basename(f), sec=int(z['sec']), k_rad=float(p[npoly]), k_tan=float(p[npoly + 1]), theta_t=float(tht),
                C=a['C'], G=a['G'], rad=rad, tan=tan, keep=keep, sblue=float(sblue), th=th)


frames = [frame(f) for f in sorted(glob.glob(f'{R}/results/x10w_*.npz'))]
kref = np.mean([fr['k_rad'] for fr in frames])
print('per-frame k_rad (mpx/mag, at target field angle):', [(fr['sec'], round(fr['k_rad'] * 1e3, 1), round(fr['k_tan'] * 1e3, 1), round(fr['theta_t'], 2)) for fr in frames])
bins = [(0.3, 0.6), (0.6, 0.75), (0.75, 0.9), (0.9, 1.05), (1.05, 1.25), (1.25, 1.5), (1.5, 1.75), (1.75, 2.0), (2.0, 2.3),
        (2.3, 2.6), (2.6, 3.0), (3.0, 3.6), (3.6, 4.5)]
allC = np.concatenate([fr['C'] for fr in frames]); allr = np.concatenate([fr['rad'] * kref / fr['k_rad'] for fr in frames])
allt = np.concatenate([fr['tan'] * kref / fr['k_rad'] for fr in frames]); allk = np.concatenate([fr['keep'] for fr in frames])
allsec = np.concatenate([np.full(len(fr['C']), fr['sec']) for fr in frames])
curve = []
print(f'\ncombined colour curve (radial, mpx at k_ref = {kref*1e3:.1f} mpx/mag; frames normalised; 2 FFIs per sector => pairs not independent):')
for lo, hi in bins:
    s = allk & (allC >= lo) & (allC < hi)
    if s.sum() < 3:
        continue
    # error: per-sector medians scatter (sectors are independent star samples only partially) -> use MAD/sqrt(n_unique/2)
    med = np.median(allr[s]); mad = 1.4826 * np.median(np.abs(allr[s] - med))
    n_eff = max(s.sum() / 2, 1)
    curve.append(dict(lo=lo, hi=hi, n=int(s.sum()), C_med=float(np.median(allC[s])), rad_mpx=float(med * 1e3),
                      err_mpx=float(1.2533 * mad / np.sqrt(n_eff) * 1e3), tan_mpx=float(np.median(allt[s]) * 1e3),
                      lin_pred_mpx=float(kref * (np.median(allC[s]) - C0) * 1e3)))
    c = curve[-1]
    print(f"  C {lo:.2f}-{hi:.2f}: n={c['n']:5d} med C {c['C_med']:.2f}  rad {c['rad_mpx']:+7.1f} +- {c['err_mpx']:4.1f} mpx "
          f"(linear pred {c['lin_pred_mpx']:+7.1f})  tan {c['tan_mpx']:+5.1f}")
to_json(dict(kref=kref, frames=[{k: v for k, v in fr.items() if k in ('file', 'sec', 'k_rad', 'k_tan', 'theta_t', 'sblue')} for fr in frames],
             curve=curve), f'{R}/results/x12_red_end.json')
np.savez_compressed(f'{R}/results/x12_points.npz', C=allC, rad=allr, tan=allt, keep=allk, sec=allsec)
