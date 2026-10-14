"""VHS differential test using a matched saturated star as the PSF template.

N228 (Gaia DR3 near 330.9009 -42.7729; Ks 10.43, J 10.96) lies on the same VIRCAM chip 3.8' from the target with almost the
same brightness (target Ks 10.42, J 11.04), so its (equally saturated) image is the best available PSF template.
For each pawprint (J1, J2, Ks1, Ks2):
  * T = sky-subtracted N228 stamp; the template is shifted to the target with cubic-spline interpolation;
  * model target = a * T(x0, y0) + b [+ f * a * T(x0 + dx, y0 + dy)], pixels saturated in either star masked
    (raw > SATLIM), plus pixels where the template itself is non-linear;
  * scan the companion over 0.5-3.5"; report best dchi2, sep, PA, dm, and the residual image;
  * symmetric control: the same fit with roles swapped (N228 modelled with the target as template): a real companion
    of the target appears there as a NEGATIVE-flux feature at the mirrored position, so we scan allowing f < 0 too;
  * third star control: N169 (Ks 11.77) modelled with N228 as template.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares
from v3lib import V3, to_json
import warnings; warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

SATLIM = float(os.environ.get('H03_SATLIM', 24000)); GAIN, RN = 4.2, 25.0
H = 14; EP = 2012.64
POS = {'target': (330.887948 + 38.5 * (EP - 2016) / 3.6e6 / np.cos(np.radians(-42.8355)), -42.835511 - 3.13 * (EP - 2016) / 3.6e6)}
gw = Table.read(V3 + '/gaia/wide_2p7deg_G18p5.fits').to_pandas()
for name, (ra, de) in {'N228': (330.9009, -42.7729), 'N169': (330.8423, -42.8027)}.items():
    d = np.hypot((gw.ra - ra) * np.cos(np.radians(de)), gw.dec - de) * 3600
    r = gw.iloc[int(np.argmin(d))]
    POS[name] = (r.ra + np.nan_to_num(r.pmra) * (EP - 2016) / 3.6e6 / np.cos(np.radians(r.dec)), r.dec + np.nan_to_num(r.pmdec) * (EP - 2016) / 3.6e6)


def get(d, w, name, n):
    x, y = w.all_world2pix([POS[name]], 0)[0]
    ix, iy = int(round(x)), int(round(y))
    big = d[iy - 45:iy + 46, ix - 45:ix + 46]
    yy, xx = np.mgrid[-45:46, -45:46]
    sky = float(np.median(big[(np.hypot(xx, yy) > 28)]))
    raw = d[iy - n:iy + n + 1, ix - n:ix + n + 1].astype(float)
    # refine centre with a centroid of the outer core (r 3-8 px) to avoid the flat top
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    return raw, sky, x - ix, y - iy


def fit(tgt, tsky, tmp, msky, mask_extra, pix, Jn, allow_neg=False, rmin=0.5, rmax=3.5, step=0.1):
    st = tgt - tsky; T = tmp - msky
    yy, xx = np.mgrid[-H:H + 1, -H:H + 1]
    sig = np.sqrt(np.clip(tgt, 1, None) / GAIN + RN ** 2 + (0.0 * T) ** 2)
    base = (tgt < SATLIM) & (np.hypot(xx, yy) <= H) & mask_extra

    def model_shift(dx, dy):
        return ndshift(T, (dy, dx), order=3, mode='constant', cval=0.0)

    def res(p):
        m = model_shift(p[0], p[1])
        mk = base & (np.abs(m) > -1)     # placeholder (template saturation handled via mask_extra)
        A = np.column_stack([m[mk], np.ones(mk.sum())]) / sig[mk][:, None]
        c, *_ = np.linalg.lstsq(A, st[mk] / sig[mk], rcond=None)
        return A @ c - st[mk] / sig[mk]
    r = least_squares(res, [0.0, 0.0], diff_step=1e-3)
    p0 = r.x
    m0 = model_shift(*p0)
    mk = base
    w = 1 / sig[mk]; y = st[mk] * w
    A0 = np.column_stack([m0[mk] * w, w]); c0, *_ = np.linalg.lstsq(A0, y, rcond=None)
    chi0 = float(((A0 @ c0 - y) ** 2).sum())
    best = (0, None, None); grid = np.arange(-rmax, rmax + 1e-9, step)
    for a in grid:
        for b in grid:
            rr = np.hypot(a, b)
            if rr < rmin or rr > rmax:
                continue
            m1 = model_shift(p0[0] + a / pix, p0[1] + b / pix)
            A = np.column_stack([m0[mk] * w, m1[mk] * w, w])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            if c[1] < 0 and not allow_neg:
                continue
            d = chi0 - float(((A @ c - y) ** 2).sum())
            if d > best[0]:
                best = (d, (a, b), c[1] / c[0])
    so = Jn @ np.array(best[1]) if best[1] else np.zeros(2)
    resid = st - (c0[0] * m0 + c0[1])
    return dict(chi2r=chi0 / (mk.sum() - 4), dchi2=best[0], sep=float(np.hypot(*best[1])) if best[1] else None,
                pa=float(np.degrees(np.arctan2(so[0], so[1])) % 360), f=best[2],
                dm=(float(-2.5 * np.log10(best[2])) if best[2] and best[2] > 0 else None), resid=resid, sig=sig, mask=mk)


out = []
fig, axs = plt.subplots(4, 3, figsize=(12, 15))
for k, f in enumerate(sorted(glob.glob(V3 + '/imaging/big/vhs_*.fits.fz'))):
    h = fits.open(f); ext = [e for e in h if e.data is not None and e.data.shape[1] > 100][0]
    d = ext.data.astype(float); w = WCS(ext.header)
    band = 'J' if '14:49' in f or '14:52:33' in f else 'Ks'
    x0, y0 = w.all_world2pix([POS['target']], 0)[0]
    r0 = w.all_pix2world([[x0, y0], [x0 + 1, y0], [x0, y0 + 1]], 0); cc = np.cos(np.radians(r0[0, 1]))
    J = np.array([[(r0[1, 0] - r0[0, 0]) * cc * 3600, (r0[2, 0] - r0[0, 0]) * cc * 3600], [(r0[1, 1] - r0[0, 1]) * 3600, (r0[2, 1] - r0[0, 1]) * 3600]])
    pix = float(np.sqrt(abs(np.linalg.det(J)))); Jn = J / pix
    S = {n: get(d, w, n, H) for n in POS}
    rows = {}
    for tname, mname in (('target', 'N228'), ('N228', 'target'), ('N169', 'N228'), ('N169', 'target')):
        tgt, tsky, tfx, tfy = S[tname]; tmp, msky, mfx, mfy = S[mname]
        # align template's sub-pixel phase to the target's: shift template by (tfx - mfx, tfy - mfy)
        tmp_al = ndshift(tmp - msky, (tfy - mfy, tfx - mfx), order=3, mode='nearest') + msky
        mask_extra = (tmp < SATLIM) & (ndshift((tmp >= SATLIM).astype(float), (tfy - mfy, tfx - mfx), order=1) < 0.01)
        r = fit(tgt, tsky, tmp_al, msky, mask_extra, pix, Jn, allow_neg=(tname == 'N228'))
        rows[f'{tname}|{mname}'] = r
        out.append(dict(file=os.path.basename(f), band=band, star=tname, template=mname, chi2r=r['chi2r'], dchi2=r['dchi2'],
                        dchi2_over_chi2r=r['dchi2'] / r['chi2r'], sep=r['sep'], pa=r['pa'], f=r['f'], dm=r['dm']))
        print(f"{os.path.basename(f)[4:30]} {band:2s} {tname:6s} tmpl {mname:6s}: chi2r {r['chi2r']:.1f} best dchi2 {r['dchi2']:.0f} "
              f"(/chi2r {r['dchi2'] / r['chi2r']:.0f}) sep {r['sep']:.2f}\" PA {r['pa']:.0f} f {r['f']:+.4f} dm {r['dm']}", flush=True)
    for j, key in enumerate(['target|N228', 'N228|target', 'N169|N228']):
        r = rows[key]; z = np.where(r['mask'], r['resid'] / r['sig'], np.nan)
        ax = axs[k, j]; ax.imshow(z, origin='lower', cmap='RdBu_r', vmin=-15, vmax=15); ax.set_title(f'{band} {key} resid/sig', fontsize=9)
        # compass: E and N directions in image axes
        inv = np.linalg.inv(Jn)
        for v, lab in ((np.array([1.0, 0]), 'E'), (np.array([0, 1.0]), 'N')):
            p = inv @ v * 5
            ax.arrow(H, H, p[0], p[1], color='k', head_width=0.6); ax.text(H + p[0] * 1.2, H + p[1] * 1.2, lab)
plt.tight_layout(); plt.savefig(V3 + '/plots/h03_vhs_template_resid.png', dpi=70)
to_json(out, V3 + '/results/h03_vhs_template.json')
