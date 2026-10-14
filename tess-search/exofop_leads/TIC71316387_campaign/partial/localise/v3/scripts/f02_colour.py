"""v3 lateral-colour analysis of the f01 field-star tables.

Residuals (fitted - Gaia/WCS) are converted to sky arcsec with the local Jacobian, then modelled per FFI as
   r = poly2(x, y) [per sky component]  +  m * (G - 12)  +  (C - C0) * s(theta) * (k_rad * u_rad + k_tan * u_tan)
with C = BP-RP, C0 = 0.786 (target), u_rad = unit sky vector from the star TOWARDS the camera boresight (tess-point),
s(theta) = theta / theta_target (lateral colour ~ field angle) or 1 (variant).  k_rad < 0: redder stars appear
displaced AWAY from the boresight.  Robust (4-sigma MAD clipping), bootstrap over stars for errors.
Then the binned colour curve (radial residual after removing poly + magnitude terms, rescaled to the target field angle)
is built per FFI and combined, to test linearity at the red end.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v3lib import V3, to_json

C0 = 0.786
SCALE = os.environ.get('F02_SCALE', 'theta')
CMAXFIT = float(os.environ.get('F02_CMAXFIT', 9.0))
TAG = os.environ.get('F01_TAG', '')
FLOOR = 0.15   # arcsec error floor per star and component


def load(f):
    z = np.load(f)
    a = z['stars']
    ok = (a['chi2r'] < 20) & (a['ruwe'] < 1.4) & (a['peak'] < 1.0e5) & (a['flux'] > 0) & (a['ex'] < 0.15) & (a['ey'] < 0.15)
    a = a[ok]
    J = z['J']
    dpx = np.vstack([a['xf'] - a['xw'], a['yf'] - a['yw']])
    dsky = J @ dpx                              # arcsec (E, N)
    esky = np.sqrt((J ** 2) @ np.vstack([a['ex'] ** 2, a['ey'] ** 2]))
    pa = np.radians(a['pa_bs'])
    ur = np.vstack([np.sin(pa), np.cos(pa)])    # towards boresight (E, N)
    ut = np.vstack([ur[1], -ur[0]])
    th = a['theta']; tht = float(z['theta_t'])
    return dict(a=a, d=dsky, e=esky, ur=ur, ut=ut, th=th, tht=tht, sec=int(z['sec']), xt=float(z['xt']), yt=float(z['yt']),
                file=os.path.basename(f), pa_t=float(z['pa_t']))


def design(F, scale=SCALE):
    a = F['a']; N = len(a)
    x = (a['xw'] - F['xt']) / 400.0; y = (a['yw'] - F['yt']) / 400.0
    poly = [np.ones(N), x, y, x * x, x * y, y * y]
    Z = np.zeros(N)
    sc = F['th'] / F['tht'] if scale == 'theta' else np.ones(N)
    dc = (a['C'] - C0) * sc
    dg = a['G'] - 12.0
    cols = [np.r_[p, Z] for p in poly] + [np.r_[Z, p] for p in poly]
    cols += [np.r_[dg, Z], np.r_[Z, dg]]
    cols += [np.r_[dc * F['ur'][0], dc * F['ur'][1]], np.r_[dc * F['ut'][0], dc * F['ut'][1]]]
    A = np.column_stack(cols)
    b = np.r_[F['d'][0], F['d'][1]]
    e = np.sqrt(np.r_[F['e'][0], F['e'][1]] ** 2 + FLOOR ** 2)
    return A, b, e, sc


def rfit(A, b, e, use, n=None):
    n = len(b) // 2 if n is None else n
    ok = use.copy()
    for _ in range(10):
        w = 1 / e[ok]
        p, *_ = np.linalg.lstsq(A[ok] * w[:, None], b[ok] * w, rcond=None)
        r = (b - A @ p) / e
        s = 1.4826 * np.median(np.abs(r[ok]))
        g = np.abs(r) < 4 * s
        g = g[:n] & g[n:]
        new = use & np.r_[g, g]
        if np.array_equal(new, ok):
            break
        ok = new
    return p, ok, s


frames = [load(f) for f in sorted(glob.glob(f'{V3}/results/f01{TAG}_tess*.npz'))]
res = []
curves = []
rng = np.random.default_rng(3)
for F in frames:
    A, b, e, sc = design(F)
    N = len(F['a'])
    use = np.r_[F['a']['C'] < CMAXFIT, F['a']['C'] < CMAXFIT]
    p, ok, s = rfit(A, b, e, use)
    boots = []
    for _ in range(400):
        i = rng.integers(0, N, N); ii = np.r_[i, i + N]
        pb, _, _ = rfit(A[ii], b[ii], e[ii], use[ii])
        boots.append(pb[-4:])
    boots = np.array(boots)
    kr, kt = p[-2], p[-1]
    ekr, ekt = boots[:, -2].std(), boots[:, -1].std()
    rr = b - A[:, :-2] @ p[:-2]
    rad = (rr[:N] * F['ur'][0] + rr[N:] * F['ur'][1]) / sc
    tan = (rr[:N] * F['ut'][0] + rr[N:] * F['ut'][1]) / sc
    okk = ok[:N]
    res.append(dict(file=F['file'], sec=F['sec'], n=N, n_used=int(okk.sum()), n_red2=int(((F['a']['C'] > 2.0) & okk).sum()),
                    theta_t=F['tht'], pa_t=F['pa_t'], k_rad=kr, k_rad_err=ekr, k_tan=kt, k_tan_err=ekt,
                    m_E=p[12], m_N=p[13], rms=s))
    curves.append(dict(C=F['a']['C'], rad=rad, tan=tan, ok=okk, sec=F['sec'], G=F['a']['G'], plx=F['a']['plx']))
    print(f"S{F['sec']:3d} theta_t={F['tht']:.2f} n={okk.sum()}/{N} (C>2: {res[-1]['n_red2']}) k_rad={kr:+.3f}+-{ekr:.3f} \"/mag "
          f"k_tan={kt:+.3f}+-{ekt:.3f}  mag-term ({p[12]:+.3f},{p[13]:+.3f}) \"/mag  clip-scale {s:.2f}", flush=True)
k = np.array([r['k_rad'] for r in res]); ek = np.array([r['k_rad_err'] for r in res])
kw = (k / ek ** 2).sum() / (1 / ek ** 2).sum()
print(f'weighted mean k_rad = {kw:+.3f} +- {1 / np.sqrt((1 / ek ** 2).sum()):.3f} "/mag; chi2 = {(((k - kw) / ek) ** 2).sum():.1f} / {len(k) - 1}')
# binned combined colour curve (rescaled to target field angle; NOT normalised by per-frame k)
bins = [(-0.1, 0.55), (0.55, 0.7), (0.7, 0.85), (0.85, 1.0), (1.0, 1.15), (1.15, 1.35), (1.35, 1.6), (1.6, 1.9), (1.9, 2.3), (2.3, 2.7),
        (2.7, 3.2), (3.2, 4.5)]
allC = np.concatenate([c['C'] for c in curves]); allR = np.concatenate([c['rad'] for c in curves])
allT = np.concatenate([c['tan'] for c in curves]); allok = np.concatenate([c['ok'] for c in curves])
allplx = np.concatenate([c['plx'] for c in curves]); allG = np.concatenate([c['G'] for c in curves])
curve = []
print('\ncombined radial colour curve (arcsec, + = towards boresight), rescaled to target field angle:')
for lo, hi in bins:
    m = allok & (allC >= lo) & (allC < hi)
    if m.sum() < 3:
        continue
    med = np.median(allR[m]); mad = 1.4826 * np.median(np.abs(allR[m] - med))
    err = 1.2533 * mad / np.sqrt(m.sum())
    # dwarfs (M_G > 5 with plx/err) vs others
    absG = allG + 5 * np.log10(np.clip(allplx, 1e-3, None) / 100.0)
    dw = m & (absG > 5.0) & (allplx > 1)
    curve.append(dict(lo=lo, hi=hi, n=int(m.sum()), C=float(np.median(allC[m])), rad=float(med), err=float(err),
                      tan=float(np.median(allT[m])), n_dwarf=int(dw.sum()),
                      rad_dwarf=float(np.median(allR[dw])) if dw.sum() >= 3 else None,
                      lin=float(kw * (np.median(allC[m]) - C0))))
    cc = curve[-1]
    print(f"  C {lo:5.2f}-{hi:4.2f} n={cc['n']:4d} medC={cc['C']:.2f} rad={cc['rad']:+.3f}+-{cc['err']:.3f}  (linear {cc['lin']:+.3f})  "
          f"tan={cc['tan']:+.3f}  dwarfs n={cc['n_dwarf']} rad={cc['rad_dwarf'] if cc['rad_dwarf'] is None else round(cc['rad_dwarf'], 3)}")
to_json(dict(scale=SCALE, frames=res, k_rad_wmean=kw, curve=curve), f'{V3}/results/f02_colour{TAG}_{SCALE}.json')
np.savez_compressed(f'{V3}/results/f02_points{TAG}.npz', C=allC, rad=allR, tan=allT, ok=allok,
                    sec=np.concatenate([np.full(len(c['C']), c['sec']) for c in curves]), G=allG, plx=allplx)
