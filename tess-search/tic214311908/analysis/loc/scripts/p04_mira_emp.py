"""Does PRF-model mismatch of the Mira bias the rejection of the Mira hypothesis?  Fit the real DIFF/AMP images with
the Mira hypothesis modelled by the Mira's EMPIRICAL image (median image minus all other stars and background,
divided by the Mira's fitted flux; r01_sim.empirical_mira) instead of the PRF model, and also with the PRF placed at
the empirical image's own centroid.  Compare chi2 with the M3V (PRF) per sector and combined.
Output: out/mira_emp_fit.json"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from scipy.optimize import least_squares


def empirical_mira(M):
    S = M.S; g = stars(); xy = S.star_xy(g)
    im = int(np.where(g.source_id == IDS['Mira'])[0][0])
    O = np.median(S.flux[S.good], 0)
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mx, my = M.mira_xy
    mask = np.hypot(xx - mx, yy - my) <= 3.5
    fixed = np.zeros_like(O); ref = None
    for j in np.where(g['T'].values < 18)[0]:
        if j == im: continue
        img = M.psf(*xy[j])
        if g.source_id.values[j] == IDS['ref']: ref = img
        else: fixed += g.F.values[j] * img
    A = np.array([M.psf(mx, my), ref, np.ones_like(O)])
    w = 1 / np.sqrt(np.abs(O[mask]) + (0.02 * O[mask]) ** 2)
    c, *_ = np.linalg.lstsq((A[:, mask].T) * w[:, None], (O - fixed)[mask] * w, rcond=None)
    E = O - fixed - c[1] * ref - c[2]
    E = np.where(np.hypot(xx - mx, yy - my) <= 4.5, E, M.psf(mx, my) * c[0])
    return E / c[0], c[0]


def fit1(M, d, col):
    w = M.w; X = np.column_stack([col, np.ones_like(col)]) * np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(X, d * np.sqrt(w), rcond=None)
    r = X @ c - d * np.sqrt(w); return float(r @ r), c[0]


out = {}
for meth in ('DIFF', 'AMP'):
    tot = dict(M3V=0.0, Mira_prf=0.0, Mira_emp=0.0, pairA=0.0)
    for sec in SECTORS:
        M = SectorModel(sec, meth, grid_half=5.0, grid_step=5.0)
        E, fM = empirical_mira(M)
        e = E.ravel()[M.mask]
        c2, a, b, sa = M.fit(M.d_real, 'H'); h = M.hyp
        i3 = int(np.where(h.name == 'M3V')[0][0]); im = int(np.where(h.name == 'Mira')[0][0]); ip = int(np.where(h.name == 'pairA')[0][0])
        ce, ae = fit1(M, M.d_real, e)
        # mismatch size: rms of (empirical - PRF) relative to the PRF peak within 2 px
        pr = M.H[im]
        mis = float(np.sqrt(np.mean((e - pr)[pr > 0.02 * pr.max()] ** 2)) / pr.max())
        out[f'{meth}_S{sec}'] = dict(chi2_M3V=float(c2[i3]), chi2_Mira_prf=float(c2[im]), chi2_Mira_emp=ce, chi2_pairA=float(c2[ip]),
                                     amp_Mira_emp=float(ae), mira_flux=float(fM), emp_minus_prf_rms_over_peak=mis)
        tot['M3V'] += c2[i3]; tot['Mira_prf'] += c2[im]; tot['Mira_emp'] += ce; tot['pairA'] += c2[ip]
        print(f'{meth} S{sec}: chi2 - chi2(M3V): Mira(PRF) {c2[im]-c2[i3]:+6.1f}  Mira(empirical PSF) {ce-c2[i3]:+6.1f}  pairA {c2[ip]-c2[i3]:+6.1f};  '
              f'empirical-PRF rms/peak {mis*100:.1f}%; amp(Mira emp) {ae:.1f}', flush=True)
    out[f'{meth}_combined'] = {k: float(v - tot['M3V']) for k, v in tot.items()}
    print(f'{meth} combined chi2 - chi2(M3V):', {k: round(v, 1) for k, v in out[f"{meth}_combined"].items()})
to_json(out, f'{BASE}/out/mira_emp_fit.json')
