"""Robustness of the real-data attribution to analysis choices (no injections): combined chi2 differences of the main
hypotheses relative to the M3V and the model-based posterior (kappa = 1, equal priors over the 93 catalogued stars).
Variants: fit radius, DIFF out-of-transit window width, AMP spline knot spacing, no background subtraction, lateral-
colour curve scaled x0.7 / x1.3, WCS offset shifted by +-0.03 px, the Mira at its Gaia position or at Gaia +
extrapolated colour curve instead of its measured apparent position, and dropping each sector.
Output: out/variants.json"""
import sys, copy
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from l1img import Data
SIDS = {}
MAINS = ['M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208']
Ss = {s: Sector(s) for s in SECTORS}
Ds = {s: Data(Ss[s]) for s in SECTORS}
def run(label, meth='DIFF', data_kw=None, model_kw=None, colscale=1.0, wcs_shift=None, drop=None):
    chis = []; frees = []; names = None
    for s in SECTORS:
        if drop is not None and s == drop: continue
        S = Ss[s]; cs0, wo0 = S.colscale, S.wcs_off.copy()
        S.colscale = cs0 * colscale
        if wcs_shift is not None: S.wcs_off = wo0 + np.array(wcs_shift)
        Dd = Ds[s] if data_kw is None else Data(S, **data_kw)
        M = SectorModel(s, meth, S=S, Dd=Dd, grid_half=40.0, **(model_kw or {}))
        c2, a, b, sa = M.fit(M.d_real, 'H'); cg, *_ = M.fit(M.d_real, 'G')
        chis.append(c2); frees.append(cg); names = M.hyp.name.values
        S.colscale, S.wcs_off = cs0, wo0
    C = np.sum(chis, 0); G = np.sum(frees, 0)
    idx = {n: i for i, n in enumerate(names)}
    Pp = np.exp(-0.5 * (C - C.min())); Pp /= Pp.sum()
    r = {n: float(C[idx[n]] - C[idx['M3V']]) for n in MAINS}
    r['best'] = str(names[int(np.argmin(C))]); r['free_minus_M3V'] = float(G.min() - C[idx['M3V']])
    r['P_M3V'] = float(Pp[idx['M3V']]); r['P_pair'] = float(Pp[idx['pairA']] + Pp[idx['pairB']]); r['P_Mira'] = float(Pp[idx['Mira']])
    r['P_faint_near_M3V'] = float(Pp[idx['5925373272827260544']] + Pp[idx['5925373272827262208']])
    print(f"{label:28s} {meth}: chi2-chi2(M3V): " + ', '.join(f'{n[-6:]} {r[n]:6.1f}' for n in MAINS[1:]) +
          f"; best {r['best'][-8:]}; free-fit gain {r['free_minus_M3V']:.1f}; P(M3V) {r['P_M3V']:.4f} P(pair) {r['P_pair']:.2g} P(Mira) {r['P_Mira']:.2g} P(faint) {r['P_faint_near_M3V']:.2g}", flush=True)
    return r
out = {}
for meth in ('DIFF', 'AMP'):
    out[f'base_{meth}'] = run('base', meth)
    out[f'colour x0.7_{meth}'] = run('colour curve x0.7', meth, colscale=0.7)
    out[f'colour x1.3_{meth}'] = run('colour curve x1.3', meth, colscale=1.3)
    out[f'mira gaia_{meth}'] = run('Mira at Gaia position', meth, model_kw=dict(mira_mode='gaia'))
    out[f'mira colourcurve_{meth}'] = run('Mira at Gaia+colour curve', meth, model_kw=dict(mira_mode='colourcurve'))
    for rf in (3.5, 4.5, 7.0):
        out[f'rfit {rf}_{meth}'] = run(f'fit radius {rf} px', meth, model_kw=dict(r_fit=rf))
    for ws in ([0.03, 0], [-0.03, 0], [0, 0.03], [0, -0.03]):
        out[f'wcs {ws}_{meth}'] = run(f'WCS shift {ws} px', meth, wcs_shift=ws)
    for d in SECTORS:
        out[f'drop S{d}_{meth}'] = run(f'drop S{d}', meth, drop=d)
for wo in (0.15, 0.40):
    out[f'w_oot {wo}_DIFF'] = run(f'OOT window {wo} d', 'DIFF', model_kw=dict(w_oot=wo))
out['no background subtraction_DIFF'] = run('no background subtraction', 'DIFF', data_kw=dict(bkg=False))
for kn in (0.3, 1.0):
    out[f'knot {kn}_AMP'] = run(f'spline knots {kn} d', 'AMP', data_kw=dict(knot=kn))
to_json(out, f'{BASE}/out/variants.json')
