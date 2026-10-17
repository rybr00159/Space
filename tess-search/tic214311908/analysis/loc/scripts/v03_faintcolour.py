"""Sensitivity of the M3V-versus-(G = 17.5 star 3.2" away, Gaia DR3 5925373272827260544) discrimination to that
star's assumed colour (its Gaia BP-RP = 1.32 may be biased by the M3V 3.2" away): combined chi2 difference when its
apparent position is computed with BP-RP = 1.0, 1.32 (catalogue), 2.0, 2.76 (as red as the M3V); real data, both
methods, with and without S39.  Also the same for the G = 17.7 star 5.7" away (..262208).
Output: out/faint_colour.json"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
g = stars()
out = {}
for meth in ('DIFF', 'AMP'):
    Ms = {s: SectorModel(s, meth, grid_half=2.0, grid_step=2.0) for s in SECTORS}
    for sid in (5925373272827260544, 5925373272827262208):
        for C in (1.0, 1.32 if sid == 5925373272827260544 else 0.96, 2.0, 2.76):
            tot = 0.0; tot_no39 = 0.0
            for s, M in Ms.items():
                r = g[g.source_id == sid].copy(); r['C'] = C
                xy = M.S.star_xy(r)[0]
                p = M.psf(*xy).ravel()[M.mask]
                w = M.w; X = np.column_stack([p, np.ones_like(p)]) * np.sqrt(w)[:, None]
                c, *_ = np.linalg.lstsq(X, M.d_real * np.sqrt(w), rcond=None); rr = X @ c - M.d_real * np.sqrt(w)
                c2, *_ = M.fit(M.d_real, 'H'); i3 = int(np.where(M.hyp.name == 'M3V')[0][0])
                d = float(rr @ rr - c2[i3]); tot += d; tot_no39 += d if s != 39 else 0.0
            out[f'{meth}_{sid}_C{C}'] = dict(dchi2_vs_M3V=tot, dchi2_vs_M3V_without_S39=tot_no39)
            print(f'{meth} {sid} assumed BP-RP {C}: chi2 - chi2(M3V) = {tot:.1f} (without S39: {tot_no39:.1f})', flush=True)
to_json(out, f'{BASE}/out/faint_colour.json')
