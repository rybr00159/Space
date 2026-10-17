"""Cube-level injection check of the linear-response shortcut used in r01_sim.py: inject a 15 e/s transit (template,
fake ephemeris T0 + 0.37 P) on the M3V's apparent position directly into the TESScut flux cube, run the full image
pipeline (background plane from the faintest pixels, DIFF and AMP), and compare with (same windows without injection)
+ A * PRF.  Also with the PRF placed on the Mira's apparent position and on a star 60" away."""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from l1img import Data
out = {}
for sec in (39, 104):
    S = Sector(sec); psf = S.prf()
    M0 = SectorModel(sec, 'DIFF', S=S)
    for nm in ('M3V', 'Mira', '5925373307210519552'):
        i = np.where(M0.hyp.name == nm)[0]
        if len(i) == 0:
            continue
        x, y = M0.hyp.loc[i[0], ['x', 'y']].values.astype(float)
        img = psf(x, y)
        A = 15.0
        tcs = real_tcs(S.t[S.good]) + 0.37 * P
        tpl = np.zeros_like(S.t); ok = np.isfinite(S.t)
        tpl[ok] = template(S.t[ok], tcs, S.exptime)
        fl = S.flux - A * tpl[:, None, None] * img[None]
        D0 = Data(S); D1 = Data(S, flux=fl)
        for meth in ('DIFF', 'AMP'):
            res = []
            for tc in tcs:
                if meth == 'DIFF':
                    p0, p1 = D0.diff_window(tc), D1.diff_window(tc)
                    if p0 is None or p1 is None: continue
                    diff = p1['D'] - p0['D']
                else:
                    p0, p1 = D0.amp_window(tc), D1.amp_window(tc)
                    if p0 is None or p1 is None: continue
                    diff = -(p1['num'] - p0['num']) / p1['den']
                res.append(diff - A * img.ravel())
            if not res: continue
            res = np.array(res)
            out[f'S{sec}_{meth}_{nm}'] = dict(max_abs_dev=float(np.abs(res).max()), rel_to_peak=float(np.abs(res).max() / (A * img.max())),
                                              sum_dev=float(res.sum(1).mean()), n=len(res))
            print(f'S{sec} {meth} {nm:20s}: windows {len(res)}; max |response - A*PRF| = {np.abs(res).max():.4f} e/s '
                  f'({np.abs(res).max() / (A * img.max()) * 100:.2f}% of peak); summed deviation {res.sum(1).mean():+.4f} e/s')
to_json(out, f'{BASE}/out/cubecheck.json')
