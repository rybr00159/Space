"""s12b: null model for the background stars' photometry: a single main-sequence star (Mamajek mass 0.08-1.4,
distance, E(B-V) in [0, 0.30]) fitted to the DECaPS2 g,r,i,z mean AB magnitudes (sigma 0.05) and the Gaia parallax,
exactly as in s12. Gives the photometric chi2 that an EB scenario must be compared with.
Output: out/bg_single.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from eem import ms
from specs import absmag, EXT
from scipy.optimize import differential_evolution
SD = 'tic214311908/analysis/shape'
OBS = dict(bg3=dict(mag=dict(g=18.164, r=17.352, i=17.100, z=16.900), plx=(0.655+0.03, 0.110)),
           bg5=dict(mag=dict(g=18.090, r=17.589, i=17.407, z=17.304), plx=(0.283+0.03, 0.115)))
out = {}
for src, o in OBS.items():
    def c2(x, detail=False):
        m, logd, ebv = x; d = 10**logd
        T = ms(m, 'Teff'); R = ms(m, 'R_Rsun')
        mags = {b: absmag(T, R, b) + 5*np.log10(d/10) + EXT[b]*ebv for b in o['mag']}
        cp = sum(((mags[b] - o['mag'][b])/0.05)**2 for b in o['mag']); cx = ((1000/d - o['plx'][0])/o['plx'][1])**2
        if detail: return dict(chi2_phot=float(cp), chi2_plx=float(cx), mass=float(m), Teff=float(T), R=float(R), d_pc=float(d), ebv=float(ebv), mags={b: float(v) for b, v in mags.items()})
        return cp + cx
    r = differential_evolution(c2, [(0.08, 1.4), (2.7, 4.3), (0, 0.3)], seed=1, popsize=30, tol=1e-10, polish=True)
    out[src] = c2(r.x, True); out[src]['total'] = float(r.fun)
    print(src, out[src])
json.dump(out, open(f'{SD}/out/bg_single.json', 'w'), indent=1)
