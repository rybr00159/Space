"""s13b: the s13 fitting functions (planet free, WD occultation by the M3V, V-shape) applied to the real INV and
REF windows, so that the observed dchi2(WD - planet) and dchi2(V - planet) can be compared with the injection
distributions of s13. Output: out/inject_real.json"""
import sys, json, numpy as np
sys.argv = ['x', 'planet']
src = open('tic214311908/analysis/shape/scripts/s13_inject.py').read()
src = src[:src.index("rows = []")]
exec(src)
out = {}
for ds in ['inv', 'ref']:
    lc = LC(f'{SD}/data/{ds}_win.npz')
    fp = fit_planet(lc, T0r); fw = fit_wd(lc, T0r); fv = fit_vshape(lc, T0r)
    out[ds] = dict(planet=fp, wd=fw, vshape=fv, dchi2_wd_minus_planet=fw['chi2'] - fp['chi2'], dchi2_v_minus_planet=fv['chi2'] - fp['chi2'],
                   dchi2_transit=float(lc.chi2(np.ones(lc.n)) - fp['chi2']))
    print(ds, out[ds], flush=True)
json.dump(out, open(f'{SD}/out/inject_real.json', 'w'), indent=1)
