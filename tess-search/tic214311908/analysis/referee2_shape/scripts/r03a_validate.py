"""r03a: validate the rwd overlap model against batman (uniform LD) for k = 0.5, 3, 5, 30 and check exposure averaging.
Output: out/validate_rwd.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from rwd import *
D = load('ref', 'white'); lc = LC(D); sup = Super(lc)
res = {}
for k, aR, b in [(0.5, 20., 0.3), (3.0, 60., 1.0), (5.0, 100., 2.0), (30.0, 1000., 10.)]:
    RA = 1.0; RB = k; a = aR*RA; cosi = b/aR
    m1 = wd_occultation(sup, T0REF, PREF, a, cosi, RA, RB, 1.0)
    m2 = lc.transit(T0REF, PREF, k, aR, b, 0., 0.)   # batman: quadratic with u=0 -> uniform
    res[str(k)] = dict(maxdiff=float(np.max(np.abs(m1 - m2))), mindepth=float(1 - m1.min()))
    print(k, res[str(k)])
json.dump(res, open(f'{RD}/out/validate_rwd.json', 'w'), indent=1)
