"""r10: sensitivity of the referee-2 baseline FPP to the fixed transit epoch: REF grid built with T0 = 3840.9160 (the REF
free-fit T0, 2.7 min later than the phase-1 T0 3840.9141 used elsewhere here; the analyst used 3840.9155).
Also re-checks ln B(single/twin) for that grid. Output: out/r10_t0_check.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
import r05_fpp as F
from rlib import jdump
out = {}
for ds in ('ref_beta', 'ref_beta_T0ref'):
    for ecc in ('kipping', 'circ'):
        p = dict(F.BASE); p.update(ds=ds, ecc=ecc); r = F.evaluate(p)
        tab = F.build(ds, ecc, 0.0)
        zs = np.mean(F.Zof(tab, 'N', 1.0, F.rng.normal(9.42, 0.88, 300), 0.381)); zt = np.mean(F.Zof(tab, 'N', 0.5, F.rng.normal(15.45, 1.45, 300), 0.289))
        out[f'{ds}|{ecc}'] = dict(FPP=r['FPP'], FPP_strict=r['FPP_strict'], PTP=r['post']['PTP_companion'], lnB_single_twin=float(np.log(zs/zt)))
        print(ds, ecc, out[f'{ds}|{ecc}'], flush=True)
p = dict(F.BASE); p.update(ds='ref_beta_T0ref', cmd='gauss', sig_s=0.30, xobs=-0.753); r = F.evaluate(p)
out['ref_beta_T0ref|kipping|analyst CMD'] = dict(FPP=r['FPP'], FPP_strict=r['FPP_strict'], PTP=r['post']['PTP_companion'])
print(out['ref_beta_T0ref|kipping|analyst CMD'])
jdump(out, 'out/r10_t0_check.json')
