"""r08: extra what-if checks on the referee-2 baseline: two-band depth-ratio precision 10% (instead of 5%), with and
without resolving pair A; the same on the INV light curve; f_BD = 0 to show the BD share. Output: out/r08_whatif.json"""
import sys, json
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
import r05_fpp as F
from rlib import jdump
HRI = ([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0], [2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0])
out = {}
for ds in ('ref_beta', 'inv_beta'):
    for name, kw in (('none', {}), ('two-band 10% + pair resolved', dict(chroma=0.10, lr_pairA=1e-8)),
                     ('two-band 5% + pair resolved', dict(chroma=0.05, lr_pairA=1e-8)), ('speckle', dict(hri=HRI)),
                     ('speckle + two-band 10%', dict(hri=HRI, chroma=0.10, lr_pairA=1e-8)),
                     ('speckle + two-band 10%, no BD class', dict(hri=HRI, chroma=0.10, lr_pairA=1e-8, f_BD=0.0))):
        p = dict(F.BASE); p.update(ds=ds, **kw); r = F.evaluate(p)
        out[f'{ds} | {name}'] = dict(FPP=r['FPP'], FPP_strict=r['FPP_strict'], NFPP=r['NFPP'], dominant=r['dominant_FP'],
                                     BD_total=r['post']['BD_on_M3V_single'] + r['post']['BD_on_primary_of_binary'] + r['post']['BDHEB_on_companion'])
        print(f'{ds} | {name:40s}', {k: (round(v, 5) if isinstance(v, float) else v) for k, v in out[f"{ds} | {name}"].items()}, flush=True)
jdump(out, 'out/r08_whatif.json')
