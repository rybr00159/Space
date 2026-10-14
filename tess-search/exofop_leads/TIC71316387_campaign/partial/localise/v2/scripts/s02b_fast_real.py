import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastfit import *
from loclib import R, SECTORS, to_json
out = {}
for sec in [int(a) for a in sys.argv[1:]] or SECTORS:
    t0 = time.time()
    F = Fast(sec); ctx = F.ctx
    k = ctx.kind
    o = {}
    ref = json.load(open(f'{R}/results/real_s{sec}_prf.json'))
    for name in ['deep', 'shallow']:
        idx = np.where(k == name)[0]
        Ds, ss, Bs = ctx.stack(idx, np.ones(len(idx)))
        rB = F.fit_direct(Bs, ss); rD = F.fit_diff(Ds, ss)
        sky, pix = F.offset_sky(rD, rB)
        cs = F.jac @ rD['cov'] @ F.jac.T
        o[name] = dict(dE=sky[0], dN=sky[1], dx=pix[0], dy=pix[1], amp=rD['amp'], F=rB['F'], n=len(idx),
                       formal_sky=np.sqrt(np.diag(cs)).tolist(), edge=rD['edge'])
        print(f"S{sec} {name:7s} fast ({sky[0]:+.2f},{sky[1]:+.2f}) amp {rD['amp']:.1f} | lsq ({ref[name]['dE']:+.2f},{ref[name]['dN']:+.2f}) amp {ref[name]['amp']:.1f}  [{time.time()-t0:.0f}s]", flush=True)
    out[sec] = o
to_json(out, f'{R}/results/real_fast.json')
