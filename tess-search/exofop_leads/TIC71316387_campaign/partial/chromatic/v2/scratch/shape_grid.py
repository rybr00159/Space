import sys, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v2/scripts')
import lib2 as L
exec(open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v2/scripts/tess_ref2.py').read().split("def shp_from")[0])
prep = {}
for sec in data:
    d = data[sec]; k, sg = windows(d['t'], (0.0, 0.5)); prep[sec] = (d['t'][k], d['f'][k], sg, d['exptime'])
# fixed sigma per sector from out-of-eclipse residual of a quick fit
sig = {}
for s in prep:
    t,f,sg,ex = prep[s]
    oot = (L.dist_to_phase(t,0)>0.06)&(L.dist_to_phase(t,0.5)>0.06)
    a,ae,r,m = L.lsq(t[oot],f[oot],np.ones(oot.sum()),sg[oot],[],1) if False else (None,None,None,None)
    # simple: robust sd of point-to-point differences / sqrt2
    sig[s] = L.mad_std(np.diff(f[oot]))/np.sqrt(2)
print({k:round(v*1e6) for k,v in sig.items()})
def chi(par, which, secs):
    T,p,q = par
    shp = {'d':dict(T=0.064,p=1.26,q=1.6,phc=0.0),'s':dict(T=0.066,p=1.2,q=1.5,phc=0.5)}
    shp[which] = dict(T=T,p=p,q=q,phc=0.0 if which=='d' else 0.5)
    tot=0; amps=[]
    for s in secs:
        t,f,sg,ex = prep[s]
        a,ae,r = L.measure(t,f,np.ones_like(t),sg,shp,ex,1)
        tot += np.sum(r**2)/sig[s]**2; amps.append(a[0] if which=='d' else a[1])
    return tot, np.mean(amps)
secs=(28,68,95,105,106)
for which in ('d','s'):
    best=(1e30,)
    for p in (0.8,1.0,1.2,1.4,1.7,2.0,2.5,3.0):
        for q in (0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.0):
            row=[]
            for T in np.arange(0.054,0.080,0.002):
                c,A = chi((T,p,q),which,secs); row.append((c,T,A))
            c,T,A = min(row)
            if c<best[0]: best=(c,p,q,T,A)
            print(which,'p=%.1f q=%.1f bestT=%.3f (%.2f h) chi2=%.1f A=%.0f'%(p,q,T,T*24,c,A*1e6),flush=True)
    print('BEST',which,best)
