# Referee: independent L3 profile of the TESS eclipse shape (own model, e=0), 200-s sectors S68, S95, S105, S106
import sys, time, json, pickle
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all')
import numpy as np
from scipy.optimize import least_squares
from rebm import pair_flux
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d=pickle.load(open(B+'lcmodel/v6/results/data_default.pkl','rb'))
keys=['SPOC68','TC95','TC105','TC106']
sel=np.isin(d['key'],keys)
t=d['t'][sel]; y=d['y'][sel]; s=d['sig'][sel]; dt=d['dt'][sel]; key=d['key'][sel]; ev=d['ev'][sel]
# event boundaries (contiguous)
o=np.argsort(ev,kind='stable'); t,y,s,dt,key,ev=[a[o] for a in (t,y,s,dt,key,ev)]
_,starts=np.unique(ev,return_index=True)
ksc=np.array([keys.index(k) for k in key])
P=1.30367228; T0=3165.670416
NS=3; texp=200/86400.
off=(np.arange(NS)+0.5)/NS-0.5
ts=(t[:,None]+off[None,:]*texp).ravel()
w=1/s**2
def ld(q1,q2):
    return (2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2))
def model_pair(p):
    r1,k,cosi,J,q1a,q2a,q1b,q2b,dT0,dts=p[:10]
    F=pair_flux(ts,T0+dT0,P,r1,k,cosi,J,ld(q1a,q2a),ld(q1b,q2b),dts)
    return F.reshape(-1,NS).mean(1)
def resid(p,f):
    sc=np.array([1.0,p[10],p[11],p[12]])[ksc]
    m=1+sc*f*(model_pair(p)-1)
    X=np.vstack([m,m*dt,m*dt**2]).T
    # per-event weighted LS for baseline coefficients
    A=np.einsum('ni,nj,n->nij',X,X,w); bvec=np.einsum('ni,n,n->ni',X,y,w)
    As=np.add.reduceat(A,starts,axis=0); bs=np.add.reduceat(bvec,starts,axis=0)
    c=np.linalg.solve(As,bs[...,None])[...,0]
    cfull=np.repeat(c,np.diff(np.append(starts,len(y))),axis=0)
    mod=(X*cfull).sum(1)
    return (y-mod)/s
lo=[0.06,0.08,0.0,0.05,0,0,0,0,-0.003,-0.003,0.8,0.8,0.8]
hi=[0.25,1.3,0.35,1.5,1,1,1,1,0.003,0.003,1.2,1.2,1.2]
def fit(f,x0):
    r=least_squares(resid,x0,bounds=(lo,hi),args=(f,),x_scale='jac',max_nfev=600,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    return r
if __name__=='__main__':
    fgrid=[float(a) for a in sys.argv[1].split(',')]
    starts_x=[
      [0.112,0.44,0.012,0.52,0.36,0.3,0.36,0.3,0,0,1,1,1],
      [0.13,0.30,0.06,0.45,0.36,0.3,0.36,0.3,0,0,1,1,1],
      [0.16,0.20,0.10,0.45,0.36,0.3,0.36,0.3,0,0,1,1,1],
      [0.11,0.60,0.03,0.55,0.25,0.3,0.25,0.3,0,0,1,1,1]]
    out=[]
    for f in fgrid:
        best=None
        for x0 in starts_x:
            t1=time.time(); r=fit(f,np.array(x0))
            c2=float((r.fun**2).sum())
            print(f'f={f:.4f} start k0={x0[1]} chi2={c2:.2f} k={r.x[1]:.3f} r1={r.x[0]:.4f} b={r.x[2]/r.x[0]:.3f} J={r.x[3]:.3f} {time.time()-t1:.0f}s',flush=True)
            if best is None or c2<best[0]: best=(c2,r.x.tolist())
        out.append(dict(f=f,chi2=best[0],x=best[1]))
        json.dump(dict(N=len(y),res=out),open(B+'referee_all/r04_profile_%s.json'%sys.argv[2],'w'),indent=1)
