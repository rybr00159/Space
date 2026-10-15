# Referee: combine my independent L3 profile with a soft main-sequence prior (two different MS relation sets)
import sys, json, glob
import numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib'); import ms
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all/'
a=[-0.642,-0.208,-8.43e-4,7.87e-3,1.42e-4,-2.13e-4]
P=1.30367228*86400; G=6.674e-8; rhosun=1.408
pts=[]
for f in sorted(glob.glob(B+'r04_profile_*.json')):
    for r in json.load(open(f))['res']: pts.append((r['f'],r['chi2'],r['x']))
pts.sort(); cmin=min(p[1] for p in pts)
# relation set A: Mann+2019 mass, Mann+2015 radius, P&M T-Ks colour; set B: lcmodel smoothed P&M table
mg=np.linspace(0.1,0.85,300); MKg=ms.MKs(mg); MTg=ms.MT(mg); o=np.argsort(MKg)
def setA(MK):
    x=MK-7.5; M=10**sum(c*x**i for i,c in enumerate(a)); R=1.9515-0.3520*MK+0.01680*MK**2
    MT=np.interp(MK,MKg[o],MTg[o]); return M,R,MT
def setB(M):
    return M,ms.radius(M),ms.MT(M)
gridA=np.linspace(4.6,9.3,400); gridB=np.linspace(0.1,0.85,400)
sig=dict(rho=0.21,k=0.10,J=0.10)
out={}
for name,fun,grid in (('A_Mann',setA,gridA),('B_PM_table',setB,gridB)):
    M,R,MT=[np.array(v) for v in fun(grid)]
    SB=10**(-0.4*MT)/R**2
    i1,i2=np.meshgrid(np.arange(len(grid)),np.arange(len(grid)),indexing='ij')
    ok=M[i2]<=M[i1]
    rows=[]
    for fEB,c2,x in pts:
        r1,k,cosi,J=x[:4]
        q=M[i2]/M[i1]
        rho_fit=3*np.pi/(G*P**2)/r1**3/(1+q)/rhosun
        rho_ms=M[i1]/R[i1]**3
        kms=R[i2]/R[i1]; Jms=SB[i2]/SB[i1]
        chi=(np.log(rho_fit/rho_ms)/sig['rho'])**2+(np.log(k/kms)/sig['k'])**2+(np.log(J/Jms)/sig['J'])**2
        chi=np.where(ok,chi,np.inf); j=np.unravel_index(np.argmin(chi),chi.shape)
        m1,m2=M[i1][j],M[i2][j]
        MTp=-2.5*np.log10(10**(-0.4*MT[i1][j])+10**(-0.4*MT[i2][j]))
        d=10**((11.358-2.5*np.log10(fEB)-MTp-0.02)/5+1)
        rows.append(dict(f=fEB,dchi2_lc=round(c2-cmin,2),chi2_ms=round(float(chi[j]),2),total=round(c2-cmin+float(chi[j]),2),M1=round(float(m1),3),M2=round(float(m2),3),
                         R1=round(float(R[i1][j]),3),k_fit=round(k,3),k_ms=round(float(kms[j]),3),J_fit=round(J,3),J_ms=round(float(Jms[j]),3),dist_pc=round(float(d),1)))
    out[name]=rows
    print(name)
    for r in rows: print(' ',r)
json.dump(out,open(B+'r09_msprior.json','w'),indent=1)
