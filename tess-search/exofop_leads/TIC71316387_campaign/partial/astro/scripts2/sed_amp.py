"""SED-only upper limits on the light of a main-sequence pair at the G star's distance: chi2 profile vs amplitude A
(A=1: pair on the Mamajek sequence at 297 pc), Teff/R/E(B-V) of the G star free. Uses sed_fit.py machinery unchanged.
Also: same with a Gaussian prior on the G-star Teff (GSP-Spec ANN 5843+-150 K or MG 6235+-160 K)."""
import sys, json, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/scripts')
from sed_fit import *
from scipy.optimize import least_squares
plx=3.3622; d=1000/plx
bands=[b for b,v in PH.items() if v[4]]
m=np.array([PH[b][0] for b in bands]); e=np.array([np.hypot(PH[b][1],PH[b][3]) for b in bands])
A_ext=np.array([TAB[b][4] for b in bands])
bbT=['TESS__TESS.Red']
def fitA(pair,A,tprior=None):
    Mp=pair_abs(pair,bands) if A>0 else None
    def res(p):
        mg=gmodel(bands,p[0],p[1],d,p[2])
        mod=mg if A==0 else S.addmag(mg,Mp+5*np.log10(d/10)+A_ext*p[2]-2.5*np.log10(A))
        r=[(m-mod)/e,[(p[2]-0.009)/0.01]]
        if tprior: r.append([(p[0]-tprior[0])/tprior[1]])
        return np.concatenate(r)
    best=None
    for t0 in (5600,5900,6300):
        rr=least_squares(res,[t0,1.1,0.01],bounds=([4800,0.5,0],[6800,2.5,0.3]),x_scale=[100,0.05,0.01])
        if best is None or rr.cost<best.cost: best=rr
    p=best.x
    if A>0:
        mg=gmodel(bbT,p[0],p[1],d,p[2])[0]; mp=pair_abs(pair,bbT)[0]+5*np.log10(d/10)+TAB[bbT[0]][4]*p[2]-2.5*np.log10(A)
        fT=1/(1+10**(-0.4*(mg-mp)))
    else: fT=0.0
    return 2*best.cost,fT,p[0]
grid=np.array([0,0.25,0.5,1,1.5,2,3,4,6,8,12,16])
out={}
for pr in [('K7V','M2V'),('M0V','M3V'),('M0V','M4.5V'),('M1V','M5V'),('M1.5V','M5V'),('M2V','M5V'),('M2.5V','M5V')]:
    for lab,tp in [('free',None),('ANN5843',(5843,150)),('MG6235',(6235,160))]:
        rows=[fitA(pr,A,tp) for A in grid]; c=np.array([r[0] for r in rows]); f=np.array([r[1] for r in rows])
        cmin=c.min(); dc=c-cmin
        # 95% upper limit on A (dchi2=3.84 above minimum, beyond minimum)
        imin=int(np.argmin(c)); up=None
        for i in range(imin,len(grid)-1):
            if dc[i]<3.84<=dc[i+1]:
                up=float(np.interp(3.84,[dc[i],dc[i+1]],[grid[i],grid[i+1]])); fup=float(np.interp(3.84,[dc[i],dc[i+1]],[f[i],f[i+1]])); break
        out['|'.join(pr)+'_'+lab]=dict(A=grid.tolist(),chi2=c.round(2).tolist(),fTESS=f.round(4).tolist(),TeffG=[round(r[2]) for r in rows],
              A_best=float(grid[imin]),A95=up,fTESS95=(fup if up else None),dchi2_A1_vs_A0=float(c[3]-c[0]))
        print(pr,lab,'A_best',grid[imin],'A95',up,'fTESS95',round(fup,3) if up else None,'dchi2(A=1)-(A=0)=%.2f'%(c[3]-c[0]),'chi2(A=0)=%.1f'%c[0],flush=True)
json.dump(out,open(D+'results/sed_amp_profiles.json','w'),indent=1)
