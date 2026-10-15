# Referee: independent regression of TESS lateral-colour effect from localise f01 PRF-fit tables
# (uses their measured PRF positions and Gaia/WCS positions, but my own detrending and regression)
import numpy as np, glob, json
from scipy.spatial import cKDTree
from scipy.stats import theilslopes
R='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3/results/'
out={}
allC=[];allR=[]
for f in sorted(glob.glob(R+'f01_tess*.npz')):
    z=np.load(f,allow_pickle=True); a=z['stars']; J=z['J']
    d=J@np.vstack([a['xf']-a['xw'],a['yf']-a['yw']])   # arcsec (E,N)
    pa=np.radians(a['pa_bs'])
    u=np.vstack([np.sin(pa),np.cos(pa)])              # PA east of north -> (E,N) unit vector toward boresight
    ok=(a['ruwe']<1.4)&np.isfinite(a['C'])&(a['chi2r']<np.nanpercentile(a['chi2r'],95))&(a['G']<14.5)
    a=a[ok]; d=d[:,ok]; u=u[:,ok]
    # nonparametric detrend: subtract median residual of 25 nearest neighbours (position) with colour within 0.3 of 0.8
    xy=np.vstack([a['xw'],a['yw']]).T; tr=cKDTree(xy)
    ref=np.abs(a['C']-0.8)<0.35
    trr=cKDTree(xy[ref]); dref=d[:,ref]
    _,idx=trr.query(xy,k=26)
    dd=d.copy()
    for j in range(2):
        nb=dref[j][idx[:,1:]]; dd[j]=d[j]-np.median(nb,axis=1)
    rad=dd[0]*u[0]+dd[1]*u[1]   # + = toward boresight
    tan=-dd[0]*u[1]+dd[1]*u[0]
    sc=a['theta']/float(z['theta_t'])
    x=(a['C']-0.786)*sc
    m=np.abs(rad)<3
    s,ic,lo,hi=theilslopes(rad[m],x[m])
    st,_,_,_=theilslopes(tan[m],x[m])
    # also OLS with magnitude term
    A=np.vstack([x[m],a['G'][m]-12,np.ones(m.sum())]).T
    p,*_=np.linalg.lstsq(A,rad[m],rcond=None)
    sec=int(z['sec'])
    out[sec]=dict(n=int(m.sum()),k_rad_theilsen=s,k_rad_ci=[lo,hi],k_tan=st,k_rad_ols_withG=p[0],G_term=p[1])
    allC.append(a['C'][m]); allR.append(rad[m]/sc[m])
    print(sec,out[sec])
C=np.concatenate(allC);Rr=np.concatenate(allR)
bins=[(0.55,0.7),(0.7,0.85),(0.85,1.0),(1.0,1.3),(1.3,1.7),(1.7,2.1),(2.1,2.4),(2.4,2.8),(2.8,3.5)]
curve=[]
for lo,hi in bins:
    mm=(C>=lo)&(C<hi)
    if mm.sum()<3: continue
    med=np.median(Rr[mm]); err=1.2533*np.std(Rr[mm])/np.sqrt(mm.sum())
    curve.append([lo,hi,int(mm.sum()),float(np.median(C[mm])),float(med),float(err)])
    print(f'C {lo}-{hi} n={mm.sum()} medC={np.median(C[mm]):.2f} rad(toward boresight)={med:+.2f}+-{err:.2f}')
out['curve']=curve
json.dump(out,open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all/r01_colour.json','w'),indent=1,default=float)
