# Referee: independent SuperWASP deep/shallow eclipse amplitude (own cleaning, own shape, linear baselines per night-camera)
import sys, json
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all')
import numpy as np
from astropy.io import fits
from rebm import pair_flux
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d=fits.open(B+'chromatic/data/wasp/1SWASP_J220333.07-425007.9.fits')[1].data
cam=np.array([int(s[:3]) for s in d['IMAGEID']])
col=sys.argv[1] if len(sys.argv)>1 else 'TAMFLUX2'
ok=(d['FLAG']==32)&np.isin(cam,[227,228])&np.isfinite(d[col])&(d[col]>0)
hjd=d['TMID'][ok]/86400.+2453005.5
t=hjd+65.184/86400.-2457000.0          # ~BJD_TDB in BTJD (HJD-BJD < 4 s ignored)
fl=d[col][ok].astype(float); cm=cam[ok]
night=np.floor(t-0.25).astype(int)     # local night (Chile/SAAO ~ UT-? ; SuperWASP-South at SAAO, UT+2)
segid=night*1000+cm
P=1.30367228; T0=3165.670416
# shape functions from referee model at lcmodel MS-like geometry, normalised to unit central depth
def shapes(tt,T0x):
    F=pair_flux(tt,T0x,P,0.112,0.441,0.0162,0.523,(0.3,0.3),(0.3,0.3))
    ph=((tt-T0x)/P)%1
    Fc=pair_flux(np.array([T0x,T0x+P/2]),T0x,P,0.112,0.441,0.0162,0.523,(0.3,0.3),(0.3,0.3))
    Sd=np.where((ph<0.25)|(ph>0.75),(1-F)/(1-Fc[0]),0.0)
    Ss=np.where((ph>0.25)&(ph<0.75),(1-F)/(1-Fc[1]),0.0)
    return Sd,Ss,1-Fc
# per-segment normalisation and clipping (iterative 4 sigma about per-segment line)
keep=np.ones(len(t),bool)
segs=np.unique(segid)
y=np.zeros(len(t)); wgt=np.zeros(len(t))
for it in range(4):
    for sg in segs:
        m=(segid==sg)&keep
        if m.sum()<10: keep[segid==sg]=False; continue
        mm=segid==sg
        med=np.median(fl[m]); yy=fl[mm]/med
        A=np.vstack([np.ones(mm.sum()),t[mm]-t[mm].mean()]).T
        c,*_=np.linalg.lstsq(A[keep[mm]],yy[keep[mm]],rcond=None)
        r=yy-A@c; sd=1.4826*np.median(np.abs(r[keep[mm]]))
        keep[mm]&=np.abs(r)<4*sd
        y[mm]=yy; wgt[mm]=1/sd**2
t,y,wgt,segid=t[keep],y[keep],wgt[keep],segid[keep]
segs,inv=np.unique(segid,return_inverse=True)
def amp(T0x):
    Sd,Ss,dep=shapes(t,T0x)
    # residualise y, Sd, Ss against per-segment linear baseline (weighted)
    out=[]
    for v in (y,Sd,Ss):
        rv=np.empty_like(v)
        for j in range(len(segs)):
            m=inv==j; tt=t[m]-t[m].mean(); W=wgt[m]
            A=np.vstack([np.ones(m.sum()),tt]).T
            c=np.linalg.solve(A.T@(A*W[:,None]),A.T@(v[m]*W)); rv[m]=v[m]-A@c
        out.append(rv)
    ry,rd,rs=out
    X=np.vstack([rd,rs]).T
    M=X.T@(X*wgt[:,None]); a=np.linalg.solve(M,X.T@(ry*wgt))
    chi=((ry-X@a)**2*wgt).sum()
    return -a, np.sqrt(np.diag(np.linalg.inv(M))), chi, ((Sd>0.5).sum(),(Ss>0.5).sum()), dep  # amplitude (depth) positive
a,ea,chi,nin,dep=amp(T0)
print(col,'N=',len(t),'segments',len(segs),'chi2/N',chi/len(t),'in-eclipse pts (deep,shallow)',nin)
print('deep depth = %.0f +- %.0f (formal) ppm; shallow = %.0f +- %.0f'%(a[0]*1e6,ea[0]*1e6,a[1]*1e6,ea[1]*1e6))
# fake-phase null: shift epoch by random phases avoiding the eclipses
rng=np.random.default_rng(1)
nul=[]
for ph in rng.uniform(0.08,0.42,150):
    sgn=rng.choice([-1,1])
    aa,_,_,_,_=amp(T0+sgn*ph*P)
    nul.append(aa)
nul=np.array(nul)
sd=nul.std(0); print('fake-phase null mean',nul.mean(0)*1e6,'sd',sd*1e6)
# timing scan
scan=[]
for dtm in np.arange(-20,21,2.0):
    aa,ee,cc,_,_=amp(T0+dtm/1440.); scan.append([dtm,aa[0]*1e6,cc])
scan=np.array(scan); print('timing scan: best dt (min)=',scan[np.argmin(scan[:,2]),0],'chi2 range',scan[:,2].min(),scan[:,2].max())
# ratios with dilution correction D=0.042 (chromatic) and TESS depth from lcmodel (5770+-100) / chromatic (5900+-172)
D=0.042
for TT,eT in ((5770,100),(5900,172)):
    r=a[0]*1e6*(1+D)/TT; er=r*np.sqrt((sd[0]/a[0])**2+(eT/TT)**2)
    print('SuperWASP/TESS deep ratio (TESS %d): %.3f +- %.3f (stat only)'%(TT,r,er))
json.dump(dict(col=col,N=int(len(t)),nseg=int(len(segs)),deep_ppm=a[0]*1e6,shallow_ppm=a[1]*1e6,formal=(ea*1e6).tolist(),
               null_sd=(sd*1e6).tolist(),null_mean=(nul.mean(0)*1e6).tolist(),ratio_deep_5770=a[0]*1e6*1.042/5770,
               scan=scan.tolist()),open(B+'referee_all/r05_wasp_%s.json'%col,'w'),indent=1)
