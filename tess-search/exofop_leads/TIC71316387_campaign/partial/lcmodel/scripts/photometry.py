"""Own aperture photometry on the 11x11 TESScut sub-cubes (all six sectors), with
 * target-specific barycentric time correction (baryfix.py; validated vs SPOC TIMECORR to <1 ms),
 * per-cadence sky from the faintest-25% pixels of the full 21x21 stamp,
 * several apertures; the one with the lowest 1-h OOE scatter is kept,
 * PRF-based contamination fraction from Gaia DR3 neighbours (TESS_PRF model) for the chosen aperture.
Output: data/lc_tc_s{sec}.npz with t (BJD_TDB-2457000), f (normalised), aperture mask, contamination."""
import numpy as np, json, glob, os, sys
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
import pandas as pd
HERE='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel'
P=1.3036712; T0=3165.6697; T14=1.55/24
def tmag_from_gaia(g, c):
    return g - 0.00522555*c**3 + 0.0891337*c**2 - 0.633923*c + 0.0324473
gaia=pd.read_csv(f'{HERE}/data/gaia_240as.csv')
gaia['bp_rp']=gaia['bp_rp'].fillna(1.0)
gaia['tmag']=tmag_from_gaia(gaia.phot_g_mean_mag,gaia.bp_rp)
def ecl_mask(t, w=0.75):
    ph=((t-T0)/P)%1.0
    d=np.minimum(np.abs(ph*P),np.abs((1-ph)*P)); d2=np.abs((ph-0.5)*P)
    return (d<w*T14)|(d2<w*T14)
def cdpp1h(t,f,m):
    # scatter of 1-h bins after removing a 0.5-d running median, OOE only
    tt,ff=t[~m],f[~m]
    out=np.empty_like(ff)
    for i in range(len(tt)):
        s=np.abs(tt-tt[i])<0.25; out[i]=np.median(ff[s])
    r=ff-out; nb=np.floor((tt-tt[0])/(1/24)).astype(int)
    b=np.array([r[nb==k].mean() for k in np.unique(nb) if (nb==k).sum()>=0.6*np.median(np.bincount(nb))])
    return 1.4826*np.median(np.abs(b-np.median(b)))
res={}
import PRF
for fn in sorted(glob.glob(f'{HERE}/data/cube_s*.npz'), key=lambda s:int(s.split('_s')[-1][:-4])):
    z=np.load(fn,allow_pickle=True); sec=int(fn.split('_s')[-1][:-4])
    t=z['t']; q=z['q']; F=z['F'].astype(float); bk=z['bk']
    bc=np.load(f'{HERE}/data/bary_tc_s{sec}.npy')
    ok=np.isfinite(t)&(q==0)&np.all(np.isfinite(F.reshape(len(t),-1)),axis=1)&np.isfinite(bk)
    # time correction: match by TIME
    corr=np.interp(t,bc[0],bc[1]); tt=t+corr/86400
    t1,F1,b1=tt[ok],F[ok],bk[ok]
    img=F1-b1[:,None,None]
    med=np.median(img,axis=0)
    w=WCS(fits.Header.fromstring(str(z['wcs_sub'])))
    x0,y0=[float(v) for v in w.all_world2pix(330.887948,-42.835511,0)]
    m=ecl_mask(t1)
    best=None; trials=[]
    yy,xx=np.mgrid[:11,:11]
    r=np.hypot(xx-x0,yy-y0)
    for kind,par in [('r',1.0),('r',1.3),('r',1.6),('r',2.0),('r',2.5),('thr',0.15),('thr',0.08),('thr',0.04)]:
        if kind=='r': ap=r<=par
        else:
            ap=(med>par*med.max())&(r<3.5)
        if ap.sum()<2: continue
        fl=img[:,ap].sum(1); f=fl/np.median(fl[~m])
        c=cdpp1h(t1[::1],f,m)
        trials.append((kind,par,int(ap.sum()),c))
        if best is None or c<best[3]: best=(kind,par,ap.copy(),c,f,fl)
    kind,par,ap,c,f,fl=best
    print(f'S{sec}: target at ({x0:.2f},{y0:.2f}); trials', [(a,b,n,round(cc*1e6)) for a,b,n,cc in trials], '-> best',kind,par,int(ap.sum()),f'{c*1e6:.0f} ppm/1h', flush=True)
    # PRF contamination: model each Gaia star in the 11x11 stamp with the TESS PRF at its position
    cam,ccd=int(z['cam']),int(z['ccd']); col0,row0=float(z['col0']),float(z['row0'])
    ix,iy=int(z['ix']),int(z['iy'])
    colc=col0+ix; rowc=row0+iy   # approximate CCD column/row of target
    try:
        prf=PRF.TESS_PRF(cam,ccd,sec,colc,rowc,localdatadir=f'{HERE}/prf_files/'+('s0001' if sec<4 else 's0004'))
    except Exception as e:
        print('PRF error',e); prf=None
    ep=2016.0
    contam=None
    if prf is not None:
        dt_yr=(np.median(t1)+2457000-2451545)/365.25+2000-ep
        ra=gaia.ra+gaia.pmra.fillna(0)/3.6e6/np.cos(np.deg2rad(gaia.dec))*dt_yr
        de=gaia.dec+gaia.pmdec.fillna(0)/3.6e6*dt_yr
        px,py=w.all_world2pix(ra.values,de.values,0)
        tot=np.zeros((11,11)); tgt=None; flux_in={}
        for i in range(len(gaia)):
            if px[i]<-8 or py[i]<-8 or px[i]>18 or py[i]>18: continue
            # PRF.locate wants positions relative to the stamp; clip far-away stars
            try:
                model=prf.locate(px[i],py[i],(11,11))
            except Exception:
                continue
            fl_=10**(-0.4*(gaia.tmag.values[i]-10))*model
            if i==0: tgt=fl_
            tot+=fl_
        contam=float(1-tgt[ap].sum()/tot[ap].sum())
        flux_frac=float(tgt[ap].sum()/tgt.sum())
        # check of model vs data shape: fraction of median image in aperture
        print(f'   PRF contamination in aperture {contam*100:.2f}%  (target flux fraction captured {flux_frac:.3f})')
    np.savez(f'{HERE}/data/lc_tc_s{sec}.npz',t=t1,f=f,flux=fl,bk=b1,ap=ap,contam=contam if contam is not None else np.nan,x0=x0,y0=y0,
             cad=np.median(np.diff(t1)), centroid=np.vstack([(img*xx).reshape(len(t1),-1)[:,ap.ravel()].sum(1)/fl,(img*yy).reshape(len(t1),-1)[:,ap.ravel()].sum(1)/fl]))
    res[sec]=dict(aperture=[kind,par,int(ap.sum())],cdpp1h_ppm=c*1e6,contam=contam,trials=trials,n=len(t1),bk_med=float(np.median(b1)))
json.dump(res,open(f'{HERE}/results/photometry.json','w'),indent=1,default=str)
