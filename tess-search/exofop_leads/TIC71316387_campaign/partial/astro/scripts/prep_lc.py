"""Prepare locally-detrended eclipse light curves (PDCSAP and QLP SAP) for TRICERATOPS and depth checks.
Local baseline: per eclipse event, a quadratic fit to out-of-eclipse data within +-0.22 d of the event,
excluding |dt|<0.055 d (half-duration 0.032 d + margin)."""
import numpy as np, json, sys
from astropy.io import fits
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/lc/'
OUT='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/tri/'
P=1.3036712; T0=3165.6697
WIN=0.22; MASK=0.055
def load(fn, kind):
    h=fits.open(D+fn); d=h[1].data
    t=np.array(d['TIME'],float)
    if kind=='pdc':
        f=np.array(d['PDCSAP_FLUX'],float); q=np.array(d['QUALITY'])
        good=np.isfinite(t)&np.isfinite(f)&(q==0)
    else:
        f=np.array(d['SAP_FLUX'],float); q=np.array(d['QUALITY'])
        good=np.isfinite(t)&np.isfinite(f)&(q==0)
    t,f=t[good],f[good]; f=f/np.nanmedian(f)
    return t,f
def events(t,f,phase0):
    """phase0=0 deep, 0.5 shallow; returns dt, normalised flux, event id"""
    n=np.round((t-T0)/P-phase0); tc=T0+(n+phase0)*P
    dt=t-tc; out=[]
    for k in np.unique(n[np.abs(dt)<WIN]):
        s=(n==k)&(np.abs(dt)<WIN)
        oot=s&(np.abs(dt)>MASK); ins=s&(np.abs(dt)<0.02)
        if oot.sum()<10 or ins.sum()<2: continue
        # require coverage both sides
        if (dt[oot]<0).sum()<4 or (dt[oot]>0).sum()<4: continue
        c=np.polyfit(dt[oot],f[oot],2); base=np.polyval(c,dt[s])
        out.append((dt[s],f[s]/base,np.full(s.sum(),k)))
    if not out: return None
    return [np.concatenate(x) for x in zip(*out)]
def depth_core(dt,fn,half=0.008):
    """mean depth in central +-half d (approx 23 min) - simple, robust"""
    m=np.abs(dt)<half; return 1-np.mean(fn[m])
res={}; store={}
for fn,kind,sec,cad in [('spoc_s0001.fits','pdc',1,1800),('spoc_s0028.fits','pdc',28,600),('spoc_s0068.fits','pdc',68,200),
                        ('qlp_s0095.fits','qlp',95,200),('qlp_s0105.fits','qlp',105,200)]:
    t,f=load(fn,kind)
    r={}
    for nm,ph in [('deep',0.0),('shal',0.5)]:
        e=events(t,f,ph)
        if e is None: continue
        dt,fn_,ev=e
        oot=np.abs(dt)>MASK; sig=1.4826*np.median(np.abs(fn_[oot]-np.median(fn_[oot])))
        d=depth_core(dt,fn_)
        # empirical null: same estimator at fake phases (window centred at phase 0.25/0.75 +- jitter)
        nulls=[]
        rng=np.random.default_rng(1)
        for ph0 in rng.uniform(0.12,0.38,200).tolist()+rng.uniform(0.62,0.88,200).tolist():
            e2=events(t,f,ph0)
            if e2 is None: continue
            nulls.append(depth_core(e2[0],e2[1]))
        nulls=np.array(nulls)
        r[nm]=dict(n_events=int(len(np.unique(ev))),depth_core_ppm=round(d*1e6,1),null_std_ppm=round(np.std(nulls)*1e6,1),
                   null_mean_ppm=round(np.mean(nulls)*1e6,1),sigma_point=float(sig))
        store[(sec,nm)]=(dt,fn_,sig,cad)
    res[f'S{sec}_{kind}']=r; print(sec,kind,r,flush=True)
json.dump(res,open(OUT+'depths_core.json','w'),indent=1)
# ---- build TRICERATOPS inputs ----
def build(keys, name, ref_sigma_cad=200):
    """fold deep+shallow together (P_CTOI ~ P_EB/2). Longer cadences are replicated (cad/200 copies) so that the runner's
    equal-weight binning is equivalent to inverse-variance weighting; per-point sigma of 200-s data is used."""
    T=[];F=[];sig200=None
    for k in keys:
        dt,fn_,sig,cad=store[k]
        rep=int(round(cad/ref_sigma_cad))
        if cad==ref_sigma_cad: sig200=sig if sig200 is None else np.sqrt(0.5*(sig200**2+sig**2))
        T.append(np.repeat(dt,rep)); F.append(np.repeat(fn_,rep))
    if sig200 is None:  # derive equivalent from longer cadence
        dt,fn_,sig,cad=store[keys[0]]; sig200=sig*np.sqrt(cad/ref_sigma_cad)
    T=np.concatenate(T);F=np.concatenate(F)
    np.savez(OUT+name+'.npz',t=T,f=F,sigma=sig200)
    print(name,len(T),'sigma',sig200)
build([(68,'deep'),(68,'shal')],'pdc_s68_fold')
build([(28,'deep'),(28,'shal'),(68,'deep'),(68,'shal')],'pdc_s28s68_fold')
build([(95,'deep'),(95,'shal'),(105,'deep'),(105,'shal')],'qlp_s95s105_fold')
build([(68,'deep')],'pdc_s68_deeponly')
build([(28,'deep'),(68,'deep')],'pdc_s28s68_deeponly')
build([(68,'shal')],'pdc_s68_shalonly')
