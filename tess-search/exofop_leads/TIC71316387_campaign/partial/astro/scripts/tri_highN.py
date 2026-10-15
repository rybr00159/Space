"""High-N TRICERATOPS evidence for the bound-companion EB (SEB, SEBx2P) vs background EB (BEB, BEBx2P) scenarios.
Same light-curve handling as the FIXED runner (/home/user/work/src/run_tri.py): fold, +-3 dur window, 5-min bins, depth from 3-bin
'valid' smoothing, calc_depths(all_ap_pixels=None). The target object is built once; then calc_probs is called repeatedly with
parallel=True (TRICERATOPS's own vectorised likelihoods) in chunks of N draws, dropping scenarios excluded by independent data:
 - TP, PTP, STP, DTP, BTP (planets): excluded by the unequal odd/even eclipses (5,900 vs 2,700 ppm) at P_CTOI;
 - EB/EBx2P, PEB/PEBx2P, DEB/DEBx2P (G star itself eclipsing): excluded by the constant Gaia RV (would need K~100 km/s);
 - resolved neighbours (NEB/NTP): excluded by TESS difference-image localisation (remove_nearby_stars()).
Chunks have equal N, so the pooled evidence is exactly log(mean(exp(lnZ_chunk))).
usage: tri_highN.py sectors P lcfile out.json nchunks Nchunk seed [dropextra]"""
import sys, json, os, warnings, time
import numpy as np
warnings.filterwarnings("ignore")
import triceratops.triceratops as tr
sectors=[int(s) for s in sys.argv[1].split(',')]; P=float(sys.argv[2]); lcfile=sys.argv[3]; out=sys.argv[4]
nch=int(sys.argv[5]); Nch=int(sys.argv[6]); seed=int(sys.argv[7]); dur=0.064; t0=0.0
np.random.seed(seed)
z=np.load(lcfile)
ph=((z["t"]-t0)/P+0.5)%1-0.5; m=np.abs(ph*P)<3*dur
time_d,flux,ferr=ph[m]*P,z["f"][m],np.full(m.sum(),float(z["sigma"]))
o=np.argsort(time_d); time_d,flux,ferr=time_d[o],flux[o],ferr[o]
edges=np.arange(time_d.min(),time_d.max()+5/1440,5/1440); idx=np.digitize(time_d,edges)
tb,fb,eb=[],[],[]
for k in np.unique(idx):
    s=idx==k
    if s.sum()>=3: tb.append(time_d[s].mean()); fb.append(flux[s].mean()); eb.append(ferr[s].mean()/np.sqrt(s.sum()))
tb,fb,eb=map(np.array,(tb,fb,eb))
sig_scale=float(os.environ.get('TRI_SIGSCALE',1.0))
target=tr.target(ID=71316387,sectors=sectors,search_radius=10)
target.stars.loc[0,"mass"]=float(os.environ.get('TRI_MASS',1.05)); target.stars.loc[0,"rad"]=float(os.environ.get('TRI_RAD',1.07)); target.stars.loc[0,"plx"]=float(os.environ.get('TRI_PLX',3.3622))
tdepth=float(1-np.min(np.convolve(fb,np.ones(3)/3,"valid")))
target.calc_depths(tdepth=tdepth,all_ap_pixels=None)
target.remove_nearby_stars()
drop=["TP","PTP","STP","DTP","BTP","EB","PEB","DEB"]+([s for s in sys.argv[8].split(',')] if len(sys.argv)>8 else [])
kw={}
if os.environ.get("TRI_CC"): kw=dict(contrast_curve_file=os.environ["TRI_CC"],filt=os.environ.get("TRI_CC_FILT","TESS"))
res=dict(tdepth=tdepth,flux_err_0=float(np.median(eb))*sig_scale,Nchunk=Nch,P=P,sectors=sectors,lcfile=lcfile,seed=seed,drop=drop,cc=kw.get('contrast_curve_file'),chunks=[])
for c in range(nch):
    t1=time.time()
    target.calc_probs(time=tb,flux_0=fb,flux_err_0=float(np.median(eb))*sig_scale,P_orb=P,parallel=True,N=Nch,verbose=0,drop_scenario=drop,**kw)
    df=target.probs; lz={s:float(l) for s,l in zip(df.scenario,target.lnZ)}
    bestp={s:dict(M_s=float(r.M_s),R_s=float(r.R_s),inc=float(r.inc),b=float(r.b),M_EB=float(r.M_EB),R_EB=float(r.R_EB),ecc=float(r.ecc)) for s,r in zip(df.scenario,df.itertuples()) if s in ('SEB','SEBx2P','BEB','BEBx2P')}
    res['chunks'].append(dict(lnZ={k:lz[k] for k in ('SEB','SEBx2P','BEB','BEBx2P')},best=bestp))
    json.dump(res,open(out,'w'),indent=1)
    print(c,round(time.time()-t1,1),{k:round(lz[k],2) for k in ('SEB','SEBx2P','BEB','BEBx2P')},flush=True)
