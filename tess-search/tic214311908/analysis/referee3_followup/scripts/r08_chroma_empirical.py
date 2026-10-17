"""r08: empirical chromatic predictions (independent of BT-Settl).
PTP (planet on an unresolved bound companion): companion/primary band contrast from empirical absolute-magnitude relations of
nearby single M dwarfs (data/mdw_seq.csv: d < 20 pc, RUWE < 1.2, Gaia XP synthetic g r i z V Ic + 2MASS Ks), as a function of
M_Ks (running median + smoothing spline). Primary at M_Ks = 6.441, companion at 6.441 + dKs, with dKs and weights from the
analyst's posterior-weighted draws (followup/out/ptp_draw_weights.npz, REF sign-fixed and INV sign-fixed); the relation is
differential, so the target's metallicity offset cancels to first order. Diluted depth in X proportional to
f_X = 1/(1 + 10^(0.4 dX)); band ratios X/Y = f_X/f_Y. Variant: shift the companion relation by the scatter of the relation.
WD occulted by the M3V: depth ratio g/TESS from observed M3V g (XP) and T (Stassun+19 relation from G, BP-RP) and a WD of
Teff 5500-7200 K given main-sequence-like colours from the Riello+21 G-g and Stassun G-T relations (blackbody check).
Output: out/r08_chroma.json"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd
from scipy.interpolate import UnivariateSpline
B='tic214311908/analysis'
R=f'{B}/referee3_followup'
d=pd.read_csv(f'{R}/data/mdw_seq.csv'); d=d[(d.MKs>5.3)&(d.MKs<9.3)]
dm=5*np.log10(d.parallax/100.)
bands={'g':'g_sdss_mag','r':'r_sdss_mag','i':'i_sdss_mag','z':'z_sdss_mag','V':'v_jkc_mag','Ic':'i_jkc_mag','G':'phot_g_mean_mag'}
rel={}; scat={}
xs=np.linspace(5.6,8.8,33)
for b,c in bands.items():
    y=(d[c]+dm).values; x=d.MKs.values; ok=np.isfinite(y)
    # running median in 0.3-mag windows of M_Ks, then smoothing spline
    med=np.array([np.median(y[ok&(np.abs(x-t)<0.25)]) for t in xs]); n=np.array([np.sum(ok&(np.abs(x-t)<0.25)) for t in xs])
    sp=UnivariateSpline(xs[n>=5], med[n>=5], k=3, s=len(xs)*0.003)
    rel[b]=sp; scat[b]=float(np.median(np.abs(y[ok]-sp(x[ok])))*1.4826)
print('N stars', len(d), 'scatter about relation', {k:round(v,3) for k,v in scat.items()})
MKp=6.441
print('relation at primary', {b:round(float(rel[b](MKp)),2) for b in bands}, ' target observed abs:', 
      {b:round(v-4.0405,2) for b,v in dict(g=16.363,r=14.983,i=13.575,z=12.815,V=15.591,Ic=12.859,G=14.281).items()})
W=np.load(f'{B}/followup/out/ptp_draw_weights.npz')
out={}
pct=lambda v,w,q=(0.025,0.16,0.5,0.84,0.975): [round(float(t),3) for t in np.interp(q, np.cumsum(w[np.argsort(v)])/w.sum(), np.sort(v))]
for case in ('ref_beta__signfixed','inv_beta__signfixed'):
    w=W[case].astype(float); dKs=W['dKs'].astype(float); keep=w>0
    w=w[keep]; dKs=dKs[keep]
    MKc=np.clip(MKp+dKs, 5.6, 8.8)
    dX={b: rel[b](MKc)-rel[b](MKp) for b in bands}
    f={b: 1/(1+10**(0.4*dX[b])) for b in bands}
    res={'dV':pct(dX['V'],w),'dIc':pct(dX['Ic'],w),'dz':pct(dX['z'],w)}
    for a,c in (('g','z'),('r','z'),('g','i'),('r','i'),('i','z'),('V','Ic')):
        res[f'{a}/{c}']=pct(f[a]/f[c],w)
    out[case]=res
    print(case, res)
# dKs relation check vs referee2 dI/dV
print('analyst draws dV/dI medians for comparison: see followup/out/ptp_chroma.json (BT-Settl g/z 0.785, r/z 0.75, V/Ic 0.80 REF)')
# WD
def T_from(G,x): return G-0.00522555*x**3+0.0891337*x**2-0.633923*x+0.0324473   # Stassun+19 (TIC v8) G->T
def g_from(G,x): return G-(0.2199-0.6365*x-0.1548*x**2+0.0064*x**3)              # Riello+21 G-g (SDSS), valid BP-RP<~2
gM, TM = 16.363, T_from(14.281,2.756)
wd={}
for Teff,x in ((5500,0.95),(6200,0.72),(7200,0.50)):
    gT_wd = g_from(0,x)-T_from(0,x)
    ratio_small = 10**(0.4*((gM-TM)-gT_wd))
    a=0.0152/(1-0.0152)   # F_WD/F_M in TESS
    dg=a*ratio_small/(1+a*ratio_small)
    wd[Teff]=dict(g_minus_T_WD=round(gT_wd,3), g_over_TESS=round(dg/0.0152,2))
print('M3V g-T', round(gM-TM,3), 'T', round(TM,3), 'WD', wd)
out['WD']=dict(M3V_T=TM, M3V_g_minus_T=gM-TM, by_Teff=wd, note='analyst (BT-Settl M3V, 0.43 mag too bright in g): 6.0 (4.9-7.4)')
out['relation_scatter']=scat
json.dump(out, open(f'{R}/out/r08_chroma.json','w'), indent=1)
