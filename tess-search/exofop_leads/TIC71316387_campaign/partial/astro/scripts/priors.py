"""Prior odds: bound hierarchical triple vs chance-aligned background EB, for the hidden EB near TIC 71316387."""
import numpy as np, json
from astropy.table import Table
from scipy.stats import norm
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
rng=np.random.default_rng(3)
T=Table.read(D+'data/gaia_counts_30arcmin.ecsv')
sep=np.array(T['sep_arcmin']); G=np.array(T['phot_g_mean_mag'],float); br=np.array(T['bp_rp'].filled(np.nan) if hasattr(T['bp_rp'],'filled') else T['bp_rp'],float)
keep=(sep>1.0); area_deg2=np.pi*(0.5**2-(1/60)**2)
def tmag(G,br):
    x=np.where(np.isfinite(br),np.clip(br,-0.2,3.5),1.0)
    return G-0.00522555*x**3+0.0891337*x**2-0.633923*x+0.0324473
Tm=tmag(G,br)[keep]; brk=br[keep]; Gk=G[keep]
T_tot=11.358; d_obs=0.0059
dT=Tm-T_tot; f=1/(1+10**(0.4*dT)); dreq=d_obs/f
# completeness check: differential counts per 0.5 mag
hist=np.histogram(Gk,bins=np.arange(10,20.5,0.5))[0]
# EB statistics
K=Table.read(D+'data/kepler_eb_kirk2016.ecsv'); P=np.array(K['Per'],float); mo=np.array(K['Morph'],float)
NT_kep=2878/0.013   # Kirk et al. 2016: EBs are 1.3% of Kepler targets
win=(P>=1.0)&(P<=2.0)&(mo<0.7)
r_EB=win.sum()/NT_kep
r_EB_all12=((P>=1.0)&(P<=2.0)).sum()/NT_kep
E=Table.read(D+'data/tess_ebs_prsa2022.ecsv')
Pe=np.array(E['Per'],float); me=np.array(E['Morph'],float)
dp=np.array(E['Dp-pf'].filled(np.nan) if hasattr(E['Dp-pf'],'filled') else E['Dp-pf'],float)
dp2=np.array(E['Dp-2g'].filled(np.nan) if hasattr(E['Dp-2g'],'filled') else E['Dp-2g'],float)
dp=np.where(np.isfinite(dp),dp,dp2)
we=(Pe>=1.0)&(Pe<=2.0)&(me<0.7)&np.isfinite(dp)&(dp>0)
dps=np.sort(dp[we])
def qdepth(x): return 1-np.searchsorted(dps,x)/len(dps)
out=dict(field=dict(n_stars=int(keep.sum()),area_deg2=area_deg2,dN_per_halfmag_G10_20=hist.tolist()),
         EB=dict(kepler_n_window=int(win.sum()),kepler_rate_per_star_P1to2_morph_lt0p7=r_EB,kepler_rate_per_star_P1to2_all=r_EB_all12,
                 tess_eb_n_window=int(we.sum()),tess_depth_quantiles={str(qq):float(np.quantile(dps,qq)) for qq in [0.1,0.25,0.5,0.75,0.9]}))
print(json.dumps(out,indent=1))
res={}
for dmax in [0.5,1.0]:
    for colour in ['any','red_bprp_ge_1.5']:
        ok=(dreq<=dmax)
        if colour!='any': ok&=(brk>=1.5)
        w=qdepth(dreq[ok])
        n_dens=ok.sum()/area_deg2/3600**2        # stars per arcsec^2 able to host the signal
        wdens=np.sum(w)/area_deg2/3600**2        # depth-weighted
        for r in [1.0,1.24,1.8,4.0]:
            A=np.pi*r**2
            res[f'dmax{dmax}_{colour}_r{r}']=dict(N_hosts=float(n_dens*A),P_bgEB_nodepth=float(n_dens*A*r_EB),P_bgEB=float(wdens*A*r_EB),
                     n_stars_ok=int(ok.sum()),Tmax=float(Tm[ok].max()) if ok.sum() else None)
for k,v in res.items(): print(k,{kk:(f'{vv:.3g}' if isinstance(vv,float) else vv) for kk,vv in v.items()})
# ---- bound prior ----
mu,sig=5.03,2.28; csf=0.33+2*0.08+3*0.03   # Raghavan+2010
dpc=1000/3.3622; Mtot=1.03+0.8
def fsep(s1,s2,fac=1.26):
    a1,a2=fac*s1*dpc,fac*s2*dpc
    P1=np.log10(np.sqrt(a1**3/Mtot)*365.25); P2=np.log10(np.sqrt(a2**3/Mtot)*365.25)
    return norm.cdf((P2-mu)/sig)-norm.cdf((P1-mu)/sig),P1,P2
bound={}
for (s1,s2) in [(0.05,1.24),(0.1,1.24),(0.1,1.8),(0.2,1.0)]:
    fs,P1,P2=fsep(s1,s2)
    for qlab,qf in [('q0.5-1.0',0.5/0.9),('q0.5-0.73',0.23/0.9)]:
        for Elab,Ef in [('E1',1.0),('E3',3.0)]:
            # bound pair supplies 2.5-9% of TESS light -> required undiluted depth 6.5-24%
            qd=[qdepth(d_obs/ff) for ff in (0.025,0.05,0.09)]
            Pb=csf*fs*qf*r_EB*Ef*np.array(qd)
            bound[f'sep{s1}-{s2}_{qlab}_{Elab}']=dict(f_companion_in_window=float(csf*fs*qf),logP_range=(round(P1,2),round(P2,2)),
                    depth_qual=[round(float(x),3) for x in qd],P_bound=[float(x) for x in Pb])
for k,v in bound.items(): print(k,v)
# ---- odds for the localisation-limited case (r<=1.24") ----
odds={}
bg_lo=res['dmax0.5_red_bprp_ge_1.5_r1.24']['P_bgEB']; bg_hi=res['dmax1.0_any_r1.24']['P_bgEB_nodepth']
b_lo=min(min(v['P_bound']) for k,v in bound.items() if k.startswith('sep0.1-1.24'))
b_hi=max(max(v['P_bound']) for k,v in bound.items() if k.startswith('sep0.1-1.24'))
odds['r1.24']=dict(P_bound_range=(b_lo,b_hi),P_bg_range=(bg_lo,bg_hi),odds_min=b_lo/bg_hi,odds_max=b_hi/bg_lo,
                   odds_central=float(bound['sep0.1-1.24_q0.5-1.0_E1']['P_bound'][1]/res['dmax0.5_any_r1.24']['P_bgEB']))
for r in [1.0,1.8,4.0]:
    odds[f'r{r}_bg_only']={'P_bg_any_dmax0.5':res[f'dmax0.5_any_r{r}']['P_bgEB'],'P_bg_any_dmax1_nodepth':res[f'dmax1.0_any_r{r}']['P_bgEB_nodepth'],'P_bg_red_dmax0.5':res[f'dmax0.5_red_bprp_ge_1.5_r{r}']['P_bgEB']}
print(json.dumps(odds,indent=1))
json.dump(dict(inputs=out,background=res,bound=bound,odds=odds),open(D+'results/priors.json','w'),indent=1)
