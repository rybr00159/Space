"""Independent arithmetic check of results/priors.json, plus same-window (IPD-informed annulus) and colour-matched variants.
Re-uses only the stored inputs (Kepler EB rate, TESS-EB depth quantile function via the stored file, Gaia counts)."""
import numpy as np, json
from scipy.stats import norm
from astropy.table import Table
S='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
P=json.load(open(S+'results/priors.json'))
r_EB=P['inputs']['EB']['kepler_rate_per_star_P1to2_morph_lt0p7']
print('r_EB',r_EB,'= 225/(2878/0.013) ->',225/(2878/0.013))
# bound: recompute f_sep
mu,sig=5.03,2.28; csf=0.33+2*0.08+3*0.03; dpc=1000/3.3622; Mtot=1.03+0.8
def fsep(s1,s2,fac=1.26,Mt=Mtot):
    lp=lambda s: np.log10(np.sqrt((fac*s*dpc)**3/Mt)*365.25)
    return norm.cdf((lp(s2)-mu)/sig)-norm.cdf((lp(s1)-mu)/sig),lp(s1),lp(s2)
fs,p1,p2=fsep(0.1,1.24); print('csf',csf,'fsep(0.1-1.24)',round(fs,4),'logP',round(p1,2),round(p2,2),'csf*fs*qf',csf*fs*0.5/0.9,'stored',P['bound']['sep0.1-1.24_q0.5-1.0_E1']['f_companion_in_window'])
# TESS-EB depth tail function (rebuild exactly as priors.py)
E=Table.read(S+'data/tess_ebs_prsa2022.ecsv')
Pe=np.array(E['Per'],float); me=np.array(E['Morph'],float)
dp=np.array(E['Dp-pf'].filled(np.nan) if hasattr(E['Dp-pf'],'filled') else E['Dp-pf'],float)
dp2=np.array(E['Dp-2g'].filled(np.nan) if hasattr(E['Dp-2g'],'filled') else E['Dp-2g'],float)
dp=np.where(np.isfinite(dp),dp,dp2); we=(Pe>=1)&(Pe<=2)&(me<0.7)&np.isfinite(dp)&(dp>0); dps=np.sort(dp[we])
q=lambda x: 1-np.searchsorted(dps,x)/len(dps)
qd=[q(0.0059/f) for f in (0.025,0.05,0.09)]; print('bound depth-quality',np.round(qd,3),'stored',P['bound']['sep0.1-1.24_q0.5-1.0_E1']['depth_qual'])
# background: rebuild host density
T=Table.read(S+'data/gaia_counts_30arcmin.ecsv')
sep=np.array(T['sep_arcmin']); G=np.array(T['phot_g_mean_mag'],float); br=np.array(T['bp_rp'].filled(np.nan) if hasattr(T['bp_rp'],'filled') else T['bp_rp'],float)
keep=sep>1.0; area=np.pi*(0.5**2-(1/60)**2)*3600**2
x=np.where(np.isfinite(br),np.clip(br,-0.2,3.5),1.0); Tm=(G-0.00522555*x**3+0.0891337*x**2-0.633923*x+0.0324473)[keep]; brk=br[keep]
f=1/(1+10**(0.4*(Tm-11.358))); dreq=0.0059/f
def bg(dmax,red,A,depthw=True):
    ok=dreq<=dmax
    if red: ok&=brk>=1.5
    w=q(dreq[ok]) if depthw else np.ones(ok.sum())
    return w.sum()/area*A*r_EB, ok.sum()/area*A
A124=np.pi*1.24**2
b_any,nh=bg(0.5,False,A124); print('bg any dmax0.5 r1.24: P=%.3g (stored %.3g) N_hosts=%.3g'%(b_any,P['background']['dmax0.5_any_r1.24']['P_bgEB'],nh))
out={}
# (a) stored central and range
Pb_c=P['bound']['sep0.1-1.24_q0.5-1.0_E1']['P_bound'][1]
out['central_check']=Pb_c/b_any
# (b) colour-matched: bound q 0.5-0.73 (M-dwarf pair 0.52-0.75 Msun) vs background BP-RP>=1.5 hosts
Pb_red=P['bound']['sep0.1-1.24_q0.5-0.73_E1']['P_bound'][1]
b_red,_=bg(0.5,True,A124); b_red1,_=bg(1.0,True,A124)
out['colour_matched_central']=Pb_red/b_red; out['colour_matched_dmax1']=Pb_red/b_red1
# (c) same-window annulus 0.48-1.2 arcsec (Gaia IPD-informed) for both
fsA,_,_=fsep(0.48,1.2); Aann=np.pi*(1.2**2-0.48**2)
for lab,qf,red in [('any',0.5/0.9,False),('colour',0.23/0.9,True)]:
    Pb=[csf*fsA*qf*r_EB*Ef*qd[1] for Ef in (1,3)]
    bgA,_=bg(0.5,red,Aann); bgA1,_=bg(1.0,red,Aann,depthw=False)
    out[f'annulus0.48-1.2_{lab}']=dict(fsep=fsA,P_bound_E1=Pb[0],P_bound_E3=Pb[1],P_bg_dmax0p5=bgA,P_bg_dmax1_nodepth=bgA1,odds_E1=Pb[0]/bgA,odds_E1_conservative=Pb[0]/bgA1,odds_E3=Pb[1]/bgA)
# (d) bound window extended down to 0.02" (PMA excludes only ~50-90% at <=30 au)
fs2,_,_=fsep(0.02,1.24); out['fsep_0.02-1.24_over_0.1-1.24']=fs2/fs
# (e) no-q-window variant (any companion mass): odds rise by 0.9/0.5
print(json.dumps(out,indent=1))
json.dump(out,open(S+'results/priors_check.json','w'),indent=1)
