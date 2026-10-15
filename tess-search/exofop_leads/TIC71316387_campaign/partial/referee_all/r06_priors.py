# Referee: independent prior-odds arithmetic (bound companion EB vs chance-aligned EB)
import numpy as np, json
from scipy.stats import norm
from astropy.table import Table
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d_pc=1000/3.3622; out={}
# --- bound: Raghavan+2010 lognormal logP(d) mu 5.03 sig 2.28; CSF 0.33+2*0.08+3*0.03
CSF=0.33+2*0.08+3*0.03
def fsep(s1,s2,Mtot):
    a1,a2=1.26*s1*d_pc,1.26*s2*d_pc
    lp=lambda a: np.log10(np.sqrt(a**3/Mtot)*365.25)
    return norm.cdf(lp(a2),5.03,2.28)-norm.cdf(lp(a1),5.03,2.28)
for Mt in (1.66,1.83):
    out[f'fsep_0.1-1.24_M{Mt}']=fsep(0.1,1.24,Mt)
out['fsep_0.48-1.24']=fsep(0.48,1.24,1.66)
qwin_any=0.5/0.9; qwin_match=(0.70-0.55)/0.9     # pair mass 0.63-0.68 over 1.03 Msun -> q 0.61-0.66; allow 0.55-0.70
Pb_any=CSF*fsep(0.1,1.24,1.66)*qwin_any
Pb_match=CSF*fsep(0.1,1.24,1.66)*qwin_match
out.update({'CSF':CSF,'Pbound_q0.5_1':Pb_any,'Pbound_qmatch':Pb_match})
# --- background from Gaia field counts (1-30 arcmin annulus)
T=Table.read(B+'astro/data/gaia_counts_30arcmin.ecsv')
sep=np.array(T['sep_arcmin']); G=np.array(T['phot_g_mean_mag'],float)
br=np.array(T['bp_rp'].filled(np.nan) if hasattr(T['bp_rp'],'filled') else T['bp_rp'],float)
k=sep>1.0; area=np.pi*(0.5**2-(1/60)**2)*3600**2   # arcsec^2
def Tmag(G,c):
    c=np.where(np.isfinite(c),c,1.0)
    return G-0.00522555*c**3+0.0891337*c**2-0.633923*c+0.0324473
Tm=Tmag(G,br)
circ=np.pi*1.24**2
Ttgt=11.358
# (i) any host bright enough that 5900 ppm needs <=50% undiluted depth
Tmax=Ttgt-2.5*np.log10(0.0059/0.5)
n_any=(k&(Tm<=Tmax)).sum(); out['n_hosts_any']=int(n_any); out['N_any_r1.24']=n_any/area*circ
# (ii) host that would give the observed pair: an M1.5-M2.5 eclipsed star (BP-RP 2.1-2.5) that supplies 2.2-6.5% of TESS light
fmin,fmax=0.022,0.065
Tlo=Ttgt-2.5*np.log10(fmax); Thi=Ttgt-2.5*np.log10(fmin)
m_match=k&(br>=2.1)&(br<=2.5)&(Tm>=Tlo)&(Tm<=Thi)
out['n_match']=int(m_match.sum()); out['N_match_r1.24']=m_match.sum()/area*circ; out['T_window']=[Tlo,Thi]
m_mdw=k&(br>=2.0)&(br<=2.7)&(G>=15)&(G<=17.8)
out['n_Mdwarf_localise_def']=int(m_mdw.sum()); out['N_Mdwarf_r1.24']=m_mdw.sum()/area*circ
# odds (EB occurrence and depth factors assumed equal for like-for-like comparisons)
out['odds_like_for_like_Mpair']=Pb_match/out['N_match_r1.24'] if out['n_match']>0 else None
# depth weighting for 'any host' comparison: bound needs undiluted depth >= 5900ppm/f; P(depth>=x) from Prsa+2022 quantiles (log-interp)
qs=np.array([0.1,0.25,0.5,0.75,0.9]); dq=np.array([0.004,0.015,0.077,0.272,0.442])
def frac_deeper(x):  # survival function by interpolation in log depth
    return float(np.interp(np.log10(x),np.log10(dq),1-qs,left=0.9,right=0.1))
for f in (0.028,0.05):
    out[f'depthfac_bound_f{f}']=frac_deeper(0.0059/f)
# background depth weighting: each host needs depth 5900ppm * Ftot/Fhost
hosts=k&(Tm<=Tmax)
need=0.0059*10**(-0.4*(Ttgt-Tm[hosts]))
wbg=np.mean([frac_deeper(x) for x in need]); out['depthfac_bg_mean']=wbg
for f in (0.028,0.05):
    out[f'odds_any_f{f}']=Pb_any*out[f'depthfac_bound_f{f}']/(out['N_any_r1.24']*wbg)
print(json.dumps(out,indent=1,default=float))
json.dump(out,open(B+'referee_all/r06_priors.json','w'),indent=1,default=float)
