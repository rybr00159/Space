# Final referee f03: power of a common-proper-motion test (field M-dwarf PM density near the G star's PM)
import numpy as np, json
from astropy.table import Table
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
t=Table.read(B+'synthesis/field_Mdwarfs_1p5deg.ecsv')
br=np.array(t['bp_rp'],float); G=np.array(t['phot_g_mean_mag'],float)
pmra=np.array(t['pmra'].filled(np.nan),float); pmde=np.array(t['pmdec'].filled(np.nan),float); rv=np.array(t['radial_velocity'].filled(np.nan),float)
Tm=G-0.00522555*br**3+0.0891337*br**2-0.633923*br+0.0324473; f=10**(-0.4*(Tm-11.358))
dmu=np.hypot(pmra-38.503,pmde+3.134)
out={}
for lab,m in {'2.1-2.5 window':(br>=2.1)&(br<=2.5)&(f>=0.022)&(f<=0.065),'1.8-2.7 window':(br>=1.8)&(br<=2.7)&(f>=0.022)&(f<=0.065),
              'all BP-RP 1.8-2.8 G14-18.5':np.ones(len(t),bool)}.items():
    m=m&np.isfinite(dmu); n=int(m.sum())
    o={'n':n}
    for R in (3,5,10,15):
        k=int((dmu[m]<R).sum()); o[f'k_lt{R}']=k
    # smooth estimate: fraction within 10 and 15 mas/yr scaled by area to 3 and 5 mas/yr (locally uniform density)
    for R in (3,5):
        o[f'frac_lt{R}_scaled_from15']=float((dmu[m]<15).mean()*(R/15)**2)
        o[f'frac_lt{R}_UL95_direct']=float(1-0.05**(1/max(n,1))) if o[f'k_lt{R}']==0 else None
    out[lab]=o
# LR with bound P(dmu<R) ~ 0.999 (R=3, soft) / ~0.9 (R=3, empirical B, 84th pct 3.1) 
for lab,o in out.items():
    o['LR_dmu<5_smooth']=0.99/o['frac_lt5_scaled_from15'] if o['frac_lt5_scaled_from15']>0 else None
    o['LR_dmu<3_smooth']=0.85/o['frac_lt3_scaled_from15'] if o['frac_lt3_scaled_from15']>0 else None
mrv=np.isfinite(rv)&(br>=1.8)
out['RV']=dict(n=int(mrv.sum()),frac_within3=float((np.abs(rv[mrv]+6.07)<3).mean()),frac_within5p5=float((np.abs(rv[mrv]+6.07)<5.5).mean()))
json.dump(out,open(B+'final_referee/f03_cpm.json','w'),indent=1)
print(json.dumps(out,indent=1))
