import sys; sys.path.insert(0,'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
import pandas as pd
d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
d = d[np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A') & (d.ks_msigcom < 0.05) & np.isfinite(d.ruwe)].copy()
d['MKs'] = d.ks_m + 5*np.log10(d.parallax/100.)
d['GK']=d.phot_g_mean_mag-d.ks_m
icpt=7.279-1.531*2.756149
d['off']=d.MKs-(icpt+1.531*d.bp_rp)
b=d[(d.ruwe<1.1)&(abs(d.off)<0.25)]
c=np.polyfit(b.bp_rp,b.GK,1); d['rGK']=d.GK-np.polyval(c,d.bp_rp)
tGK=14.281043-10.482; trGK=tGK-np.polyval(c,2.756149)
print('GK slope vs bp_rp', c, 'target rGK', trGK)
cl=d[(abs(d.bp_rp-2.756)<0.1)&(d.ruwe<1.15)&(d.ipd_frac_multi_peak==0)&(d.ipd_gof_harmonic_amplitude<0.04)&(d.non_single_star==0)]
for lo,hi in [(0.4,1.5),(0.15,0.4),(-0.15,0.15),(-0.4,-0.15),(-0.6,-0.4),(-0.75,-0.6),(-0.9,-0.75),(-1.2,-0.9)]:
    s=cl[(cl.off>lo)&(cl.off<=hi)]
    print(f'off {lo:5.2f}..{hi:5.2f} n {len(s):4d} rGK median {np.median(s.rGK):+.3f} robust sd {1.4826*np.median(abs(s.rGK-np.median(s.rGK))):.3f}  ruwe med {np.median(s.ruwe):.3f}')
m1=0.253; qg=np.array([1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3])
k1=mks_of_mass(m1); k2=mks_of_mass(qg*m1)
f=lambda x:10**(-0.4*x)
def comp(col):
    a1=k1+colour_offsets(m1,col); a2=k2+colour_offsets(qg*m1,col); return -2.5*np.log10(f(a1)+f(a2)), a1
G,G1=comp('M_G'); BP,BP1=comp('M_BP'); RP,RP1=comp('M_RP'); K=-2.5*np.log10(f(k1)+f(k2))
dGK=(G-K)-(G1-k1); dbr=(BP-RP)-(BP1-RP1)
dx=(K-k1)-1.531*dbr; drGK=dGK-c[0]*dbr
for q,a,b_,c_ in zip(qg,dx,drGK,dbr): print(f'q {q:.1f}: offset {a:+.3f}  rGK {b_:+.3f}  dBP-RP {c_:+.3f}')
