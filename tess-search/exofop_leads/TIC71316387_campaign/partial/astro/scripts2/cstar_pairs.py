"""Predicted Gaia BP/RP flux-excess change dC for MS pairs at 297 pc (A=1), if the pair is inside the BP/RP windows but
fraction g of its G flux is captured by the G-band PSF fit. Compare with observed C*_target - median C*_twins = +0.0242."""
import sys, numpy as np, json
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/scripts')
from sed_fit import *
d=1000/3.3622; B=['GAIA__GAIA3.Gbp','GAIA__GAIA3.G','GAIA__GAIA3.Grp','TESS__TESS.Red','Generic__Cousins.I','2MASS__2MASS.Ks']
I={'BP':10**(-0.4*(12.066102-25.3385)),'G':10**(-0.4*(11.775332-25.6874)),'RP':10**(-0.4*(11.280263-24.7479))}
obs=0.0204+0.0038
mG=gmodel(B,5800,1.148,d,0.010)
out={}
for pr in [('K7V','M2V'),('M0V','M3V'),('M0V','M4.5V'),('M1V','M5V'),('M1.5V','M5V'),('M2V','M5V'),('M2.5V','M5V'),('M3V','M5V')]:
    mp=pair_abs(pr,B)+5*np.log10(d/10)+np.array([TAB[b][4] for b in B])*0.010
    f=10**(-0.4*(mp-mG))   # flux relative to G star
    fBP,fG,fRP,fT,fI,fK=f
    C0=(I['BP']+I['RP'])/I['G']
    dC={g:(fBP*I['BP']+fRP*I['RP'])/I['G']-C0*g*fG for g in (0.0,0.5,1.0)}
    # amplitude A (pair luminosity scaling) needed to give obs with g=0 and g=0.5
    Aneed={g:obs/dC[g] if dC[g]>0 else None for g in (0.0,0.5)}
    out['|'.join(pr)]=dict(fTESS_frac=float(fT/(1+fT)),fRP=float(fRP),fBP=float(fBP),dC_g0=float(dC[0.0]),dC_g05=float(dC[0.5]),dC_g1=float(dC[1.0]),
        fTESS_needed_g0=float(fT*Aneed[0.0]/(1+fT*Aneed[0.0])),fTESS_needed_g05=(float(fT*Aneed[0.5]/(1+fT*Aneed[0.5])) if Aneed[0.5] else None))
    print(pr,{k:round(v,4) if isinstance(v,float) else v for k,v in out['|'.join(pr)].items()})
json.dump(out,open(D+'results/cstar_pairs.json','w'),indent=1)
