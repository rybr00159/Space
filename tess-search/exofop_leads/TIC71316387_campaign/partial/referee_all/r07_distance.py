# Referee: photometric distance / contrast of the M-dwarf pair using an independent mass-luminosity relation (Mann+2019 M_Ks)
import numpy as np, json
from scipy.optimize import brentq
a=[-0.642,-0.208,-8.43e-4,7.87e-3,1.42e-4,-2.13e-4]   # Mann et al. 2019, log10 M vs (M_Ks-7.5), n=5, no [Fe/H]
def logM(MK): x=MK-7.5; return sum(c*x**i for i,c in enumerate(a))
def MKs_of(M): return brentq(lambda mk: logM(mk)-np.log10(M),4.0,11.0)
out={}
Ks_blend=10.341; T_blend=11.358; mu_G=5*np.log10(1000/3.3622/10)
# T-Ks colours of dwarfs (empirical, approx.): M2V ~2.40, M4.5V ~2.95 (Mamajek: V-Ks, V-Ic and T~Ic+0.25) -- referee's adopted values
TK={0.44:2.40,0.19:2.95}
for (M1,M2,fEB,lab) in [(0.44,0.19,0.0282,'msfree'),(0.47,0.207,0.0270,'msbound')]:
    k1,k2=MKs_of(M1),MKs_of(M2)
    MK=-2.5*np.log10(10**(-0.4*k1)+10**(-0.4*k2))
    MT=-2.5*np.log10(10**(-0.4*(k1+TK[0.44]))+10**(-0.4*(k2+TK[0.19])))
    Tpair=T_blend-2.5*np.log10(fEB)
    d=10**((Tpair-MT-0.02)/5+1)
    Kpair=MK+mu_G
    Kg=-2.5*np.log10(10**(-0.4*Ks_blend)-10**(-0.4*Kpair))
    out[lab]=dict(MKs1=k1,MKs2=k2,MKs_pair=MK,MT_pair=MT,T_pair=Tpair,dist_pc=d,dKs_at_297pc=Kpair-Kg)
print(json.dumps(out,indent=1))
json.dump(out,open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all/r07_distance.json','w'),indent=1)
