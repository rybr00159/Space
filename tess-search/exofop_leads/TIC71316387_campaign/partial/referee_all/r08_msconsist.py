# Referee: which pair light fraction f makes the (independent) geometric profile solutions consistent with main-sequence M dwarfs?
# Mass(M_Ks): Mann+2019; Radius(M_Ks): Mann+2015 eq.4; T-Ks colour vs M_Ks: Pecaut&Mamajek sequence (lcmodel table, colour only)
import sys, json, glob
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib'); import ms
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all/'
a=[-0.642,-0.208,-8.43e-4,7.87e-3,1.42e-4,-2.13e-4]
def Mass(MK): x=MK-7.5; return 10**sum(c*x**i for i,c in enumerate(a))
def Rad(MK): return 1.9515-0.3520*MK+0.01680*MK**2
mg=np.linspace(0.1,0.85,300); MKg=ms.MKs(mg); MTg=ms.MT(mg)
o=np.argsort(MKg); MKg,MTg=MKg[o],MTg[o]
def MT_of(MK): return np.interp(MK,MKg,MTg)
def SB(MK): return 10**(-0.4*MT_of(MK))/Rad(MK)**2
P=1.30367228*86400; G=6.674e-8; rhosun=1.408
res=[]
pts=[]
for f in sorted(glob.glob(B+'r04_profile_*.json')):
    for r in json.load(open(f))['res']: pts.append((r['f'],r['chi2'],r['x']))
pts.sort()
cmin=min(p[1] for p in pts)
for fEB,c2,x in pts:
    r1,k,cosi,J=x[:4]
    # solve MK1, MK2: J_MS(MK1,MK2)=J ; density of star1 (given q) = rho from r1
    def solve():
        MK2_of=lambda MK1: brentq(lambda MK2: SB(MK2)/SB(MK1)-J, MK1, 9.3)
        def g(MK1):
            MK2=MK2_of(MK1); q=Mass(MK2)/Mass(MK1)
            rho=3*np.pi/(G*P**2)/r1**3/(1+q)/rhosun
            return Mass(MK1)/Rad(MK1)**3-rho
        MK1=brentq(g,4.6,7.0); MK2=MK2_of(MK1)
        return MK1,MK2
    try:
        MK1,MK2=solve()
    except Exception as e:
        res.append(dict(f=fEB,dchi2=c2-cmin,err=str(e))); continue
    kMS=Rad(MK2)/Rad(MK1)
    MTp=-2.5*np.log10(10**(-0.4*MT_of(MK1))+10**(-0.4*MT_of(MK2)))
    Tpair=11.358-2.5*np.log10(fEB); d=10**((Tpair-MTp-0.02)/5+1)
    res.append(dict(f=fEB,dchi2=round(c2-cmin,2),k_fit=round(k,3),k_MS=round(kMS,3),M1=round(Mass(MK1),3),R1=round(Rad(MK1),3),M2=round(Mass(MK2),3),
                    R2=round(Rad(MK2),3),J=round(J,3),b=round(cosi/r1,3),dist_pc=round(d,1)))
for r in res: print(r)
json.dump(res,open(B+'r08_msconsist.json','w'),indent=1)
