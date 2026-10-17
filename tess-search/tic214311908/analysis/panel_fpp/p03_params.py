# Independent recomputation of planet/star derived quantities (own code; constants IAU 2015 nominal)
import numpy as np, json
G=6.67430e-11; Msun=1.98840987e30; Rsun=6.957e8; au=1.495978707e11; RE=6.3781e6; ME=5.9722e24; MJ=1.89813e27
day=86400.
P=5.691506; k0=0.1185; b0=0.52
Rs=0.383; Ms=0.370; Teff=3350.; J=11.36; Ks=10.48
def derived(Rs,Ms,P,k,b,Teff):
    a=(G*Ms*Msun*(P*day)**2/(4*np.pi**2))**(1/3)
    aR=a/(Rs*Rsun); inc=np.degrees(np.arccos(b/aR))
    T14=P/np.pi*np.arcsin(np.sqrt((1+k)**2-b**2)/(aR*np.sin(np.radians(inc))))*24
    T23=P/np.pi*np.arcsin(np.sqrt((1-k)**2-b**2)/(aR*np.sin(np.radians(inc))))*24
    Teq=Teff*np.sqrt(Rs*Rsun/(2*a)); L=Rs**2*(Teff/5772.)**4; S=L/(a/au)**2
    rho=Ms*Msun/(4/3*np.pi*(Rs*Rsun)**3)/1000.
    Rp=k*Rs*Rsun/RE
    return dict(a_au=a/au,aR=aR,inc=inc,T14_h=T14,T23_T14=T23/T14,Teq=Teq,L=L,S=S,rho_cgs=rho,Rp=Rp)
out={}
for R in [0.381,0.383,0.385]:
    out['R%.3f'%R]=derived(R,Ms,P,k0,b0,Teff)
def K_ms(Mp_ME,Ms,P,inc=88.8,e=0):
    return (2*np.pi*G/(P*day))**(1/3)*Mp_ME*ME*np.sin(np.radians(inc))/((Ms*Msun+Mp_ME*ME)**(2/3))/np.sqrt(1-e**2)
for M in [16.6,21.0,21.34,28.3]:
    out['K_M%.2f'%M]=K_ms(M,Ms,P)
def tsm(Rp,Mp,Teq,Rs,J):
    sc=0.190 if Rp<1.5 else (1.26 if Rp<2.75 else (1.28 if Rp<4.0 else 1.15))
    return sc*Rp**3*Teq/(Mp*Rs**2)*10**(-J/5)
def planck(lam,T):
    h=6.62607015e-34;c=2.99792458e8;kB=1.380649e-23
    return 2*h*c**2/lam**5/(np.exp(h*c/(lam*kB*T))-1)
def esm(k,Teq,Teff,Kmag):
    lam=7.5e-6
    return 4.29e6*planck(lam,1.10*Teq)/planck(lam,Teff)*k**2*10**(-Kmag/5)
d=out['R0.383']
for Mp in [16.6,21.0,21.34,28.3]:
    out['TSM_M%.2f'%Mp]=tsm(d['Rp'],Mp,d['Teq'],0.383,J)
out['ESM']=esm(k0,d['Teq'],Teff,Ks)
# Monte Carlo with CK17-like mass scatter not attempted; propagate star + k
rng=np.random.default_rng(1); N=200000
Rs_s=rng.normal(0.383,0.012,N); Ms_s=rng.normal(0.370,0.011,N); T_s=rng.normal(3350,55,N)
dep=rng.normal(15.8,0.9,N); k_s=k0*np.sqrt(dep/15.8)
b_s=rng.uniform(0.49,0.54,N)
dd=derived(Rs_s,Ms_s,P,k_s,b_s,T_s)
for key in ['Rp','Teq','aR','a_au','S','T14_h','rho_cgs']:
    q=np.percentile(dd[key],[16,50,84]); out['MC_'+key]=list(q)
# twin host: Mann+15 R(MKs), each component 0.753 mag fainter in Ks
dm=5*np.log10(64.31)-5; MK=Ks-dm
R15=lambda m:1.9515-0.3520*m+0.01680*m**2
out['MKs']=MK; out['R_Mann15_single']=R15(MK); out['R_Mann15_twincomp']=R15(MK+0.753)
# depth doubling for equal twin -> k*sqrt(2) approx
out['Rp_twin_approx']=k0*np.sqrt(2)*R15(MK+0.753)*Rsun/RE
# insolation / Teq check, distance
json.dump(out,open('out/p03_params.json','w'),indent=1,default=float)
for kk,v in out.items(): print(kk, v)
