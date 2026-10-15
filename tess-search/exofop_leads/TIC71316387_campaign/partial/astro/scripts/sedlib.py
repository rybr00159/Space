"""SED utilities: BT-Settl CIFIST synthetic surface magnitudes (SVO filters), extinction (CCM89, R_V=3.1),
Mamajek dwarf sequence for companion pairs."""
import numpy as np, glob, json, re
M='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/models/'
EEM='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/ext/EEM_dwarf.txt'
FILT=np.load(M+'filters.npy',allow_pickle=True).item()
RSUN=6.957e10; PC=3.0857e18
def _load(logg):
    T=[];mags=[];wl=None;fl=[]
    for f in sorted(glob.glob(M+f'bts_*_{logg:.1f}.npz')):
        z=np.load(f); T.append(int(z['teff'])); mags.append(json.loads(str(z['mags']))); wl=z['wl']; fl.append(z['flux'])
    return np.array(T),mags,wl,np.array(fl)
G45=_load(4.5); G50=_load(5.0)
def surf_mag(band,teff,system='Vega'):
    """surface magnitude interpolated linearly in Teff; logg 4.5 for Teff>=4500, 5.0 below"""
    teff=np.atleast_1d(np.asarray(teff,float)); out=np.empty_like(teff)
    for i,t in enumerate(teff):
        T,mags,_,_=G45 if t>=4500 else G50
        v=np.array([m[band][system] for m in mags])
        out[i]=np.interp(t,T,v)
    return out
def surf_spec(teff):
    T,_,wl,fl=G45 if teff>=4500 else G50
    j=np.searchsorted(T,teff); j=min(max(j,1),len(T)-1); w=(teff-T[j-1])/(T[j]-T[j-1])
    return wl,(1-w)*fl[j-1]+w*fl[j]
def dist_mod_R(R_rsun,d_pc): return -5*np.log10(R_rsun*RSUN/(d_pc*PC))
def ccm_AlamAv(lam_A,Rv=3.1):
    x=1e4/np.asarray(lam_A,float)  # 1/um
    a=np.zeros_like(x); b=np.zeros_like(x)
    ir=x<1.1; a[ir]=0.574*x[ir]**1.61; b[ir]=-0.527*x[ir]**1.61
    op=(x>=1.1)&(x<3.3); y=x[op]-1.82
    a[op]=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
    b[op]=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    uv=(x>=3.3)&(x<8); xx=x[uv]; Fa=np.where(xx>5.9,-0.04473*(xx-5.9)**2-0.009779*(xx-5.9)**3,0); Fb=np.where(xx>5.9,0.2130*(xx-5.9)**2+0.1207*(xx-5.9)**3,0)
    a[uv]=1.752-0.316*xx-0.104/((xx-4.67)**2+0.341)+Fa; b[uv]=-3.090+1.825*xx+1.206/((xx-4.62)**2+0.263)+Fb
    return a+b/Rv
def A_band(band,ebv,Rv=3.1): return ccm_AlamAv(FILT[band]['leff'])*Rv*ebv
def mag(band,teff,R,d,ebv,system='Vega'):
    return surf_mag(band,teff,system)+dist_mod_R(R,d)+A_band(band,ebv)
def addmag(*ms): return -2.5*np.log10(np.sum([10**(-0.4*np.asarray(m)) for m in ms],axis=0))
def eem():
    rows=[]
    for line in open(EEM):
        if line.startswith('#') or not line.strip(): continue
        p=line.split()
        if not re.match(r'^[OBAFGKMLTY]\d(\.\d)?V$',p[0]): continue
        try: rows.append((p[0],float(p[1]),float(p[6]) if p[6]!='...' else np.nan,float(p[-2]) if p[-2]!='...' else np.nan))
        except Exception: pass
    return rows
