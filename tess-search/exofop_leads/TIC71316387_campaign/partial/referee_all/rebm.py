# Referee's own minimal detached-EB model (spheres, quadratic LD, circular orbit, ring integration)
import numpy as np
NR=300
_r=(np.arange(NR)+0.5)/NR; _dr=1.0/NR
def occ_frac(z,p,u1,u2):
    """fraction of flux of unit-radius quadratic-LD star hidden by disk radius p at centre distance z (arrays z)"""
    z=np.atleast_1d(z)[:,None]; r=_r[None,:]
    mu=np.sqrt(1-_r**2)
    I=1-u1*(1-mu)-u2*(1-mu)**2
    w=I*2*np.pi*_r*_dr; tot=w.sum()
    with np.errstate(invalid='ignore',divide='ignore'):
        c=(r**2+z**2-p**2)/(2*r*np.maximum(z,1e-12))
    frac=np.arccos(np.clip(c,-1,1))/np.pi
    frac=np.where(z+r<=p,1.0,frac)          # ring entirely inside occulter
    frac=np.where((r+p<=z)|(z+p<=r),0.0,frac)
    return (frac*w[None,:]).sum(1)/tot
def pair_flux(t,T0,P,r1,k,cosi,J,ua,ub,dts=0.0):
    """normalised flux of the pair (=1 out of eclipse). r1=R1/a, k=R2/R1, J=disk-averaged surface brightness ratio"""
    ph=2*np.pi*(t-T0)/P
    # secondary eclipse time shift dts (days): shift phase locally for points near phase 0.5
    near_sec=np.cos(ph)<0
    ph=np.where(near_sec,2*np.pi*(t-T0-dts)/P,ph)
    sini=np.sqrt(1-cosi**2)
    d=np.sqrt(np.sin(ph)**2+(cosi*np.cos(ph))**2)/r1   # in units of R1
    L1=1.0; L2=k**2*J
    F=np.full(t.shape,L1+L2)
    m1=(~near_sec)&(d<1+k)
    if m1.any(): F[m1]-=L1*occ_frac(d[m1],k,*ua)
    m2=near_sec&(d<1+k)
    if m2.any(): F[m2]-=L2*occ_frac(d[m2]/k,1.0/k,*ub)
    return F/(L1+L2)
