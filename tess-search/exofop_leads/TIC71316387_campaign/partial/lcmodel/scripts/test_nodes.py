import numpy as np, sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
from ebmodel import occ_frac, eb_flux
import ebmodel
def brute(z,p,u1,u2,N=4000):
    x=np.linspace(-1,1,N); X,Y=np.meshgrid(x,x); R2=X**2+Y**2; inn=R2<1
    mu=np.sqrt(np.clip(1-R2,0,1)); I=(1-u1*(1-mu)-u2*(1-mu)**2)*inn
    cov=((X-z)**2+Y**2)<p*p
    return (I*cov).sum()/I.sum()
for NG in [16,24,32,64]:
    xg,wg=np.polynomial.legendre.leggauss(NG)
    worst=0
    for p in [0.3,0.6,0.8,1.25,1.6]:
        for z in np.linspace(abs(1-p)+1e-3, 1+p-1e-3, 25):
            ref=occ_frac(z,p,0.3,0.3,*np.polynomial.legendre.leggauss(200))
            v=occ_frac(z,p,0.3,0.3,xg,wg); worst=max(worst,abs(v-ref))
    t=np.linspace(-0.1,0.1,30000)
    eb_flux(t,0.,1.3,0.2,0.7,0.12,0.45,0.9,0.,0.,0.4,0.25,0.3,0.3,200/86400,3,xg,wg)
    t0=time.time()
    for i in range(5): eb_flux(t,0.,1.3,0.2,0.7,0.12,0.45,0.9,0.,0.,0.4,0.25,0.3,0.3,200/86400,3,xg,wg)
    print(NG,'max abs err vs 200-node',worst,'time',(time.time()-t0)/5)
print('brute check p=0.6', occ_frac(0.9,0.6,0.3,0.3,*np.polynomial.legendre.leggauss(200)), brute(0.9,0.6,0.3,0.3))
