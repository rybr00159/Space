import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/referee_all')
import numpy as np, batman
from rebm import occ_frac
z=np.linspace(0,1.5,61); k=0.44; u=(0.3,0.3)
mine=occ_frac(z,k,*u)
pm=batman.TransitParams(); pm.t0=0; pm.per=1.; pm.rp=k; pm.a=10.; pm.inc=90.; pm.ecc=0; pm.w=90; pm.u=list(u); pm.limb_dark='quadratic'
# map z to time: z = a*sin(2pi t)
t=np.arcsin(z/10.)/(2*np.pi)
bm=batman.TransitModel(pm,t).light_curve(pm)
print('max |mine-batman| =',np.max(np.abs((1-mine)-bm)),' depth', 1-bm.min())
# big occulter p=1/k on small star check: total occultation should give 1
print(occ_frac(np.array([0.,0.5,1.0,1/k-1-0.01, 1/k+1+0.01]),1/k,0.3,0.3))
