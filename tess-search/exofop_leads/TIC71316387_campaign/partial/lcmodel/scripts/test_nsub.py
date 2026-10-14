import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit2
from ebmodel import eb_flux
t=np.linspace(-0.06,0.06,2001)
for par in [(0.2,0.7,0.16,0.45,0.9),(0.25,0.9,0.22,0.45,0.95),(0.15,0.5,0.11,0.45,0.5)]:
    rsum,k,cosi,J,L3=par
    for texp,ns in [(200/86400,1),(200/86400,2),(200/86400,3),(600/86400,3),(600/86400,5),(1800/86400,9),(1800/86400,15)]:
        ref=eb_flux(t,0.,1.3036712,rsum,k,cosi,J,L3,0.,0.,0.3,0.3,0.3,0.3,texp,41,fit2.XG,fit2.WG)
        m=eb_flux(t,0.,1.3036712,rsum,k,cosi,J,L3,0.,0.,0.3,0.3,0.3,0.3,texp,ns,fit2.XG,fit2.WG)
        print(par, round(texp*86400), ns, 'depth %.0f ppm'%((1-ref.min())*1e6), 'max err %.1f ppm'%(np.abs(m-ref).max()*1e6))
