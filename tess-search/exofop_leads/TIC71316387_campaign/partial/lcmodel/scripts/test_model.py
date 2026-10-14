import numpy as np, batman, time, sys
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
from ebmodel import occ_frac, _xg, _wg, eb_flux
# 1) occ_frac vs batman for p<1
u1,u2=0.45,0.2
worst=0
for p in [0.05,0.3,0.6,0.9]:
    zz=np.linspace(0,1+p,301)
    pr=batman.TransitParams(); pr.t0=0; pr.per=1.; pr.rp=p; pr.a=1e6; pr.inc=90; pr.ecc=0; pr.w=90; pr.limb_dark='quadratic'; pr.u=[u1,u2]
    # batman needs times; emulate z via a huge a: use internal function
    from batman import _quadratic_ld
    fb=_quadratic_ld._quadratic_ld(zz, p, u1, u2, 1)
    fm=np.array([1-occ_frac(z,p,u1,u2,_xg,_wg) for z in zz])
    d=np.max(np.abs(fb-fm)); worst=max(worst,d)
    print(f'p={p}: max |batman-mine| = {d:.2e}  (max depth {1-fb.min():.4f})')
# 2) p>1 vs brute force 2D grid
def brute(z,p,u1,u2,N=3000):
    x=np.linspace(-1,1,N); X,Y=np.meshgrid(x,x); R2=X**2+Y**2; inn=R2<1
    mu=np.sqrt(np.clip(1-R2,0,1)); I=(1-u1*(1-mu)-u2*(1-mu)**2)*inn
    cov=((X-z)**2+Y**2)<p*p
    return (I*cov).sum()/I.sum()
for p,z in [(1.5,1.2),(1.5,2.0),(2.0,1.5),(1.2,0.5),(1.2,2.1)]:
    print(f'p={p} z={z}: mine {occ_frac(z,p,u1,u2,_xg,_wg):.6f} brute {brute(z,p,u1,u2):.6f}')
# 3) speed
t=np.linspace(0,1.3,20000)
t0=time.time(); f=eb_flux(t,0.,1.3036712,0.2,0.7,0.1,0.45,0.9,0.,0.,0.45,0.2,0.3,0.3,200/86400,3,_xg,_wg); 
t0=time.time(); f=eb_flux(t,0.,1.3036712,0.2,0.7,0.1,0.45,0.9,0.,0.,0.45,0.2,0.3,0.3,200/86400,3,_xg,_wg); print('20000 pts x3 sub: %.3f s'%(time.time()-t0))
# eccentric sanity: secondary displaced by ~ 2P/pi e cos w
for ecw in [0.0,0.01]:
    t=np.linspace(0.3,1.0,70001); f=eb_flux(t,0.,1.3036712,0.2,0.7,0.1,0.45,0.0,ecw,0.0,0.45,0.2,0.3,0.3,1e-5,1,_xg,_wg)
    ts=t[np.argmin(f)]; print('ecosw',ecw,'secondary at phase %.5f expected %.5f'%(ts/1.3036712, 0.5+2/np.pi*ecw/1.0*1))
