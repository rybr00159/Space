"""r04: own Monte-Carlo population of possible unresolved bound companions of TIC 214311908 (and their constraints).
Prior draws (re-weighted in r05 for other hyper-parameters): log10 a ~ N(1.30, 1.16) (Winters+19-like; reweighted),
e ~ U(0, 0.8), isotropic orientation, uniform mean anomaly -> projected separation; q ~ U(0.075/m1, 1) (reweighted).
m1 solved so that the system total M_Ks = 6.4406 (Mann+19); m2 = q m1. Radii Mann+15. Band offsets (T, I, G, BP, RP
minus Ks) from the Mamajek table at each mass.
Constraint likelihoods stored separately (r05 switches them):
  L_ipd     : star/out/ipd_limits.json 'cells' (all calibration stars; n >= 5), LR per (sep, Delta I)
  L_ipd_red : 'red_variant' cells (calibrated on red stars only; coarser)
  L_ipd_far : as L_ipd but no constraint inside 0.25" (the calibration pairs are mostly nearer than 64 pc, so at a given
              angle they are at smaller physical separation and shorter period: their RUWE excess from orbital motion
              does not transfer to this target)
  L_ast     : star/out/astrometric_limits.json L_combined, mapped by the photocentre wobble |B - beta_G| onto the grid's
              0.37-Msun-primary branch (as fpp/f03 did; this is the only practical use of that grid)
  L_dec     : DECaPS2 residual-image 5-sigma limits (referee-corrected) in T-ish band: companion brighter than the limit ->
              0.02 within 1.5 mag of the limit, 1e-3 if brighter still (a >= 25-sigma source cannot be missed)
  L_wide    : sep > 2.5" and Delta G < 6 -> 1e-3 (Gaia DR3 would list it; no such co-moving source exists)
  CMD       : eps = x_obs - delta, delta = dM_bin(Ks) - slope x shift of BP-RP (r02 ridge); evaluated in r05
Outputs: out/r04_companions.npz"""
import numpy as np, sys, json
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
from scipy.interpolate import RegularGridInterpolator
rng = np.random.default_rng(31)
N = 300000
loga = rng.normal(1.30, 1.16, N); a = 10**loga
e = rng.uniform(0, 0.8, N)
cosi = rng.uniform(-1, 1, N); Om = rng.uniform(0, 2*np.pi, N); om = rng.uniform(0, 2*np.pi, N); M = rng.uniform(0, 2*np.pi, N)
E = M.copy()
for _ in range(40): E = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
nu = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
r = a*(1 - e*np.cos(E))
x = r*(np.cos(Om)*np.cos(om + nu) - np.sin(Om)*np.sin(om + nu)*cosi)
y = r*(np.sin(Om)*np.cos(om + nu) + np.cos(Om)*np.sin(om + nu)*cosi)
sep = np.hypot(x, y)/DIST
u = rng.uniform(0, 1, N)
m1 = np.full(N, 0.37)
for it in range(60):
    qmin = 0.075/m1; q = qmin + u*(1 - qmin); m2 = q*m1
    k1 = mks_of_mass(m1); k2 = mks_of_mass(m2)
    tot = -2.5*np.log10(10**(-0.4*k1) + 10**(-0.4*k2))
    m1 = mass_of_mks(k1 + (MKS_OBS - tot))
qmin = 0.075/m1; q = qmin + u*(1 - qmin); m2 = q*m1
k1 = mks_of_mass(m1); k2 = mks_of_mass(m2)
print('max |total MKs - obs|', float(np.max(np.abs(-2.5*np.log10(10**(-0.4*k1) + 10**(-0.4*k2)) - MKS_OBS))))
Rad1 = radius_of_mks(k1); Rad2 = radius_of_mks(k2)
rho1 = m1*MSUN/(4/3*np.pi*(Rad1*RSUN)**3); rho2 = m2*MSUN/(4/3*np.pi*(Rad2*RSUN)**3)
band = {}
for b_ in ('M_T', 'M_I', 'M_G', 'M_BP', 'M_RP'):
    band[b_] = (k1 + colour_offsets(m1, b_), k2 + colour_offsets(m2, b_))
dT = band['M_T'][1] - band['M_T'][0]; dI = band['M_I'][1] - band['M_I'][0]; dG = band['M_G'][1] - band['M_G'][0]
Dc = 1/(1 + 10**(0.4*dT))
dKs = k2 - k1
dMbin = -2.5*np.log10(1 + 10**(-0.4*dKs))
f = lambda mm: 10**(-0.4*mm)
BPs = -2.5*np.log10(f(band['M_BP'][0]) + f(band['M_BP'][1])); RPs = -2.5*np.log10(f(band['M_RP'][0]) + f(band['M_RP'][1]))
dcol = (BPs - RPs) - (band['M_BP'][0] - band['M_RP'][0])
# ---- IPD cells
ipdj = jload('star/out/ipd_limits.json')
def ipd_lr(cells, key='LR', min_sep=0.0):
    out = np.ones(N)
    for c in cells:
        if c.get('n', 0) < 5 or not (c[key] == c[key]): continue
        if c['sep'][1] <= min_sep: continue
        mm = (sep >= max(c['sep'][0], min_sep)) & (sep < c['sep'][1]) & (dI >= c['dm'][0]) & (dI < c['dm'][1])
        out[mm] = min(1.0, c[key])
    return out
L_ipd = ipd_lr(ipdj['cells']); L_ipd84 = ipd_lr(ipdj['cells'], 'LR_hi84')
L_ipd_red = ipd_lr(ipdj['red_variant']['cells']); L_ipd_far = ipd_lr(ipdj['cells'], min_sep=0.25)
# ---- astrometric grid, wobble-mapped
ast = jload('star/out/astrometric_limits.json')
AM = np.array(ast['M2']); AA = np.log10(np.array(ast['a_au'])); AL = np.array(ast['L_combined'])
def lightG(m): return f(mks_of_mass(m) + colour_offsets(m, 'M_G'))
dg = AM/(0.37 + AM) - lightG(AM)/(lightG(0.37) + lightG(AM))
ipk = int(np.argmax(dg))
fi = RegularGridInterpolator((AM, AA), AL, bounds_error=False, fill_value=None)
Bm = m2/(m1 + m2); bG = f(band['M_G'][1])/(f(band['M_G'][0]) + f(band['M_G'][1])); wob = np.abs(Bm - bG)
M2eq = np.where(wob > dg[ipk], AM[ipk], np.interp(wob, dg[:ipk+1], AM[:ipk+1]))
L_ast = fi(np.column_stack([np.clip(M2eq, AM[0], AM[-1]), np.clip(loga, AA[0], AA[-1])]))
L_ast[loga > AA[-1]] = 1.0; L_ast[loga < AA[0]] = 1.0
L_ast = np.clip(L_ast, 1e-4, 1.0)
# ---- DECaPS2 (referee-corrected i/z max at 1-2", investigator beyond 2.5"; weak 2.0 at 0.75")
ds_ = np.array([0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])
dv_ = np.array([2.0, 4.24, 5.07, 5.80, 5.58, 6.17, 6.89, 7.53, 7.88, 8.59])
lim = np.where(sep >= 0.75, np.interp(sep, ds_, dv_), -99.)
L_dec = np.where(dT < lim - 1.5, 1e-3, np.where(dT < lim, 0.02, 1.0))
L_wide = np.where((sep > 2.5) & (dG < 6), 1e-3, 1.0)
# Holman & Wiegert (1999) S-type stability for a planet at P = 5.69 d around each component
def hw(mu, ee): return 0.464 - 0.380*mu - 0.631*ee + 0.586*mu*ee + 0.150*ee**2 - 0.198*mu*ee**2
ap1 = (G*m1*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/AU; ap2 = (G*m2*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/AU
stab1 = ap1 < hw(m2/(m1 + m2), e)*a; stab2 = ap2 < hw(m1/(m1 + m2), e)*a
# chromatic: companion light fraction in V and I (for what-if tests)
MV1 = k1 + colour_offsets(m1, 'Mv'); MV2 = k2 + colour_offsets(m2, 'Mv')
DcV = 1/(1 + 10**(0.4*(MV2 - MV1))); DcI = 1/(1 + 10**(0.4*dI))
np.savez(f'{R2}/out/r04_companions.npz', a=a, loga=loga, e=e, sep=sep, u=u, q=q, m1=m1, m2=m2, R1=Rad1, R2=Rad2, rho1=rho1,
         rho2=rho2, dT=dT, dI=dI, dG=dG, dKs=dKs, Dc=Dc, dMbin=dMbin, dcol=dcol, L_ipd=L_ipd, L_ipd84=L_ipd84, L_ipd_red=L_ipd_red,
         L_ipd_far=L_ipd_far, L_ast=L_ast, L_dec=L_dec, L_wide=L_wide, stab1=stab1, stab2=stab2, DcV=DcV, DcI=DcI)
W = L_ipd*L_ast*L_dec*L_wide
print('surviving prior mass: ipd %.3f ast %.3f dec %.3f wide %.3f all %.3f; red-ipd all %.3f' % (L_ipd.mean(), L_ast.mean(), L_dec.mean(), L_wide.mean(), W.mean(),
      (L_ipd_red*L_ast*L_dec*L_wide).mean()))
for lo, hi in [(0, 0.3), (0.3, 0.75), (0.75, 1.5), (1.5, 2.5), (2.5, 4.5), (4.5, 99)]:
    mm = (dT >= lo) & (dT < hi)
    print(f'dT {lo}-{hi}: prior {mm.mean():.3f}, surviving {np.mean(W*mm):.4f}')
for lo, hi in [(0, 0.5), (0.5, 4), (4, 6), (6, 65), (65, 200), (200, 1e6)]:
    mm = (a >= lo) & (a < hi)
    print(f'a {lo}-{hi} au: prior {mm.mean():.3f}, surviving {np.mean(W*mm):.4f}')
