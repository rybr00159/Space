"""r09: check of the analyst's effective localisation area for uncatalogued sources: A_eff = sum over the referee's real
combined Delta-chi2 sky map (0.5" grid, co-moving with the M3V) of exp(-Delta chi2/(2 kappa)) x 0.25 arcsec^2, and the
expected number of uncatalogued T <= 17.47 / T <= 16 stars inside it before any archival-imaging veto
(densities from star/out/blend_density.json). Output: out/r09_aeff.json"""
import numpy as np, sys, json
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
z = np.load(f'{B}/referee_loc/out/inj_base.npz')
mp = z['real_gsum'].reshape(len(z['GN']), len(z['GE'])).astype(float)
EE, NN = np.meshgrid(z['GE'], z['GN']); RR = np.hypot(EE, NN)
mp = mp - mp[np.unravel_index(np.argmin(RR), RR.shape)]   # relative to the M3V position
bd = jload('star/out/blend_density.json')['gaia_density_per_arcsec2']
out = {}
for kap in (1.0, 2.2, 4.8):
    L = np.exp(-np.clip(mp, -50, None)/(2*kap))
    A = float(L.sum()*0.25); A3 = float(L[RR < 3].sum()*0.25)
    out[f'kappa_{kap}'] = dict(A_eff=A, A_eff_r_lt_3=A3, N_T_le_17p47=A*bd['T<=17.47'], N_T_le_16=A*bd['T<=16'])
    print(kap, out[f'kappa_{kap}'])
out['map_min_rel_M3V'] = float(mp.min())
jdump(out, 'out/r09_aeff.json')
