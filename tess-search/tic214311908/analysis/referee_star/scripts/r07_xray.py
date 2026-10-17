"""Referee: eRASS1 X-ray source near the M3V (VizieR J/A+A/682/A34, main catalogue), own query; positional
association odds with the M3V, the pair A+B, the background star 3.1" NE and the Mira; L_X and L_X/L_bol.
eRASS1 RADEC_ERR is the combined 1-sigma positional error (sqrt(sig_RA^2+sig_Dec^2)); we use per-axis
sigma = RADEC_ERR/sqrt(2) as the nominal and also add the recommended systematic.
Output: out/r07_xray.json, data/erass1_cone.csv"""
import json, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
R = 'tic214311908/analysis/referee_star'
v = Vizier(columns=['**', '_r'], row_limit=-1)
c0 = SkyCoord(260.969674, -51.323821, unit='deg')
res = v.query_region(c0, radius=120 * u.arcsec, catalog='J/A+A/682/A34')
print(res)
t = res[0].to_pandas(); t.to_csv(f'{R}/data/erass1_cone.csv', index=False)
print(t.T.to_string()[:6000])
# ---------------- association
r = t.iloc[0]
ep = 2000.0 + (r.MJD + 2400000.5 - 2451545.0) / 365.25
g = pd.read_csv(f'{R}/data/gaia_cone45.csv')
def at_ep(row):
    dt = ep - 2016.0
    return row.ra + np.nan_to_num(row.pmra) * dt / 3.6e6 / np.cos(np.radians(row.dec)), row.dec + np.nan_to_num(row.pmdec) * dt / 3.6e6
cands = {'M3V': 5925373272850784768, 'pair A': 5925373272850782208, 'pair B': 5925373272827253632,
         'bkg 3.1" NE': 5925373272827260544, 'star 5.7"': 5925373272827262208}
xr = SkyCoord(r.RA_ICRS, r.DE_ICRS, unit='deg')
out = dict(IAUName=str(r.IAUName), epoch=ep, DET_LIKE=float(r.DetLike0), counts=[float(r.MLcts1), float(r.s_MLcts1)], flux=[float(r.MLFlux1), float(r.s_MLFlux1)],
           POS_ERR=float(r.posErr), RADEC_ERR=float(r.RADEErr), APE_POIS=float(r.APEPois1), bkg_counts=float(r.MLBkg1))
seps = {}
for k, sid in cands.items():
    row = g[g.source_id == sid].iloc[0]; a, d = at_ep(row); s = xr.separation(SkyCoord(a, d, unit='deg')).arcsec; seps[k] = s
# the Mira is outside the 45" Gaia cone? it is at 17.8": include
row = g[g.source_id == 5925373268528684288].iloc[0]; a, d = at_ep(row); seps['Mira'] = xr.separation(SkyCoord(a, d, unit='deg')).arcsec
out['separations_arcsec'] = seps
for lab, sig in (('per-axis sigma = POS_ERR/sqrt2', r.posErr / np.sqrt(2)), ('per-axis sigma = POS_ERR (investigator)', r.posErr)):
    L = {k: float(np.exp(-0.5 * (s / sig) ** 2)) for k, s in seps.items()}
    out[lab] = dict(sigma=float(sig), rel_likelihood=L, norm_sep={k: float(s / sig) for k, s in seps.items()})
    print(lab, f'sigma {sig:.2f}"', {k: f'{s:.2f}" ({s/sig:.2f} sig, L={L[k]:.3g})' for k, s in seps.items()})
d_pc = 1000 / (15.4973 + 0.0521)
LX = 4 * np.pi * (d_pc * 3.0857e18) ** 2 * r.MLFlux1
Lbol = 0.0165 * 3.828e33
out['logLX'] = float(np.log10(LX)); out['LX_Lbol'] = float(LX / Lbol)
# Wright+11 / Wright+18 rotation-activity: Rx = Rx_sat (Ro/Ro_sat)^beta for Ro > Ro_sat; Rx_sat = 10^-3.05 (W18), Ro_sat 0.14, beta -2.3 (W18 fully conv.)?
for lab, (Rsat, Rosat, beta) in {'Wright11': (10 ** -3.13, 0.13, -2.70), 'Wright18': (10 ** -3.05, 0.14, -2.3)}.items():
    Rx = LX / Lbol
    Ro = Rosat * (Rx / Rsat) ** (1 / beta) if Rx < Rsat else Rosat
    tau = 10 ** (2.33 - 1.5 * 0.37 + 0.31 * 0.37 ** 2)     # Wright+18 tau(M)
    out[f'Ro_{lab}'] = float(Ro); out[f'Prot_{lab}'] = float(Ro * tau)
    print(lab, 'Ro', round(Ro, 3), 'tau', round(tau, 1), 'Prot', round(Ro * tau, 1))
print('log LX', round(np.log10(LX), 3), 'LX/Lbol', LX / Lbol)
json.dump(out, open(f'{R}/out/r07_xray.json', 'w'), indent=1, default=float)
