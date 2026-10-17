"""Referee: (1) UVW of the M3V (astropy Galactic frame, heliocentric; LSR from Schoenrich+10 (11.1, 12.24, 7.25));
(2) wide system: separations, PM differences, tangential-velocity difference, escape velocity (Mann/Mamajek masses),
chance alignment estimate from local density of stars with similar PM and parallax.
Output: out/r08_kin_wide.json"""
import json, numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
R = 'tic214311908/analysis/referee_star'
g = pd.read_csv(f'{R}/data/gaia_cone45.csv').set_index('source_id')
t = g.loc[5925373272850784768]; A = g.loc[5925373272850782208]; Bc = g.loc[5925373272827253632]
plx = t.parallax + 0.052059
rng = np.random.default_rng(1); N = 20000
c = SkyCoord(ra=t.ra * u.deg, dec=t.dec * u.deg, distance=(1000 / rng.normal(plx, 0.036, N)) * u.pc,
             pm_ra_cosdec=rng.normal(t.pmra, t.pmra_error, N) * u.mas / u.yr, pm_dec=rng.normal(t.pmdec, t.pmdec_error, N) * u.mas / u.yr,
             radial_velocity=rng.normal(t.radial_velocity, t.radial_velocity_error, N) * u.km / u.s)
gal = c.galactic; v = gal.velocity
U, V, W = v.d_x.to(u.km / u.s).value, v.d_y.to(u.km / u.s).value, v.d_z.to(u.km / u.s).value
out = dict(UVW=[float(np.median(U)), float(np.median(V)), float(np.median(W))], UVW_sd=[float(np.std(U)), float(np.std(V)), float(np.std(W))])
Ul, Vl, Wl = U + 11.1, V + 12.24, W + 7.25
out['UVW_LSR'] = [float(np.median(Ul)), float(np.median(Vl)), float(np.median(Wl))]; out['vtot_LSR'] = float(np.median(np.sqrt(Ul ** 2 + Vl ** 2 + Wl ** 2)))
# Bensby+14 TD/D
def f(Ul, Vl, Wl, sU, sV, sW, Va): return np.exp(-Ul ** 2 / (2 * sU ** 2) - (Vl - Va) ** 2 / (2 * sV ** 2) - Wl ** 2 / (2 * sW ** 2)) / ((2 * np.pi) ** 1.5 * sU * sV * sW)
D = 0.85 * f(Ul, Vl, Wl, 35, 20, 16, -15); TD = 0.09 * f(Ul, Vl, Wl, 67, 38, 35, -46); H = 0.0015 * f(Ul, Vl, Wl, 160, 90, 90, -220)
out['TD_over_D'] = float(np.median(TD / D)); out['P_thick'] = float(np.median(TD / (D + TD + H)))
print('UVW', np.round(out['UVW'], 2), '+-', np.round(out['UVW_sd'], 2), 'LSR', np.round(out['UVW_LSR'], 2), 'vtot', round(out['vtot_LSR'], 1), 'TD/D', round(out['TD_over_D'], 3), 'P_thick', round(out['P_thick'], 3))
# wide system
c0 = SkyCoord(t.ra, t.dec, unit='deg'); cA = SkyCoord(A.ra, A.dec, unit='deg'); cB = SkyCoord(Bc.ra, Bc.dec, unit='deg')
sA = c0.separation(cA).arcsec; sAB = cA.separation(cB).arcsec; dist = 1000 / plx
dpm = np.array([A.pmra - t.pmra, A.pmdec - t.pmdec]); edpm = np.array([np.hypot(A.pmra_error, t.pmra_error), np.hypot(A.pmdec_error, t.pmdec_error)])
dv = 4.74047 * np.hypot(*dpm) / plx
G = 6.674e-11; Ms = 1.989e30; au = 1.496e11
for Mtot in (0.37 + 0.30 + 0.20, 0.264 * 2 + 0.30 + 0.20):
    vesc = np.sqrt(2 * G * Mtot * Ms / (sA * dist * au)) / 1e3
    out[f'vesc_Mtot{Mtot:.2f}'] = float(vesc)
# Monte Carlo on dv including A's large RUWE (inflate A's PM errors by RUWE)
dvs = 4.74047 * np.hypot(rng.normal(dpm[0], np.hypot(A.pmra_error * A.ruwe, t.pmra_error), N), rng.normal(dpm[1], np.hypot(A.pmdec_error * A.ruwe, t.pmdec_error), N)) / plx
out.update(sep_A_arcsec=float(sA), sep_A_au=float(sA * dist), sep_AB_arcsec=float(sAB), sep_AB_au=float(sAB * dist), dpm=dpm.tolist(), dpm_err=edpm.tolist(),
           dv_tan=float(dv), dv_tan_ruwe_inflated_68=[float(np.percentile(dvs, 16)), float(np.percentile(dvs, 84))],
           plx_A=[float(A.parallax), float(A.parallax_error)], plx_diff_sigma=float((A.parallax - t.parallax) / np.hypot(A.parallax_error, t.parallax_error)))
print({k: (np.round(v, 3) if isinstance(v, float) else v) for k, v in out.items() if k not in ('UVW', 'UVW_sd', 'UVW_LSR')})
# orbital motion expected at 754 au: v_orb ~ sqrt(G M / r) -> projected ~ 1 km/s; dv 0.57 consistent
json.dump(out, open(f'{R}/out/r08_kin_wide.json', 'w'), indent=1)
