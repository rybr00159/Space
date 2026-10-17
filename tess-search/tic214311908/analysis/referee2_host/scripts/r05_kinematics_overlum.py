"""r05: Does the frequency of over-luminous (metal-rich) single M dwarfs depend on kinematics, and what does that do
to the CMD-only twin probability for the M3V (which has hot/old-thin-disk kinematics, P_thick ~0.3)?
Label-free: uses only the clean field sample (host/data/field_mdwarfs.csv raw Gaia rows) and r01's local ridges.
(1) Target UVW (astropy, Gaia DR3 astrometry + RV) and v_tan.
(2) For clean field stars at the M3V colour (|dBP-RP|<0.08), in kinematic bins (v_tan; and |W|, v_tot for stars with
    Gaia RV): fraction with dMKs <= target, KDE density ratio n(x_t)/n(x_t + 0.753) (twin-component position), and the
    median r(G-Ks) of stars in the over-luminous window (x_t +- 0.1).
(3) CMD-only P(twin) = f_T * n(x_t+0.753)/n(x_t) for f_T = 0.8, 1.5, 2.35% per kinematic group.
Output: out/r05_kinematics.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
from scipy.stats import gaussian_kde
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['dMK'] = F.ks_m + dm - (R['MK'][0] * (F.bp_rp - C0) + R['MK'][1])
F['dGK'] = F.phot_g_mean_mag - F.ks_m - (R['GK'][0] * (F.bp_rp - C0) + R['GK'][1])
F['vtan'] = 4.74047 * np.hypot(F.pmra, F.pmdec) / F.parallax
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
D = F[clean].copy()
def uvw(ra, dec, plx, pmra, pmdec, rv):
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=(1000 / plx) * u.pc, pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                 radial_velocity=rv * u.km / u.s, frame='icrs').transform_to(Galactic())
    v = c.velocity.d_xyz.to(u.km / u.s).value
    # heliocentric U (to GC), V, W; add solar motion (Schoenrich+10) for LSR
    return v[0] + 11.1, v[1] + 12.24, v[2] + 7.25
tU, tV, tW = uvw(260.9696737, -51.3238212, 15.4973, -113.6712, 85.5436, 16.2837)
tvt = 4.74047 * np.hypot(-113.6712, 85.5436) / 15.4973
out = dict(target=dict(U_LSR=float(tU), V_LSR=float(tV), W_LSR=float(tW), vtot_LSR=float(np.sqrt(tU ** 2 + tV ** 2 + tW ** 2)), vtan_helio=float(tvt)))
print('target UVW (LSR):', out['target'])
hr = np.isfinite(D.radial_velocity)
U, V, W = uvw(D.ra[hr].values, D.dec[hr].values, D.parallax[hr].values, D.pmra[hr].values, D.pmdec[hr].values, D.radial_velocity[hr].values)
D.loc[hr, 'U'] = U; D.loc[hr, 'V'] = V; D.loc[hr, 'W'] = W; D['vtot'] = np.sqrt(D.U ** 2 + D.V ** 2 + D.W ** 2)
tx = r01['target']['cmd']['raw_plx']['dMK']; ty = r01['target']['colour_residuals']['dGK']; dx = 0.753
def stats(d, label):
    xv = d.dMK.values; k = gaussian_kde(xv, bw_method=0.06 / np.std(xv))
    nt, nc = float(k(tx)[0]), float(k(tx + dx)[0])
    win = d[(d.dMK - tx).abs() < 0.1]
    s = dict(label=label, n=int(len(d)), frac_as_bright=float(np.mean(xv <= tx)), dens_target=nt, dens_component=nc, ratio_comp_over_target=nc / nt,
             n_window=int(len(win)), window_median_rGK=float(np.median(win.dGK)) if len(win) else None,
             frac_window_rGK_ge_target=float(np.mean(win.dGK >= ty)) if len(win) else None)
    for fT in (0.008, 0.015, 0.0235):
        o = fT * nc / nt; s[f'P_twin_cmd_only_fT{fT}'] = o   # P = twin density / total density
    # bootstrap of the ratio
    rng = np.random.default_rng(3); bs = []
    for _ in range(200):
        xb = rng.choice(xv, len(xv)); kb = gaussian_kde(xb, bw_method=0.06 / np.std(xb)); bs.append(kb(tx + dx)[0] / kb(tx)[0])
    s['ratio_16_84'] = [float(v) for v in np.percentile(bs, [16, 84])]
    print(f"{label:38s} n={s['n']:5d} frac<=tx {s['frac_as_bright']:.4f} ratio {s['ratio_comp_over_target']:.2f} {s['ratio_16_84']} nwin {s['n_window']} "
          f"med r(G-Ks) {s['window_median_rGK']} frac>=t {s['frac_window_rGK_ge_target']} P_cmd(fT=1.5%) {s['P_twin_cmd_only_fT0.015']:.3f}")
    return s
res = [stats(D, 'all clean')]
for lo, hi in ((0, 25), (25, 40), (30, 60), (40, 70), (35, 55), (70, 400)):
    res.append(stats(D[D.vtan.between(lo, hi)], f'vtan {lo}-{hi}'))
dr = D[hr]
res.append(stats(dr, 'with Gaia RV'))
res.append(stats(dr[dr.W.abs().between(25, 70)], '|W_LSR| 25-70 (target %.0f)' % tW))
res.append(stats(dr[dr.vtot.between(35, 75)], 'v_tot 35-75 (target %.0f)' % out['target']['vtot_LSR']))
res.append(stats(dr[(dr.W.abs() > 20) & (dr.vtot > 35)], '|W|>20 & vtot>35'))
out['groups'] = res
json.dump(out, open(f'{B}/out/r05_kinematics.json', 'w'), indent=1)
