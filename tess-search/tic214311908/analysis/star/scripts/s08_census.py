"""Neighbour census for blend scenarios around TIC 214311908 (M3V) and the Mira TIC 214311906.
Every Gaia DR3 source within 30" of the M3V or of the Mira, plus DECaPS2 sources within 30" of the M3V that have no
Gaia counterpart within 1.0" (positions compared at the DECaPS2 epoch, ~2018.7).
For each source: G, BP-RP, parallax, PM, TESS magnitude, TESS flux (15000 e-/s at T = 10), the eclipse depth it would
need to remove the observed 15.4 +- 0.8 e-/s (deblended M3V flux loss, phot/ README), and the maximum depth it could
produce (100% eclipse) in apertures of radius 1.0 and 2.0 px centred on the M3V, for the S12 (cam1 ccd2) and S104
(cam1 ccd4) PRFs, with all catalogued stars (including the Mira at its TIC mean T) contributing light.
T magnitudes: Stassun et al. (2019) Gaia colour relation; T = G - 0.43 without colour; DECaPS-only sources from a
fit T(i, i-z) to Gaia-matched DECaPS2 stars in the field (scatter quoted). For the M3V the TIC value 13.0135 is used.
Outputs: out/neighbour_census.csv (all sources), out/census_summary.json"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import PRF
B = 'tic214311908/analysis/star'
PRFDIR = '<local>/tic713/localise/prf_files/s0004/'
M3V, MIRA, PAIRA, PAIRB = 5925373272850784768, 5925373268528684288, 5925373272850782208, 5925373272827253632
DF_OBS, DF_ERR = 15.4, 0.8
g = pd.read_csv(f'{B}/data/gaia_120.csv')
r0 = g[g.source_id == M3V].iloc[0]; rm = g[g.source_id == MIRA].iloc[0]
cosd = np.cos(np.radians(r0.dec))
def offs(ra, dec, ref): return (ra - ref.ra) * cosd * 3600, (dec - ref.dec) * 3600
g['dE'], g['dN'] = offs(g.ra, g.dec, r0)
g['sep_mira'] = np.hypot(*offs(g.ra, g.dec, rm))
def tmag(G, c):
    c = np.asarray(c, float); G = np.asarray(G, float)
    T = G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), T, G - 0.430)
g['T'] = tmag(g.phot_g_mean_mag, g.bp_rp)
g['T_src'] = np.where(np.isfinite(g.bp_rp), 'Stassun19 Gaia colour', 'G-0.43')
g.loc[g.source_id == M3V, 'T'] = 13.0135; g.loc[g.source_id == M3V, 'T_src'] = 'TIC v8.2'
g.loc[g.source_id == MIRA, 'T'] = 11.2566; g.loc[g.source_id == MIRA, 'T_src'] = 'TIC v8.2 (variable, x2)'
# pair B (2-parameter solution, no colour): an M4-M4.5 dwarf at the pair's distance (M_G = 11.67) -> assume BP-RP = 3.06
g.loc[g.source_id == PAIRB, 'T'] = float(tmag(15.740896, 3.06)); g.loc[g.source_id == PAIRB, 'T_src'] = 'Stassun19 with assumed BP-RP 3.06 (M4V companion)'
g['catalog'] = 'GaiaDR3'
sel = g[(g.sep_m3v <= 30) | (g.sep_mira <= 30)].copy()
# --- DECaPS2-only sources
d = pd.read_csv(f'{B}/data/decaps2_90.csv', usecols=['ra', 'dec', 'epochmean', 'ndet', 'ndet_ok', 'mean_mag_g', 'mean_mag_r', 'mean_mag_i',
                                                     'mean_mag_z', 'mean_mag_y', 'nmag_ok_i', 'nmag_ok_z', 'fracflux_avg_i', 'posstdev'])
ep = 2000.0 + (d.epochmean.median() + 2400000.5 - 2451545.0) / 365.25
dt = ep - 2016.0
g['ra_ep'] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / cosd; g['de_ep'] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
from scipy.spatial import cKDTree
tg = cKDTree(np.column_stack([g.ra_ep * cosd * 3600, g.de_ep * 3600]))
dd, ii = tg.query(np.column_stack([d.ra * cosd * 3600, d.dec * 3600]))
d['gaia_match_dist'] = dd; d['gaia_id'] = g.source_id.values[ii]
# calibrate T from DECaPS i, z on Gaia-matched stars with colours (field within 90")
m = (dd < 0.5) & np.isfinite(d.mean_mag_i) & np.isfinite(d.mean_mag_z) & (d.nmag_ok_i > 0) & (d.nmag_ok_z > 0)
cal = d[m].copy(); cal['T'] = g.set_index('source_id').loc[cal.gaia_id, 'T'].values
cal['bprp'] = g.set_index('source_id').loc[cal.gaia_id, 'bp_rp'].values
cal = cal[np.isfinite(cal.bprp) & (cal['T'] < 19) & (cal.mean_mag_i > 14.5) & np.isfinite(cal['T'])]
A = np.column_stack([np.ones(len(cal)), cal.mean_mag_i, cal.mean_mag_i - cal.mean_mag_z])
coef, *_ = np.linalg.lstsq(A, cal['T'], rcond=None)
rs = cal['T'] - A @ coef; sc = 1.4826 * np.median(np.abs(rs - np.median(rs)))
print(f'T = {coef[0]:.3f} + {coef[1]:.3f} i + {coef[2]:.3f} (i-z); n = {len(cal)}, robust scatter {sc:.3f} mag')
cosd0 = cosd
d['dE'] = (d.ra - (r0.ra + r0.pmra * dt / 3.6e6 / cosd)) * cosd * 3600; d['dN'] = (d.dec - (r0.dec + r0.pmdec * dt / 3.6e6)) * 3600
d['sep_m3v'] = np.hypot(d.dE, d.dN)
dn = d[(d.gaia_match_dist > 1.0) & (d.sep_m3v <= 30)].copy()
iz = np.where(np.isfinite(dn.mean_mag_z), dn.mean_mag_i - dn.mean_mag_z, 0.3)
dn['T'] = np.where(np.isfinite(dn.mean_mag_i), coef[0] + coef[1] * dn.mean_mag_i + coef[2] * iz,
                   np.where(np.isfinite(dn.mean_mag_z), dn.mean_mag_z - 0.2, np.where(np.isfinite(dn.mean_mag_y), dn.mean_mag_y - 0.2, np.nan)))
dn['T_src'] = np.where(np.isfinite(dn.mean_mag_i), 'DECaPS2 i,z (cal)', 'DECaPS2 z or y - 0.2')
dn['decaps_reliable'] = (dn.ndet_ok >= 2) & np.isfinite(dn['T'])
dn['catalog'] = 'DECaPS2-only'
dn['sep_mira'] = np.hypot((dn.ra - rm.ra) * cosd * 3600, (dn.dec - rm.dec) * 3600)
dn = dn.rename(columns={'mean_mag_i': 'decaps_i', 'mean_mag_z': 'decaps_z', 'mean_mag_r': 'decaps_r'})
cols = ['catalog', 'source_id', 'ra', 'dec', 'dE', 'dN', 'sep_m3v', 'sep_mira', 'phot_g_mean_mag', 'bp_rp', 'parallax', 'parallax_error',
        'pmra', 'pmdec', 'ruwe', 'astrometric_params_solved', 'T', 'T_src']
sel = sel.assign(source_id=sel.source_id.astype('int64').astype(str))   # keep 64-bit ids exact (no float conversion)
allc = pd.concat([sel[cols], dn[['catalog', 'ra', 'dec', 'dE', 'dN', 'sep_m3v', 'sep_mira', 'T', 'T_src', 'decaps_r', 'decaps_i', 'decaps_z', 'ndet', 'ndet_ok', 'decaps_reliable']]],
                 ignore_index=True)
allc['F_tess'] = 15000 * 10 ** (-0.4 * (allc['T'] - 10))
allc['need_ecl_depth'] = DF_OBS / allc.F_tess
allc['need_ecl_depth_lo'] = (DF_OBS - 2 * DF_ERR) / allc.F_tess
# --- aperture dilution with the PRF (all Gaia sources within 120" + DECaPS-only within 30" contribute light)
light = pd.concat([g[['dE', 'dN', 'T']], dn[['dE', 'dN', 'T']]], ignore_index=True)
light = light[np.isfinite(light['T'])]
light['F'] = 15000 * 10 ** (-0.4 * (light['T'] - 10))
for lab, (cam, ccd, sec, col, row) in {'s12': (1, 2, 12, 1550 + 10, 260 + 10), 's104': (1, 4, 104, 1310 + 10, 1465 + 10)}.items():
    p = PRF.TESS_PRF(cam, ccd, sec, col, row, localdatadir=PRFDIR)
    # pixel scale ~21"; orientation: ignore rotation (aperture is circular, PRF nearly symmetric at this level)
    pix = 21.0
    N = 41; c0 = 20.0
    def img(dx, dy):
        return p.locate(c0 + dx / pix, c0 + dy / pix, (N, N))
    yy, xx = np.mgrid[0:N, 0:N]
    tot = {1.0: 0.0, 2.0: 0.0}; fr = {}
    for i, r in light.iterrows():
        if abs(r.dE) > 200 or abs(r.dN) > 200: continue
        im = img(r.dE, r.dN)
        for rad in (1.0, 2.0):
            ap = np.hypot(xx - c0, yy - c0) <= rad
            f = float(im[ap].sum()); fr[(i, rad)] = f; tot[rad] += f * r.F
    for rad in (1.0, 2.0):
        fracs = []
        for _, r in allc.iterrows():
            im = img(r.dE, r.dN); ap = np.hypot(xx - c0, yy - c0) <= rad
            fracs.append(float(im[ap].sum()))
        allc[f'apfrac_r{rad:.0f}_{lab}'] = fracs
        allc[f'maxdepth_r{rad:.0f}_{lab}'] = allc.F_tess * allc[f'apfrac_r{rad:.0f}_{lab}'] / tot[rad]
    m3 = allc.source_id == str(M3V)
    print(lab, 'total aperture flux r1/r2:', {k: round(v, 1) for k, v in tot.items()},
          ' M3V share r1:', float(allc.loc[m3, f'maxdepth_r1_{lab}'].iloc[0]), ' r2:', float(allc.loc[m3, f'maxdepth_r2_{lab}'].iloc[0]))
    # observed depth in that aperture if the M3V is the host: DF_OBS*frac_M3V/tot
    obs_r1 = DF_OBS * float(allc.loc[m3, f'apfrac_r1_{lab}'].iloc[0]) / tot[1.0]
    allc.attrs[f'obs_depth_r1_{lab}'] = obs_r1
allc['can_supply_flux'] = allc.need_ecl_depth_lo <= 1.0
names = {M3V: 'TIC 214311908 (M3V, candidate host)', MIRA: 'TIC 214311906 (Mira, distant AGB star)', PAIRA: 'TIC 214311910 A (bound M3-M4 dwarf, 757 au)',
         PAIRB: 'Gaia companion of TIC 214311910 A (0.83", bound M4 dwarf; 2-parameter solution)', 5925373272827260544: 'background star 3.1" NE (1.5 kpc)'}
allc['note'] = allc.source_id.map({str(k): v for k, v in names.items()}).fillna('')
allc['plx_over_err'] = allc.parallax / allc.parallax_error
allc['distance_class'] = np.where(allc.plx_over_err > 5, np.where(allc.parallax > 10, 'nearby (<100 pc)', 'background (>100 pc)'), 'parallax not significant')
allc = allc.sort_values('sep_m3v').reset_index(drop=True)
allc.to_csv(f'{B}/out/neighbour_census.csv', index=False)
pd.set_option('display.width', 250)
show = allc[(allc.sep_m3v < 30)][['catalog', 'source_id', 'dE', 'dN', 'sep_m3v', 'phot_g_mean_mag', 'bp_rp', 'parallax', 'T', 'F_tess', 'need_ecl_depth',
                                   'maxdepth_r1_s104', 'maxdepth_r2_s104', 'maxdepth_r1_s12']]
print(show.round(4).to_string())
summ = dict(T_decaps_cal=dict(coef=coef.tolist(), scatter=sc, n=len(cal)), dF_obs=[DF_OBS, DF_ERR],
            n_gaia_30_m3v=int(((allc.catalog == 'GaiaDR3') & (allc.sep_m3v <= 30)).sum()),
            n_gaia_30_union=int((allc.catalog == 'GaiaDR3').sum()),
            n_decaps_only_30=int((allc.catalog == 'DECaPS2-only').sum()),
            capable=allc[allc.can_supply_flux & (allc.source_id != str(M3V))][['catalog', 'source_id', 'sep_m3v', 'sep_mira', 'T', 'need_ecl_depth', 'parallax']].round(3).to_dict('records'))
json.dump(summ, open(f'{B}/out/census_summary.json', 'w'), indent=1, default=str)
print(json.dumps(summ['capable'], indent=0, default=str))
