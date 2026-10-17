"""f02a: Gaia DR3 query (read-only ADQL) for every source within 25" of the M3V, joined with
gaiadr3.synthetic_photometry_gspc (XP synthetic photometry: JKC U B V R I, SDSS u g r i z, PS1 y),
plus the Mira's gaiadr3.vari_long_period_variable and vari_summary rows.
Outputs: data/gaia_25arcsec_gspc.csv, data/gaia_mira_lpv.csv, data/gaia_mira_vari_summary.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
q = f"""SELECT g.source_id, g.ra, g.dec, g.ref_epoch, g.pmra, g.pmdec, g.parallax, g.parallax_error, g.ruwe,
g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp, g.has_xp_continuous, g.has_xp_sampled,
g.phot_variable_flag, g.ipd_frac_multi_peak, g.astrometric_params_solved, g.phot_bp_rp_excess_factor,
s.u_jkc_mag, s.b_jkc_mag, s.v_jkc_mag, s.r_jkc_mag, s.i_jkc_mag, s.u_sdss_mag, s.g_sdss_mag, s.r_sdss_mag, s.i_sdss_mag,
s.z_sdss_mag, s.y_ps1_mag,
s.v_jkc_flag, s.i_jkc_flag, s.g_sdss_flag, s.z_sdss_flag,
DISTANCE(POINT('ICRS', g.ra, g.dec), POINT('ICRS', {TARGET['ra']}, {TARGET['dec']}))*3600 AS sep_arcsec
FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN gaiadr3.synthetic_photometry_gspc AS s ON g.source_id = s.source_id
WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {TARGET['ra']}, {TARGET['dec']}, {25/3600.}))"""
r = Gaia.launch_job(q).get_results().to_pandas()
r = r.sort_values('sep_arcsec')
r.to_csv(f'{F}/data/gaia_25arcsec_gspc.csv', index=False)
print(r[['source_id', 'sep_arcsec', 'phot_g_mean_mag', 'bp_rp', 'v_jkc_mag', 'i_jkc_mag', 'g_sdss_mag', 'z_sdss_mag']].to_string())
mira = 5925373268528684288
l = Gaia.launch_job(f"SELECT * FROM gaiadr3.vari_long_period_variable WHERE source_id = {mira}").get_results().to_pandas()
l.to_csv(f'{F}/data/gaia_mira_lpv.csv', index=False); print(l.T)
v = Gaia.launch_job(f"SELECT * FROM gaiadr3.vari_summary WHERE source_id = {mira}").get_results().to_pandas()
v.to_csv(f'{F}/data/gaia_mira_vari_summary.csv', index=False)
print(v[[c for c in v.columns if any(s in c for s in ('min_mag', 'max_mag', 'mean_mag', 'std_dev', 'num_selected'))]].T)
