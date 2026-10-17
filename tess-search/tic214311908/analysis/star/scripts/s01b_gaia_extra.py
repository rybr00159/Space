"""Gaia DR3 auxiliary tables for the key stars (NSS, variability, DR2 astrometry, 2MASS/AllWISE best neighbours).
Run after s01_gaia.py (which writes data/gaia_120.csv and data/gaia_ap.csv). Uses synchronous TAP queries.
Outputs: data/gaia_nss_vari.json, data/gaia_dr2.csv, data/gaia_xm.csv"""
import json, warnings, sys
warnings.filterwarnings('ignore')
import pandas as pd
from astroquery.gaia import Gaia
D = 'tic214311908/analysis/star/data'
ids = [5925373272850784768, 5925373272850782208, 5925373272827253632, 5925373268528684288, 5925373272827260544, 5925373204131295872]
idl = ','.join(str(i) for i in ids)
def q(adql):
    return Gaia.launch_job(adql).get_results().to_pandas()
out = {}
for tab in ['nss_two_body_orbit', 'nss_acceleration_astro', 'nss_non_linear_spectro', 'vari_classifier_result',
            'vari_summary', 'vari_long_period_variable']:
    try:
        t = q(f"select * from gaiadr3.{tab} where source_id in ({idl})")
        out[tab] = t.astype(object).where(pd.notnull(t), None).to_dict('records')
        print(tab, len(t), list(t.source_id) if len(t) else '', flush=True)
    except Exception as e:
        out[tab] = f'error {e}'; print(tab, 'error', e, flush=True)
json.dump(out, open(f'{D}/gaia_nss_vari.json', 'w'), indent=1, default=str)
dr2 = q(f"select n.dr3_source_id, n.dr2_source_id, n.angular_distance, n.magnitude_difference, d.ra, d.dec, d.parallax, "
        f"d.parallax_error, d.pmra, d.pmra_error, d.pmdec, d.pmdec_error, d.astrometric_excess_noise, d.phot_g_mean_mag, "
        f"d.bp_rp, d.radial_velocity, d.radial_velocity_error, d.rv_nb_transits, d.astrometric_params_solved "
        f"from gaiadr3.dr2_neighbourhood as n join gaiadr2.gaia_source as d on n.dr2_source_id = d.source_id "
        f"where n.dr3_source_id in ({idl})")
dr2.to_csv(f'{D}/gaia_dr2.csv', index=False); print(dr2.T.to_string(), flush=True)
xm = []
for t in ['tmass_psc_xsc_best_neighbour', 'allwise_best_neighbour']:
    x = q(f"select * from gaiadr3.{t} where source_id in ({idl})"); x['table'] = t; xm.append(x)
pd.concat(xm).to_csv(f'{D}/gaia_xm.csv', index=False); print(pd.concat(xm).to_string())
