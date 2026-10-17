"""Referee: independent Gaia DR3/DR2 queries.
- all DR3 sources within 45" of the M3V (TIC 214311908) with photometric flags
- DR2 counterparts (dr2_neighbourhood) of the M3V, pair A, pair B, Mira
- DR3 astrophysical_parameters for the M3V
Outputs: data/gaia_cone45.csv, data/gaia_dr2_xm.csv, data/gaia_ap_target.csv"""
import warnings, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
D = 'tic214311908/analysis/referee_star/data'
RA, DEC = 260.969674, -51.323821
def run(q):
    for k in range(3):
        try: return Gaia.launch_job_async(q).get_results().to_pandas()
        except Exception as e: print('retry', e); time.sleep(15)
q = f"""select g.*, distance(point({RA},{DEC}), point(g.ra,g.dec))*3600 as sep
from gaiadr3.gaia_source as g where 1=contains(point(g.ra,g.dec), circle({RA},{DEC},{45/3600}))"""
t = run(q); t.to_csv(f'{D}/gaia_cone45.csv', index=False); print('cone', len(t))
ids = '5925373272850784768,5925373272850782208,5925373272827253632,5925373268528684288'
q = f"""select n.*, d.ra as ra2, d.dec as dec2, d.pmra as pmra2, d.pmdec as pmdec2, d.pmra_error as pmra2_err, d.pmdec_error as pmdec2_err,
d.parallax as plx2, d.parallax_error as plx2_err, d.phot_g_mean_mag as g2, d.astrometric_excess_noise as aen2, d.visibility_periods_used as vpu2
from gaiadr3.dr2_neighbourhood as n join gaiadr2.gaia_source as d on n.dr2_source_id=d.source_id where n.dr3_source_id in ({ids})"""
d = run(q); d.to_csv(f'{D}/gaia_dr2_xm.csv', index=False); print(d.to_string())
q = f"select * from gaiadr3.astrophysical_parameters where source_id=5925373272850784768"
a = run(q); a.to_csv(f'{D}/gaia_ap_target.csv', index=False)
