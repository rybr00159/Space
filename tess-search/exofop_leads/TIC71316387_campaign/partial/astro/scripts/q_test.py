import warnings; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
qs=["SELECT * FROM gaiadr3.tmass_psc_xsc_best_neighbour WHERE source_id=6571343840265900544",
"SELECT * FROM gaiadr3.allwise_best_neighbour WHERE source_id=6571343840265900544",
"SELECT TOP 1 * FROM gaiadr1.tmass_original_valid",
"SELECT TOP 1 * FROM gaiadr1.allwise_original_valid",
"SELECT TOP 1 * FROM gaiadr3.tmass_psc_xsc_join"]
for q in qs:
    try:
        t=Gaia.launch_job(q).get_results(); print(q[:60], len(t), t.colnames[:25]); print(t[:1])
    except Exception as e: print(q[:60],'ERR',str(e)[:200])
