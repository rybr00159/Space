import warnings; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
for q in ["SELECT * FROM gaiadr3.tycho2tdsc_merge_best_neighbour WHERE source_id=6571343840265900544",
          "SELECT TOP 2 * FROM gaiadr1.tgas_source WHERE source_id=6571343840265900544",
          "SELECT TOP 1 * FROM public.tycho2"]:
    try:
        t=Gaia.launch_job(q).get_results(); print(t.colnames); print(t)
    except Exception as e: print('ERR',str(e)[:150])
