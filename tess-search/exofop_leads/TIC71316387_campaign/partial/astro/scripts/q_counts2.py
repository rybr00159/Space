import warnings, sys, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astropy.table import vstack
Gaia.ROW_LIMIT=-1
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/'
tabs=[]
for g1,g2 in [(0,15),(15,16),(16,17),(17,17.6),(17.6,18.2),(18.2,18.8),(18.8,19.3),(19.3,19.8),(19.8,20.3)]:
    q=f"""SELECT source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,parallax,parallax_over_error,ruwe,
    DISTANCE(POINT(ra,dec),POINT(330.887948,-42.835511))*60 AS sep_arcmin
    FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE(330.887948,-42.835511,0.5)) AND phot_g_mean_mag>={g1} AND phot_g_mean_mag<{g2}"""
    for a in range(3):
        try:
            t=Gaia.launch_job(q).get_results(); tabs.append(t); print(g1,g2,len(t),flush=True); break
        except Exception as e: print('retry',str(e)[:100]); time.sleep(5)
T=vstack(tabs); T.write(D+'gaia_counts_30arcmin.ecsv',overwrite=True); print('total',len(T))
