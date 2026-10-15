import warnings, sys, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
Gaia.ROW_LIMIT=-1
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/'
q="""SELECT source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,parallax,parallax_over_error,ruwe,
DISTANCE(POINT(ra,dec),POINT(330.887948,-42.835511))*60 AS sep_arcmin
FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE(330.887948,-42.835511,0.5)) AND phot_g_mean_mag<19.5"""
for a in range(3):
    try:
        t=Gaia.launch_job(q).get_results(); t.write(D+'gaia_counts_30arcmin.ecsv',overwrite=True); print('gaia',len(t)); break
    except Exception as e: print('retry',str(e)[:100]); time.sleep(5)
Vizier.TIMEOUT=200
v=Vizier(columns=['KIC','Per','Morph','Kpmag','Teff','GLAT'],row_limit=-1)
r=v.get_catalogs('J/AJ/151/68/catalog'); r[0].write(D+'kepler_eb_kirk2016.ecsv',overwrite=True); print('kepler EBs',len(r[0]))
