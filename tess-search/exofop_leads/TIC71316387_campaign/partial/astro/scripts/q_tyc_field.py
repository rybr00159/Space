"""Tycho-2 / TGAS / Gaia DR3 astrometry for field stars around the target (empirical null for proper-motion anomaly)."""
import warnings, sys, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astropy.table import vstack, unique
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
Gaia.ROW_LIMIT=-1
c0=SkyCoord(330.887948,-42.835511,unit='deg')
cen=[(0,0)]+[(r*np.cos(a),r*np.sin(a)) for r in [3.0] for a in np.linspace(0,2*np.pi,6,endpoint=False)]
tabs=[]
for dx,dy in cen:
    c=c0.spherical_offsets_by(dx*u.deg,dy*u.deg)
    q=f"""SELECT g.source_id,g.ra,g.dec,g.pmra,g.pmdec,g.pmra_error,g.pmdec_error,g.ra_error,g.dec_error,g.parallax,g.phot_g_mean_mag,g.bp_rp,g.ruwe,
    g.astrometric_excess_noise,g.ipd_frac_multi_peak,g.radial_velocity,g.radial_velocity_error,g.rv_nb_transits,g.rv_chisq_pvalue,g.rv_amplitude_robust,g.non_single_star,
    t.ra_mdeg,t.de_mdeg,t.pm_ra,t.pm_de,t.e_pm_ra,t.e_pm_de,t.e_ra_mdeg,t.e_de_mdeg,t.ep_ra_m,t.ep_de_m,t.ra_deg,t.de_deg,t.e_ra_deg,t.e_de_deg,t.ep_ra1990,t.ep_de1990,t.vt_mag,t.bt_mag,t.num,t.q_ra_mdeg,t.q_de_mdeg,
    s.pmra AS tgas_pmra,s.pmdec AS tgas_pmdec,s.pmra_error AS tgas_pmra_error,s.pmdec_error AS tgas_pmdec_error,s.pmra_pmdec_corr AS tgas_corr
    FROM gaiadr3.gaia_source AS g JOIN gaiadr3.tycho2tdsc_merge_best_neighbour AS x ON g.source_id=x.source_id
    JOIN public.tycho2 AS t ON x.original_ext_source_id=t.id
    LEFT OUTER JOIN gaiadr1.tgas_source AS s ON s.tycho2_id=t.id
    WHERE 1=CONTAINS(POINT(g.ra,g.dec),CIRCLE({c.ra.deg},{c.dec.deg},1.6)) AND t.hip IS NULL"""
    for a in range(3):
        try:
            t=Gaia.launch_job(q).get_results(); tabs.append(t); print(dx,dy,len(t),flush=True); break
        except Exception as e: print('retry',str(e)[:100],flush=True); time.sleep(5)
T=unique(vstack(tabs),keys='source_id'); T.write(sys.argv[1],overwrite=True); print('total',len(T))
