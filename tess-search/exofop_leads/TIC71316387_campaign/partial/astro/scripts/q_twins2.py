import sys, warnings, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astropy.table import vstack
Gaia.ROW_LIMIT=-1
cols="""source_id,ra,dec,l,b,parallax,parallax_error,pmra,pmdec,ruwe,astrometric_excess_noise,astrometric_excess_noise_sig,
ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,phot_bp_rp_excess_factor,
phot_bp_n_blended_transits,phot_bp_n_obs,radial_velocity,radial_velocity_error,rv_nb_transits,rv_chisq_pvalue,rv_amplitude_robust,
teff_gspphot,logg_gspphot,mh_gspphot,ebpminrp_gspphot,non_single_star,has_xp_sampled,astrometric_n_good_obs_al,visibility_periods_used,grvs_mag,rv_expected_sig_to_noise"""
out=sys.argv[1]; R=float(sys.argv[2])
tabs=[]
# tile in small circles to keep sync queries short
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
c0=SkyCoord(330.887948,-42.835511,unit='deg')
centres=[(0,0)]+[(r*np.cos(a),r*np.sin(a)) for r in [2.6,5.2] for a in np.linspace(0,2*np.pi,int(6*r/2.6),endpoint=False)]
for dx,dy in centres:
    c=c0.spherical_offsets_by(dx*u.deg,dy*u.deg)
    q=(f"SELECT {cols} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({c.ra.deg},{c.dec.deg},{R})) "
       "AND bp_rp BETWEEN 0.66 AND 0.92 AND phot_g_mean_mag BETWEEN 10.5 AND 13.0 AND parallax_over_error>30")
    for attempt in range(3):
        try:
            t=Gaia.launch_job(q).get_results(); tabs.append(t); print(dx,dy,len(t),flush=True); break
        except Exception as e:
            print('retry',str(e)[:80]); time.sleep(5)
T=vstack(tabs)
from astropy.table import unique
T=unique(T,keys='source_id'); T.write(out,overwrite=True); print('total',len(T))
