"""Independent Gaia DR3 neighbour queries (referee). Saves CSVs to ../data/."""
import sys, os
import numpy as np
from astroquery.gaia import Gaia

OUT = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/data'
TARGETS = {
    90086266: (276.27770, -42.64403, 150),
    61332759: (287.61974, -44.73357, 200),
    214311906: (260.97765, -51.32369, 150),
    55829034: (272.03163, -35.69140, 100),
    61141412: (286.96474, -43.00969, 150),
}
COLS = ("source_id, ra, dec, ref_epoch, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, "
        "phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe, radial_velocity, radial_velocity_error, "
        "rv_nb_transits, rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof, rv_method_used, "
        "grvs_mag, rv_template_teff, vbroad, phot_variable_flag, non_single_star, ipd_frac_multi_peak, "
        "ipd_gof_harmonic_amplitude, visibility_periods_used, teff_gspphot, classprob_dsc_combmod_star")

for tic, (ra, dec, rad) in TARGETS.items():
    fn = f'{OUT}/gaia_{tic}.csv'
    if os.path.exists(fn):
        continue
    q = (f"SELECT {COLS}, DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec}))*3600. AS sep "
         f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {rad/3600.}))")
    j = Gaia.launch_job_async(q)
    t = j.get_results()
    t.sort('sep')
    t.write(fn, format='csv', overwrite=True)
    print(tic, len(t))
