"""s07: Gaia DR3 parameters of the two catalogued background candidates (3.1" bg3 = 5925373272827260544,
5.5" bg5 = 5925373272827262208) and of the M3V: astrometry, photometry, GSP-Phot (Teff, logg, A_G, E(BP-RP),
distance), DSC class probabilities (incl. white dwarf and binary), XP availability; plus the SFD/Planck
extinction toward the field (IRSA dust service) for the luminosity/colour checks of the EB scenarios.
Output: data/gaia_bg.csv, out/gaia_bg.json"""
import json, numpy as np
from astroquery.gaia import Gaia
ids = '5925373272827260544, 5925373272827262208, 5925373272850784768'
q = f"""SELECT g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.ruwe, g.phot_g_mean_mag, g.phot_bp_mean_mag,
g.phot_rp_mean_mag, g.bp_rp, g.phot_bp_rp_excess_factor, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.has_xp_continuous, g.l, g.b,
g.classprob_dsc_combmod_star, a.classprob_dsc_combmod_whitedwarf, a.classprob_dsc_combmod_binarystar, g.classprob_dsc_combmod_quasar, g.classprob_dsc_combmod_galaxy,
a.teff_gspphot, a.teff_gspphot_lower, a.teff_gspphot_upper, a.logg_gspphot, a.mh_gspphot, a.distance_gspphot, a.distance_gspphot_lower, a.distance_gspphot_upper,
a.ag_gspphot, a.ebpminrp_gspphot, a.radius_gspphot, a.classprob_dsc_allosmod_star, a.classprob_dsc_specmod_whitedwarf, a.classprob_dsc_specmod_binarystar,
a.libname_gspphot, a.teff_esphs, a.spectraltype_esphs, a.flags_esphs, a.teff_msc1, a.teff_msc2, a.logg_msc1, a.logg_msc2, a.distance_msc, a.ag_msc
FROM gaiadr3.gaia_source AS g LEFT JOIN gaiadr3.astrophysical_parameters AS a ON g.source_id = a.source_id
WHERE g.source_id IN ({ids})"""
t = Gaia.launch_job(q).get_results().to_pandas()
t.to_csv('data/gaia_bg.csv', index=False)
pd_ = t.set_index('source_id').T
print(pd_.to_string())
out = {str(r.source_id): {k: (None if (isinstance(v, float) and not np.isfinite(v)) else (v if not hasattr(v, 'item') else v.item())) for k, v in r.items()} for _, r in t.iterrows()}
# extinction
try:
    from astroquery.ipac.irsa.irsa_dust import IrsaDust
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord(260.9697, -51.3238, unit='deg')
    tab = IrsaDust.get_query_table(c, section='ebv')
    out['irsa_dust'] = {k: float(tab[k][0]) for k in tab.colnames if 'ebv' in k.lower() or 'mean' in k.lower()}
    print(out['irsa_dust'])
except Exception as ex:
    out['irsa_dust'] = str(ex); print('dust failed', ex)
json.dump(out, open('out/gaia_bg.json', 'w'), indent=1, default=str)
