"""Identifier-level checks that were first run interactively, now reproducible:
  - ExoFOP target JSON for the three TICs and the pair's extra TIC entries -> ext/exofop/target_<TIC>.json
  - TIC 8.2 sources within 30" of the M3V -> out/tic_30arcsec.ecsv
  - Gaia DR3 flags (RUWE, IPD, epoch photometry, RVs) for the stars of interest -> out/gaia_dr3_flags.json
  - SIMBAD bibliography of every SIMBAD object within 20" of the Mira -> out/simbad_bib.json
  - arXiv API search for each identifier (title/abstract/comments) -> out/arxiv_id_search.txt
"""
import json, time, urllib.parse, warnings
warnings.filterwarnings("ignore")
import requests
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

for tic in ("214311906", "214311908", "214311910", "1515709322", "1515709367"):
    r = requests.get(f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json", timeout=120)
    if r.ok:
        open(f"ext/exofop/target_{tic}.json", "w").write(r.text)
tic = Catalogs.query_region(SkyCoord(260.9696737, -51.3238212, unit="deg"), radius=30 * u.arcsec, catalog="TIC")
tic.write("out/tic_30arcsec.ecsv", overwrite=True)
ids = (5925373268528684288, 5925373272827262208, 5925373272827253632, 5925373272850784768, 5925373272850782208, 5925373272827260544)
q = ("SELECT source_id, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe, phot_variable_flag, has_epoch_photometry, has_rvs, "
     "non_single_star, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, phot_bp_rp_excess_factor, phot_bp_n_blended_transits, "
     "phot_rp_n_blended_transits, phot_g_n_obs, phot_bp_n_obs, phot_g_mean_flux_over_error, matched_transits, visibility_periods_used, "
     "radial_velocity, radial_velocity_error, rv_nb_transits, teff_gspphot, classprob_dsc_combmod_star "
     f"FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(str(i) for i in ids)})")
g = Gaia.launch_job(q).get_results().to_pandas().fillna(0)
json.dump(g.to_dict("records"), open("out/gaia_dr3_flags.json", "w"), indent=1, default=str)
# SIMBAD bibliography (TAP; explicit join on has_ref.oidref)
bib = {}
r = Simbad.query_tap("SELECT b.main_id, b.ra, b.dec FROM basic AS b WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', 260.97765, -51.32369, 0.00556)) = 1")
for row in r:
    mid = str(row["main_id"])
    q = ("SELECT r.bibcode, r.title FROM ref AS r JOIN has_ref AS h ON h.oidbibref = r.oidbib JOIN basic AS b ON b.oid = h.oidref "
         f"WHERE b.main_id = '{mid}'")
    try:
        t = Simbad.query_tap(q); bib[mid] = [f"{x['bibcode']}: {x['title']}" for x in t]
    except Exception as e:
        bib[mid] = f"error {e}"
json.dump(bib, open("out/simbad_bib.json", "w"), indent=1)
lines = []
for ident in ("214311908", "214311906", "214311910", "5925373272850784768", "5925373268528684288", "5925373272850782208",
              "J17235464-5119254", "J17235292-5119271"):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": f'all:"{ident}"', "max_results": 5})
    txt = requests.get(url, timeout=60).text
    n = txt.count("<entry>"); lines.append(f'arXiv all:"{ident}" -> {n} entries'); time.sleep(3.5)
open("out/arxiv_id_search.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines)); print(json.dumps(bib, indent=1)[:1500]); print(len(tic), "TIC rows within 30 arcsec")
