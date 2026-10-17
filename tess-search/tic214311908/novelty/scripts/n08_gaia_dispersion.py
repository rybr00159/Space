"""n08: context for the Maiz Apellaniz et al. (2023, J/A+A/677/A137) Gaia DR3 photometric-dispersion entry of the M3V
(flag VVV): (1) the M3V's Gaia DR3 blending/contamination counters; (2) the distribution of S_G and of the VarF flag for
comparable nearby M dwarfs (same G and BP-RP range, parallax > 8 mas) within 20 deg, via VizieR TAP. Read-only.
Output: out/gaia_dispersion_context.json
Usage: python -I n08_gaia_dispersion.py <novelty_dir>
"""
import sys, os, io, json, time, datetime, urllib.parse
import numpy as np, pandas as pd, requests
from astroquery.gaia import Gaia

D = sys.argv[1]
OUT, RAW = os.path.join(D, "out"), os.path.join(D, "raw")
res = {"queried_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
q = ("SELECT source_id, phot_g_n_obs, phot_bp_n_obs, phot_rp_n_obs, phot_bp_n_blended_transits, phot_rp_n_blended_transits, "
     "phot_bp_n_contaminated_transits, phot_rp_n_contaminated_transits, phot_bp_rp_excess_factor, phot_g_mean_flux_over_error, "
     "phot_bp_mean_flux_over_error, phot_rp_mean_flux_over_error FROM gaiadr3.gaia_source WHERE source_id = 5925373272850784768")
for k in range(3):
    try:
        g = Gaia.launch_job(q).get_results().to_pandas(); break
    except Exception as e:
        g = pd.DataFrame({"error": [str(e)[:200]]}); time.sleep(10)
res["m3v_gaia_photometry_counters"] = g.astype(str).to_dict("records")
print(g.T.to_string())

adql = ('SELECT "Source", "Gmag", "BP-RP", "plx", "SG", "SSG", "SBP", "SRP", "VarF" FROM "J/A+A/677/A137/catalog" '
        "WHERE 1 = CONTAINS(POINT('ICRS', RA_ICRS, DE_ICRS), CIRCLE('ICRS', 260.97, -51.32, 20)) "
        'AND "Gmag" BETWEEN 13.8 AND 14.8 AND "BP-RP" BETWEEN 2.6 AND 2.95 AND "plx" > 8')
# TAPVizieR returned HTTP 503 on 2026-10-17, so the same selection is made with the VizieR cone service (15 deg radius)
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
try:
    v = Vizier(columns=["Source", "Gmag", "BP-RP", "plx", "SG", "SSG", "SBP", "SRP", "VarF"], row_limit=-1, timeout=900,
               column_filters={"Gmag": "13.8..14.8", "BP-RP": "2.6..2.95", "plx": ">8"})
    t = v.query_region(SkyCoord(260.97, -51.32, unit="deg"), radius=15 * u.deg, catalog="J/A+A/677/A137/catalog")[0].to_pandas()
    t.to_csv(os.path.join(RAW, "maiz23_comparison_mdwarfs.csv"), index=False)
    t["VarG"] = t["VarF"].astype(str).str[0]
    m3v = 8.000615
    res["comparison"] = {"selection": "VizieR cone 15 deg around the M3V; 13.8 < G < 14.8, 2.6 < BP-RP < 2.95, plx > 8 mas",
                         "n": int(len(t)), "frac_G_variable": float((t.VarG == "V").mean()),
                         "frac_all_three_V": float((t["VarF"].astype(str) == "VVV").mean()),
                         "SG_percentiles_mmag": {p: float(np.nanpercentile(t.SG, p)) for p in (16, 50, 84, 95)},
                         "m3v_SG_percentile": float((t.SG < m3v).mean() * 100)}
except Exception as e:
    res["comparison"] = {"error": str(e)[:300]}
print(json.dumps(res["comparison"], indent=1))
json.dump(res, open(os.path.join(OUT, "gaia_dispersion_context.json"), "w"), indent=1)
