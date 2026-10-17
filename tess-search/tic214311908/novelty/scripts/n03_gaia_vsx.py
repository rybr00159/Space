"""n03: variability and binarity flags.
  - Gaia DR3 (ESA archive TAP): gaia_source flags, all vari_* tables, all nss_* tables for the M3V, the pair, the Mira and
    the 3-6" neighbours; any vari/NSS entry within 60" of the M3V.
  - VSX: VizieR B/vsx/vsx cone (5 arcmin) and the AAVSO VSX API (5 arcmin).
Read-only queries. Outputs: out/gaia_vari_nss.json, out/vsx_status.json, raw/gaia_*.csv, raw/vsx_*.{csv,json}
Usage: python -I n03_gaia_vsx.py <novelty_dir>
"""
import sys, os, io, json, time, datetime
import numpy as np, pandas as pd, requests
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
RA16, DEC16 = 260.96967369918, -51.3238212118
SRC = {"M3V": 5925373272850784768, "pairA": 5925373272850782208, "pairB": 5925373272827253632,
       "Mira": 5925373268528684288, "star_3p1": 5925373272827260544, "star_4p3": 5925373272827258752,
       "star_5p8": 5925373272827262208}
ids = ",".join(str(v) for v in SRC.values())
G = {"queried_utc": NOW}


def gq(q, tag):
    for k in range(3):
        try:
            t = Gaia.launch_job(q).get_results().to_pandas()
            t.to_csv(os.path.join(RAW, f"gaia_{tag}.csv"), index=False)
            return t
        except Exception as e:
            err = e; time.sleep(10)
    return pd.DataFrame({"error": [str(err)[:300]]})


gs = gq("SELECT source_id, phot_g_mean_mag, bp_rp, parallax, ruwe, phot_variable_flag, non_single_star, has_epoch_photometry, "
        "has_epoch_rv, has_rvs, radial_velocity, radial_velocity_error, rv_nb_transits, rv_chisq_pvalue, rv_renormalised_gof, "
        "rv_amplitude_robust, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, duplicated_source, classprob_dsc_combmod_star "
        f"FROM gaiadr3.gaia_source WHERE source_id IN ({ids})", "source_flags")
lab = {v: k for k, v in SRC.items()}
gs["label"] = gs.source_id.map(lab)
G["gaia_source_flags"] = gs.astype(str).to_dict("records")
tabs = ["vari_summary", "vari_classifier_result", "vari_rotation_modulation", "vari_eclipsing_binary", "vari_planetary_transit",
        "vari_short_timescale", "vari_long_period_variable", "vari_rrlyrae", "vari_cepheid", "vari_agn", "vari_compact_companion",
        "vari_ms_oscillator", "vari_microlensing",
        "nss_two_body_orbit", "nss_acceleration_astro", "nss_non_linear_spectro", "nss_vim_fl"]
G["tables"] = {}
for t in tabs:
    r = gq(f"SELECT * FROM gaiadr3.{t} WHERE source_id IN ({ids})", t)
    if "error" in r.columns:
        G["tables"][t] = {"error": r["error"][0]}
    else:
        G["tables"][t] = {"n": int(len(r)), "sources": [lab.get(int(s), str(s)) for s in r.source_id] if len(r) else [],
                          "rows": r.head(5).astype(str).to_dict("records")}
    print(t, G["tables"][t].get("n", G["tables"][t].get("error")), G["tables"][t].get("sources"))
# any vari_summary / NSS entry within 60" of the M3V
near = gq("SELECT g.source_id, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.non_single_star, g.phot_variable_flag, "
          f"DISTANCE(POINT('ICRS', g.ra, g.dec), POINT('ICRS', {RA16}, {DEC16}))*3600 AS sep "
          "FROM gaiadr3.gaia_source AS g WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec), "
          f"CIRCLE('ICRS', {RA16}, {DEC16}, {60/3600})) AND (g.phot_variable_flag = 'VARIABLE' OR g.non_single_star > 0) "
          "ORDER BY sep", "flagged_within_60arcsec")
G["variable_or_nss_within_60arcsec"] = near.astype(str).to_dict("records")
print("flagged within 60\":\n", near.to_string())
ec = gq("SELECT e.source_id, e.frequency, e.global_ranking, e.model_type FROM gaiadr3.vari_eclipsing_binary AS e "
        f"WHERE e.source_id IN ({','.join(str(s) for s in near.source_id)})" if len(near) and "source_id" in near else
        "SELECT TOP 0 * FROM gaiadr3.vari_eclipsing_binary", "ecl_within_60arcsec")
if len(ec) and "frequency" in ec:
    ec["period_d"] = 1 / ec.frequency
G["ecl_within_60arcsec"] = ec.astype(str).to_dict("records")
print("Gaia EBs within 60\":\n", ec.to_string())
json.dump(G, open(os.path.join(OUT, "gaia_vari_nss.json"), "w"), indent=1)

# ---------------- VSX ----------------
V = {"queried_utc": NOW}
c = SkyCoord(260.970484934795, -51.3242024657333, unit="deg")  # J2000
v = Vizier(columns=["**", "_r"], row_limit=200)
res = v.query_region(c, radius=5 * u.arcmin, catalog="B/vsx/vsx")
if len(res):
    t = res[0].to_pandas().sort_values("_r")
    t.to_csv(os.path.join(RAW, "vsx_vizier_5arcmin.csv"), index=False)
    V["vizier_B_vsx_5arcmin"] = t[[c_ for c_ in ("_r", "OID", "Name", "Type", "Period", "max", "min", "n_max") if c_ in t]].astype(str).to_dict("records")
else:
    V["vizier_B_vsx_5arcmin"] = []
print("VSX (VizieR) within 5':", V["vizier_B_vsx_5arcmin"])
try:
    r = requests.get("https://www.aavso.org/vsx/index.php", params={"view": "api.list", "ra": 260.970485, "dec": -51.324202,
                                                                      "radius": 5 / 60, "format": "json"}, timeout=120)
    open(os.path.join(RAW, "vsx_api_5arcmin.json"), "w").write(r.text[:100000])
    V["aavso_api_5arcmin"] = {"http_status": r.status_code, "content": r.text[:3000]}
except Exception as e:
    V["aavso_api_5arcmin"] = {"error": str(e)[:300]}
print("VSX API:", V["aavso_api_5arcmin"])
json.dump(V, open(os.path.join(OUT, "vsx_status.json"), "w"), indent=1)
