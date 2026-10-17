"""s08: Gaia DR3 variability information for all sources within 4 arcmin of the M3V (network: Gaia archive TAP).
gaia_source (G, BP-RP, parallax, phot_variable_flag) + vari_classifier_result + vari_rrlyrae + vari_eclipsing_binary
+ vari_rotation_modulation, each queried separately and merged.
Output: data/gaia/gaia_vari_4arcmin.csv, out/gaia_variables.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from astroquery.gaia import Gaia
RA, DEC = 260.9696737, -51.3238212
R = 4 / 60
base = f"""SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.parallax_error, g.phot_variable_flag,
 g.pmra, g.pmdec, DISTANCE(POINT(g.ra, g.dec), POINT({RA}, {DEC}))*3600 AS sep
FROM gaiadr3.gaia_source AS g WHERE 1 = CONTAINS(POINT(g.ra, g.dec), CIRCLE({RA}, {DEC}, {R}))"""
t = Gaia.launch_job_async(base).get_results().to_pandas()
print(len(t), "sources")
def sub(table, cols):
    q = f"""SELECT v.source_id, {', '.join('v.' + c for c in cols)} FROM gaiadr3.{table} AS v JOIN gaiadr3.gaia_source AS g ON g.source_id = v.source_id
    WHERE 1 = CONTAINS(POINT(g.ra, g.dec), CIRCLE({RA}, {DEC}, {R}))"""
    try:
        r = Gaia.launch_job_async(q).get_results().to_pandas(); print(table, len(r)); return r
    except Exception as e:
        print(table, "failed", str(e)[:200]); return None
for table, cols in [("vari_classifier_result", ["best_class_name", "best_class_score"]),
                    ("vari_rrlyrae", ["pf", "p1_o", "peak_to_peak_g", "best_classification"]),
                    ("vari_eclipsing_binary", ["frequency", "global_ranking"]),
                    ("vari_rotation_modulation", ["best_rotation_period", "max_activity_index"])]:
    r = sub(table, cols)
    if r is not None and len(r): t = t.merge(r, on="source_id", how="left")
t = t.sort_values("sep")
t.to_csv(f"{ROOT}/data/gaia/gaia_vari_4arcmin.csv", index=False)
cols = [c for c in ["best_class_name", "pf", "p1_o", "frequency", "best_rotation_period"] if c in t.columns]
v = t[(t.phot_variable_flag == "VARIABLE") | t[cols].notna().any(axis=1)] if cols else t[t.phot_variable_flag == "VARIABLE"]
print(len(v), "with variability info")
print(v.drop(columns=["ra", "dec", "pmra", "pmdec"]).to_string())
save_json(dict(n_sources=len(t), variables=v.to_dict(orient="records")), f"{ROOT}/out/gaia_variables.json")
