"""Neighbours able to mimic each signal, and catalogued variables/EBs nearby with a matching period."""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
Gaia.ROW_LIMIT = -1
T = {90086266: (276.27770, -42.64403, 2.339469, 0.031, 11.11),
     61332759: (287.61974, -44.73357, 0.740299, 0.0011, 11.16),
     214311906: (260.97765, -51.32369, 5.69149, 0.0024, 11.26),
     55829034: (272.03163, -35.69140, 0.616146, 0.0016, 11.60),
     61141412: (286.96474, -43.00969, 6.5285, 0.0010, 10.02)}
def tmag(g, c):
    c = np.where(np.isfinite(c), c, 1.0)
    return g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
def pmatch(p1, p2, tol=0.01):
    if not (np.isfinite(p1) and np.isfinite(p2)) or p2 <= 0: return ""
    r = p1 / p2
    for k, lab in ((1, "P"), (2, "P/2"), (0.5, "2P"), (3, "P/3"), (1 / 3, "3P")):
        if abs(r - k) < tol * k: return lab
    return ""
out = {}
if __name__ != "__main__": T = {}
for tic, (ra, dec, P, depth, Tt) in T.items():
    q = f"""SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag AS G, g.bp_rp, g.parallax, g.pmra, g.pmdec, g.ruwe,
           g.phot_variable_flag, DISTANCE(POINT({ra},{dec}), POINT(g.ra,g.dec))*3600 AS sep
           FROM gaiadr3.gaia_source g WHERE CONTAINS(POINT(g.ra,g.dec), CIRCLE({ra},{dec},{240/3600}))=1 AND g.phot_g_mean_mag < 20"""
    g = Gaia.launch_job(q).get_results().to_pandas()
    g["T"] = tmag(g.G, g.bp_rp)
    tgt = g.sep.idxmin()
    g["dT"] = g["T"] - g.loc[tgt, "T"]
    # can a 100% (50%) eclipse of this star produce the observed diluted depth?
    g["max_depth_100"] = 10 ** (-0.4 * g.dT) / (1 + 10 ** (-0.4 * g.dT))
    g["can_mimic_50"] = 0.5 * g.max_depth_100 >= depth
    ids = ",".join(str(int(x)) for x in g.source_id)
    eb = Gaia.launch_job(f"SELECT source_id, frequency, global_ranking, derived_primary_ecl_depth, model_type FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN ({ids})").get_results().to_pandas()
    vc = Gaia.launch_job(f"SELECT source_id, best_class_name, best_class_score FROM gaiadr3.vari_classifier_result WHERE source_id IN ({ids})").get_results().to_pandas()
    g = g.merge(eb, on="source_id", how="left").merge(vc, on="source_id", how="left")
    g["P_eb"] = 1 / g.frequency
    g["eb_match"] = [pmatch(P, p) for p in g.P_eb]
    mim = g[g.can_mimic_50].sort_values("sep")
    print(f"\n=== TIC {tic}: P {P} d, diluted depth {depth*1e6:.0f} ppm; Gaia sources within 4' (G<20): {len(g)}; able to mimic with a 50% eclipse: {len(mim)}")
    print(mim[["source_id", "sep", "G", "bp_rp", "T", "dT", "parallax", "ruwe", "phot_variable_flag", "best_class_name", "P_eb", "eb_match"]].head(15).round(3).to_string(index=False))
    ebs = g[g.P_eb.notna() | g.best_class_name.notna()]
    print("  Gaia variables within 4':", ebs[["source_id", "sep", "G", "best_class_name", "P_eb", "eb_match"]].round(4).to_string(index=False) if len(ebs) else "none")
    # VSX and ASAS-SN within 5'
    c = SkyCoord(ra, dec, unit="deg")
    for cat, pcol, name in (("B/vsx/vsx", "Period", "VSX"), ("II/366/catalog", "Per", "ASAS-SN")):
        try:
            r = Vizier(columns=["**", "_r"], row_limit=200).query_region(c, radius=5 * u.arcmin, catalog=cat)
            if len(r):
                t = r[0].to_pandas()
                cols = [x for x in ("_r", "Name", "ID", "Type", "VarType", pcol, "max", "Vmag", "Amp") if x in t.columns]
                t["match"] = [pmatch(P, float(p)) if p == p else "" for p in t[pcol]]
                print(f"  {name} within 5':\n", t[cols + ["match"]].to_string(index=False))
            else:
                print(f"  {name} within 5': none")
        except Exception as e:
            print(f"  {name} query failed: {e}")
    g.to_csv(f"out/{tic}/gaia_4arcmin.csv", index=False)
