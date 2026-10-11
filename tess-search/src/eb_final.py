"""Final characterisation table for confirmed new eclipsing binaries."""
import os, sys, json, warnings
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
import stage2, known
from eb_analysis import eclipse_times

W = "/home/user/work"


def primary_ephemeris(tic, P0, T00, secs):
    ev = eclipse_times(secs)
    prim = []
    for e in ev:
        ph = ((e["t"] - T00) / P0) % 1.0
        if min(ph, 1 - ph) < 0.04: prim.append(e["t"])
    prim = np.array(sorted(prim))
    if len(prim) < 3: return None
    n = np.round((prim - T00) / P0)
    A = np.vstack([np.ones_like(n), n]).T
    (T0, P), *_ = np.linalg.lstsq(A, prim, rcond=None)
    resid = prim - (T0 + n * P)
    cov = np.linalg.inv(A.T @ A) * max(np.var(resid, ddof=2), (1 / 1440) ** 2)
    return dict(P=float(P), P_err=float(np.sqrt(cov[1, 1])), T0=float(T0), T0_err=float(np.sqrt(cov[0, 0])),
                n=int(len(prim)), rms_min=float(np.std(resid) * 1440), cycles=int(n.max() - n.min()))


def gaia(ra, dec):
    from astroquery.gaia import Gaia
    q = f"""SELECT TOP 1 source_id, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe, radial_velocity, rv_nb_transits,
           rv_amplitude_robust, rv_chisq_pvalue, non_single_star, teff_gspphot,
           DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep
           FROM gaiadr3.gaia_source WHERE CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{3 / 3600}))=1 ORDER BY phot_g_mean_mag"""
    r = Gaia.launch_job(q).get_results().to_pandas()
    return r.iloc[0].to_dict() if len(r) else {}


def novelty(ra, dec):
    out = {}
    out["eb_catalogues_21"] = [f"{h['cat']} ({h['sep']:.0f}\")" for h in known.vizier_cone(ra, dec, 21)]
    out["kostov_unvetted"] = bool(known.kostov_unvetted(ra, dec))
    v, n = known.vizier_all(ra, dec, 15)
    out["vizier_keyword_hits"] = [x for x in v if not x.startswith(("IV/38", "IV/39", "J/AJ/155/39", "I/355/epphot", "I/358/varisum"))]
    out["n_vizier_tables_15"] = n
    out["simbad"] = [f"{x['main_id']} ({x['otype']}, {float(x['sep']):.1f}\")" for x in known.simbad(ra, dec, 20)]
    return out


if __name__ == "__main__":
    rows = []
    for tic in [int(x) for x in sys.argv[1:]]:
        r = json.load(open(f"{W}/vet/EB_TIC{tic}/eb.json")); s = r["star"]
        secs = stage2.load_all(tic)
        eph = primary_ephemeris(tic, r["P"], r["t0"], secs)
        g = gaia(s["ra"], s["dec"])
        nov = novelty(s["ra"], s["dec"])
        rows.append(dict(tic=tic, ra=s["ra"], dec=s["dec"], Tmag=s["Tmag"], Teff_TIC=s["Teff"], R_TIC=s["rad"],
                         P=(eph or {}).get("P", r["P"]), P_err=(eph or {}).get("P_err"), T0_BTJD=(eph or {}).get("T0", r["t0"]),
                         T0_err=(eph or {}).get("T0_err"), n_primaries=(eph or {}).get("n"), timing_rms_min=(eph or {}).get("rms_min"),
                         depth_primary=r["depth_prim"], depth_secondary=r["depth_sec"], phase_secondary=r["phase_sec"],
                         ecosw=float(np.pi / 2 * (r["phase_sec"] - 0.5)), sectors=r["sectors"],
                         loc=[(L["sector"], round(L["offset_arcsec"], 1), round(L["offset_sigma"], 1)) for L in r.get("localisation", [])],
                         gaia=g, **nov))
        print(tic, "done", flush=True)
    out = os.environ.get("EB_FINAL_OUT", f"{W}/vet/eb_final.json")
    json.dump(rows, open(out, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
