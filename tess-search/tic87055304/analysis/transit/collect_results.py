"""Collect all results of this folder into results.json (numbers quoted in README.md)."""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(HERE, n))) if os.path.exists(os.path.join(HERE, n)) else None


def pv(p, k, nd=5):
    return dict(med=round(p[k]["med"], nd), minus=round(p[k]["lo"], nd), plus=round(p[k]["hi"], nd))


res = dict(target="TIC 87055304", area="transit model, shape, jackknife, odd/even, depths, timing (+ pixel-level depth/localisation cross-check)")
runs = {}
for nm in ("free", "dwarf", "subgiant", "dil051", "nbhost", "oddeven", "drop13", "drop39", "drop66", "drop93", "drop104"):
    js = J(f"mcmc_{nm}.json")
    if js is None:
        continue
    p = js["params"]
    r = dict(rho_prior=js["rho_prior"], dilution=js["dil"], sectors=js["sectors"], converged=js["converged"], n_steps=js["nsteps"],
             n_walkers=js["nwalkers"], max_tau=max(js["tau"].values()), n_over_tau=js["n_over_tau"], acceptance=js["acceptance"],
             best_lnL=js["best_lnL"], best_chi2=js["best_chi2"], ndata=js["ndata"])
    for k in ("P", "T0", "k", "b", "rho", "a", "inc", "depth_ppm", "k2_ppm", "T14_h", "T23_h", "u1", "u2"):
        r[k] = pv(p, k, 7 if k in ("P", "T0") else 5)
    r["T0_BJD_TDB"] = round(p["T0"]["med"] + 2457000, 6)
    if nm == "oddeven":
        r["k_even"] = pv(p, "k_even"); r["k_odd_minus_even"] = pv(p, "k_odd_minus_even", 6)
        r["depth_ratio_odd_even"] = pv(p, "depth_ratio_odd_even"); r["diff_sigma"] = js["oddeven_diff_sigma"]
    runs[nm] = r
res["mcmc"] = runs
if "free" in runs and "dwarf" in runs:
    res["density_prior_test"] = dict(dlnL_free_minus_dwarf=runs["free"]["best_lnL"] - runs["dwarf"]["best_lnL"],
                                     dchi2_equiv=2 * (runs["free"]["best_lnL"] - runs["dwarf"]["best_lnL"]))
    if "subgiant" in runs:
        res["density_prior_test"]["dlnL_free_minus_subgiant"] = runs["free"]["best_lnL"] - runs["subgiant"]["best_lnL"]
for n in ("vshape.json", "ld_test.json", "depths.json", "ttv_raw.json", "ttv_raw_w040.json", "sector_timing.json", "jackknife.json",
          "detrend_sys.json", "aperture_depth.json", "pixel_depth.json", "prf_localize.json", "localize_checks.json",
          "localize_validate.json", "localize_twosource.json", "localization_ratio_maps.json", "nb_lc_check.json", "secondary.json",
          "ld_prior.json"):
    js = J(n)
    if js is None:
        continue
    key = n.replace(".json", "")
    if key == "pixel_depth":
        for s in js["sectors"].values():
            s.pop("dF_map", None); s.pop("dF_err_map", None)
    if key in ("ttv_raw", "ttv_raw_w040"):
        js = {k: v for k, v in js.items()}
    if key == "depths":
        js.pop("per_transit_detrended", None)
    res[key] = js
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
print("wrote results.json with keys", list(res.keys()))
