"""Collect all numbers into ../results.json."""
import os, json, glob
import numpy as np
from common import ROOT, SECTORS


def load(p):
    p = os.path.join(ROOT, p)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    R = {}
    R["static"] = {str(s): {k: v for k, v in load(f"static/static_s{s}.json").items() if k != "field_star_offsets"} for s in SECTORS}
    nom = load("lc/nominal/depths.json")
    R["nominal"] = nom
    R["nominal_with_regressors"] = load("lc/nominal/depths_regs.json")
    var = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "lc", "*", "depths.json"))):
        d = json.load(open(f)); c = d["combined"]
        var[d["tag"]] = dict(f_target=c["f_target"]["f_best"], f_lo=c["f_target"]["f_lo"], f_hi=c["f_target"]["f_hi"],
                             chi2_f1=c["f_target"]["chi2_f1"], chi2_f0=c["f_target"]["chi2_f0"],
                             dchi2_target_minus_nb=c["delta_chi2_target_minus_neighbour"],
                             depth_t_psf=c["depth_t_psfnorm"], depth_n_cat=c["depth_n_catnorm"], depth_pair=c["depth_pair"],
                             per_sector_f=[d["sectors"][str(s)]["f_target"] for s in SECTORS],
                             per_sector_f_err=[d["sectors"][str(s)]["f_target_err"] for s in SECTORS])
    R["variants"] = var
    inj = load("injection/injection_results.json")
    if inj:
        summ = {}
        for host in ("target", "neighbour"):
            for c in sorted(set(r["inj_prf"] for r in inj)):
                rr = [r for r in inj if r["host"] == host and r["inj_prf"] == c]
                summ[f"{host}|{c}"] = dict(resp_t_mean=float(np.mean([r["resp_t"] for r in rr])), resp_n_mean=float(np.mean([r["resp_n"] for r in rr])),
                                           f_target_min=float(np.min([r["f_target_resp"] for r in rr])), f_target_max=float(np.max([r["f_target_resp"] for r in rr])),
                                           n=len(rr))
        R["injection_summary"] = summ
    R["oot_validation"] = load("lc/nominal/oot_validation.json")
    R["centroid_test"] = load("centroid/centroid_test.json")
    R["secondary_oddeven"] = (load("lc/nominal/secondary_oddeven.json") or {}).get("combined")
    R["secondary_oddeven_per_sector"] = (load("lc/nominal/secondary_oddeven.json") or {}).get("per_sector")
    R["neighbour_own_lcs"] = load("neighbour_lcs/neighbour_lcs.json")
    ap = {}
    tot = {}
    for s in SECTORS:
        a = load(f"aperture/aperture_s{s}.json")
        if a is None: continue
        ap[str(s)] = dict(named={k: v for k, v in a["apertures"].items() if not k.startswith("pix_")}, tests=a["tests"],
                          free_position=a["free_position"], other_hosts=a["other_hosts_pixelmap_eps0.05"][:6])
        for setn, tt in a["tests"].items():
            for e, q in tt.items():
                key = f"{setn}|{e}"
                tot.setdefault(key, dict(chi2_target=0, chi2_nb=0, dof=0))
                tot[key]["chi2_target"] += q["chi2_target"]; tot[key]["chi2_nb"] += q["chi2_nb"]; tot[key]["dof"] += q["dof"]
    R["aperture"] = ap; R["aperture_chi2_totals"] = tot
    R["comparison_with_transit_folder"] = load("comparison_transit.json")
    json.dump(R, open(os.path.join(ROOT, "results.json"), "w"), indent=1, default=float)
    print("written", os.path.join(ROOT, "results.json"))


if __name__ == "__main__":
    main()
