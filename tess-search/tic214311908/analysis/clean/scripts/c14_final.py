"""c14: final adopted depth (stat + sys), planet parameters and odd/even from the cleaned light curves.

Statistical: posterior medians of the batman fits with the single-star density prior (c04) for INV and REF, combined
by GLS with their noise correlation rho (pooled INV-REF correlation of no-injection depths at the same fake transits,
c05). Odd/even: fixed-shape group fits (c07 shape) with beta-scaled errors, per pipeline and combined.
Systematic budget (each a 1-sigma half-range, added in quadrature):
  (1) method: half the INV-REF difference of the global depth;
  (2) PRF template error: rms over sectors of (response - 1) for the 0.07-px shifted and 0.15-px blurred injections
      (both pipelines), times the depth;
  (3) normalisation F_M3V = 979 e-/s: +-2 % (XP stat) and the route spread (966-1035 -> -1.3 %/+5.7 %), taken as +-3 %;
  (4) baseline / detrending version (INV raw vs det, c11): half the difference;
  (5) window selection (REF with / without the scattered-light event E58, c07): half the difference;
  (6) residual RRab signal (largest residual first-harmonic amplitude 0.8 ppt, S93) projected on the transit: the rms
      of the per-transit depth change it causes divided by sqrt(N transits) -> taken as 0.2 ppt (upper bound).
Output: out/final_depth.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from c04_global import load
from c07_epoch import shape_from
from tfit import Model, mid_depth, fit_depth_fixed_shape, R_STAR, R_EARTH_RSUN
import pandas as pd


def gls2(a, ea, b, eb, rho):
    C = np.array([[ea ** 2, rho * ea * eb], [rho * ea * eb, eb ** 2]]); Ci = np.linalg.inv(C); o = np.ones(2)
    v = 1 / (o @ Ci @ o); w = v * (Ci @ o)
    return float(w @ [a, b]), float(np.sqrt(v)), w.tolist()


def odd_even(pipe):
    W = load(pipe); sh = shape_from(pipe)
    D0 = mid_depth(sh["k"], sh["b"], sh["aR"], sh["u1"], sh["u2"], sh["P"])
    frac = (1 - Model(W.t, W.ex).flux(sh["T0"], sh["P"], sh["k"], sh["b"], sh["aR"], sh["u1"], sh["u2"])) / D0
    g = (W.epoch % 2 != 0).astype(int)
    dep, err, chi2, ug = fit_depth_fixed_shape(W, frac, groups=g)
    return dict(even=[float(dep[0] * 1e3), float(err[0] * 1e3)], odd=[float(dep[1] * 1e3), float(err[1] * 1e3)],
                odd_minus_even=[float((dep[1] - dep[0]) * 1e3), float(np.hypot(err[0], err[1]) * 1e3)])


def main():
    G = {p: json.load(open(f"{CL}/out/global_{p}_single.json")) for p in ("inv", "ref")}
    GF = {p: json.load(open(f"{CL}/out/global_{p}_free.json")) for p in ("inv", "ref") if os.path.exists(f"{CL}/out/global_{p}_free.json")}
    # pooled noise correlation
    a, b = [], []
    for s in SECTORS:
        d = pd.read_csv(f"{CL}/out/inject_pix_s{s}.csv")
        x = d[d.variant == "invprf"][["tc", "inv_raw_depth_fix_noinj", "ref_depth_fix_noinj"]].dropna()
        sx = x.inv_raw_depth_fix_noinj.std(); sy = x.ref_depth_fix_noinj.std()
        a += list(x.inv_raw_depth_fix_noinj / sx); b += list(x.ref_depth_fix_noinj / sy)
    rho = float(np.corrcoef(a, b)[0, 1])
    def pm(q): return q["med"], 0.5 * (q["lo"] + q["hi"])
    di, ei = pm(G["inv"]["depth_mid_ppt"]); dr, er = pm(G["ref"]["depth_mid_ppt"])
    dc, ec, w = gls2(di, ei, dr, er, rho)
    out = dict(rho_inv_ref=rho, inv=[di, ei], ref=[dr, er], combined_stat=[dc, ec], weights=w)
    # systematics
    I = pd.read_csv(f"{CL}/out/inject_summary.csv")
    mis = I[I.variant.isin(["sx+", "sx-", "sy+", "sy-", "blur"])]
    prf_rms = float(np.sqrt(np.mean(np.concatenate([(mis.inv_resp_mean - 1) ** 2, (mis.ref_resp_mean - 1) ** 2]))))
    M = json.load(open(f"{CL}/out/modes_compare.json"))
    DE = json.load(open(f"{CL}/out/depth_epoch.json"))
    sys_ = {"method_inv_vs_ref": abs(di - dr) / 2,
            "prf_template": prf_rms * dc,
            "normalisation_979": 0.03 * dc,
            "baseline_version": abs(M["inv_all_raw"]["depth"] - M["inv_all_det"]["depth"]) / 2,
            "window_selection_E58": abs(DE["ref"]["all"][0] - DE["ref_noE58"]["all"][0]) / 2,
            "residual_rrab": 0.2}
    tot = float(np.sqrt(sum(v ** 2 for v in sys_.values())))
    out.update(sys_budget=sys_, sys_total=tot, prf_template_rms_frac=prf_rms,
               adopted=dict(depth_ppt=dc, stat=ec, sys=tot, total=float(np.hypot(ec, tot)),
                            flux_loss_e_s=dc * 0.979, note="mid-transit depth of the batman model (density prior) in ppt of the M3V flux (979 e-/s, PRF units)"))
    # radius ratio & radius (combine k the same way; sys from depth sys: dk/k = 0.5 d(depth)/depth)
    ki, eki = pm(G["inv"]["k"]); kr, ekr = pm(G["ref"]["k"])
    kc, ekc, _ = gls2(ki, eki, kr, ekr, rho)
    ksys = 0.5 * tot / dc * kc
    Rp = kc * R_STAR[0] / R_EARTH_RSUN
    eRp_stat = Rp * np.hypot(ekc / kc, R_STAR[1] / R_STAR[0]); eRp_sys = Rp * ksys / kc
    out["k"] = dict(inv=[ki, eki], ref=[kr, ekr], combined=[kc, ekc], sys=ksys)
    out["Rp_RE"] = dict(value=float(Rp), stat_incl_Rstar=float(eRp_stat), sys=float(eRp_sys))
    for par in ("b", "aR", "inc", "T14_h", "flat_frac", "Teq_K", "rho"):
        out[par] = {p: pm(G[p][par]) for p in G}
    out["free_density"] = {p: {par: pm(GF[p][par]) for par in ("rho", "k", "b", "T14_h", "flat_frac", "depth_mid_ppt")} | {"frac_grazing": GF[p]["frac_grazing"]} for p in GF}
    oe = {p: odd_even(p) for p in ("inv", "ref")}
    oc, eoc, _ = gls2(oe["inv"]["odd_minus_even"][0], oe["inv"]["odd_minus_even"][1], oe["ref"]["odd_minus_even"][0], oe["ref"]["odd_minus_even"][1], rho)
    oe["combined"] = [oc, eoc]
    out["odd_even"] = oe
    jd(out, f"{CL}/out/final_depth.json")
    print(json.dumps({k: out[k] for k in ("rho_inv_ref", "inv", "ref", "combined_stat", "sys_budget", "sys_total", "adopted", "k", "Rp_RE", "odd_even")}, indent=1))
    for p in GF: print(p, "free:", out["free_density"][p])


if __name__ == "__main__":
    main()
