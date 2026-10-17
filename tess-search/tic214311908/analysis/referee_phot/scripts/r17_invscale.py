"""Referee r17: is the investigator's depth normalisation (M3V = 935 e-/s, TIC Tmag 13.0135 at 15000 e/s for T=10)
on the same scale as their PRF-fitted transit flux loss (e-/s of a unit-sum SPOC PRF)?
Uses the investigator's OWN field-star PRF fluxes (phot/out/wcs_s*.json, 'free_fluxes') for the isolated field stars
selected by the referee (r01b: no brighter-than-T+1 neighbour within 1 px, >1.5 px inside the cutout, TIC T < 14),
and compares them with (a) the TIC catalogue flux and (b) the referee's own PRF-fit fluxes (r01b) for the same stars.
Then converts the flux loss into a depth with the M3V flux in PRF units, from
  (1) TIC Tmag x investigator's median PRF/TIC ratio,
  (2) Gaia-XP synthetic TESS photon rate calibrated on the PRF fluxes of field stars (r09),
  (3) the 'star' agent's SED-synthetic Tmag 13.086 x PRF/TIC ratio.
Output: out/invscale.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
FS = pd.read_csv(f"{R}/out/fluxscale_stars.csv")
XP = json.load(open(f"{R}/out/xpflux.json"))
out = {}
rat_all = []
for s in SECTORS:
    W = json.load(open(f"{R}/../phot/out/wcs_s{s}.json"))
    ff = {int(k): v for k, v in W["free_fluxes"].items()}
    d = FS[(FS.sector == s) & FS.iso & FS.inside & (FS.Ttic < 14.0)].copy()
    d["inv"] = d.source_id.map(ff)
    d = d[np.isfinite(d.inv)]
    d["inv_tic"] = d.inv / tflux(d.Ttic); d["inv_ref"] = d.inv / d.fit
    out[s] = dict(n=len(d), inv_over_tic_median=float(np.median(d.inv_tic)), inv_over_tic_err=float(1.2533 * rstd(d.inv_tic) / np.sqrt(len(d))),
                  inv_over_ref_median=float(np.median(d.inv_ref)), inv_over_ref_rstd=float(rstd(d.inv_ref)),
                  investigator_faint_scale=float(W["faint_scale"]))
    rat_all.append(d.assign(sector=s))
    print(f"S{s}: n {len(d)}; investigator PRF flux / TIC flux = {out[s]['inv_over_tic_median']:.3f} +- {out[s]['inv_over_tic_err']:.3f}; "
          f"investigator / referee PRF flux = {out[s]['inv_over_ref_median']:.3f} (rstd {out[s]['inv_over_ref_rstd']:.3f}); investigator's own faint-star scale {W['faint_scale']:.3f}")
A = pd.concat(rat_all)
k = float(np.median([out[s]["inv_over_tic_median"] for s in SECTORS])); ke = float(np.std([out[s]["inv_over_tic_median"] for s in SECTORS]) / np.sqrt(6))
F_tic = tflux(13.0135)
F1 = k * F_tic
F2 = XP["ext"]["F_M3V_prf_mean"]; F2b = XP["noext"]["F_M3V_prf_mean"]
F3 = k * tflux(13.0855)
dF_inv, dF_inv_err = 15.4, 0.8
dF_ref = json.load(open(f"{R}/out/transit_fit.json"))["per_sector_mean"]
res = dict(per_sector=out, prf_over_tic_overall=k, prf_over_tic_err=ke, F_M3V_tic_catalogue=F_tic,
           F_M3V_prf_from_TIC=F1, F_M3V_prf_from_XP=F2, F_M3V_prf_from_XP_noext=F2b, F_M3V_prf_from_SED_T=F3)
for lab, F in (("TIC_cat_935", 935.0), ("TIC_x_scale", F1), ("XP", F2), ("XP_noext", F2b), ("SED_T", F3), ("Gaia_colour_837", 837.0)):
    res[f"depth_ppt_inv_{lab}"] = (dF_inv / F * 1e3, dF_inv_err / F * 1e3)
    res[f"depth_ppt_ref_{lab}"] = (dF_ref["dF"] / F * 1e3, dF_ref["err"] / F * 1e3)
save_json(res, f"{R}/out/invscale.json")
print(f"overall investigator PRF/TIC flux ratio {k:.3f} +- {ke:.3f}; M3V in PRF units: TIC x ratio {F1:.0f}, XP {F2:.0f} ({F2b:.0f} no tail), SED-T x ratio {F3:.0f} e/s")
for lab in ("TIC_cat_935", "TIC_x_scale", "XP", "XP_noext", "SED_T", "Gaia_colour_837"):
    print(f"  depth with F={lab}: investigator flux loss -> {res['depth_ppt_inv_'+lab][0]:.2f}+-{res['depth_ppt_inv_'+lab][1]:.2f} ppt; referee -> {res['depth_ppt_ref_'+lab][0]:.2f}+-{res['depth_ppt_ref_'+lab][1]:.2f} ppt")
