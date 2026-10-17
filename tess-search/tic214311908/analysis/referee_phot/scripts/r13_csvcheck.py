"""Referee r13: audit of the investigator's deliverable lc_M3V_all.csv (and lc_pair_all.csv / lc_Mira_all.csv):
columns, time correction vs the referee's (r02), normalisation, quality fractions per sector, in-transit coverage,
and a simple depth measurement (box at the referee global ephemeris/shape) on flux_norm (model B) and flux_norm_modelA.
Output: out/csvcheck.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import trap_shape
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]
M = pd.read_csv(f"{R}/../phot/lc_M3V_all.csv")
out = dict(columns=list(M.columns), n_rows=len(M), n_quality_ok=int(M.quality_ok.sum()))
per = {}
for s in SECTORS:
    d = M[M.sector == s]
    z = np.load(f"{R}/data/tcorr_s{s}.npy")
    corr_ref = np.interp(d.time_btjd_tesscut.values, z[0], z[1])
    corr_inv = (d.time_btjd_tdb.values - d.time_btjd_tesscut.values) * 86400
    ok = (d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_norm)
    t = d.time_btjd_tdb.values[ok]; f = d.flux_norm.values[ok]; ex = d.cadence_s.values[0] / 86400
    y = trap_shape(t, T0 + np.round((t - T0) / P) * P, T14, FL, ex)
    A = np.column_stack([np.ones(len(t)), -y]); c, *_ = np.linalg.lstsq(A, f, rcond=None)
    okA = (d.quality_ok == 1) & np.isfinite(d.flux_norm_modelA)
    tA = d.time_btjd_tdb.values[okA]; fA = d.flux_norm_modelA.values[okA]
    yA = trap_shape(tA, T0 + np.round((tA - T0) / P) * P, T14, FL, ex)
    cA, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(tA)), -yA]), fA, rcond=None)
    oot = y == 0
    # residual trend check: mean of flux_norm in +-0.1..0.3 d around transits (should be 1)
    ph = np.abs(((t - T0) / P + 0.5) % 1 - 0.5) * P
    near = (ph > 0.06) & (ph < 0.3)
    per[s] = dict(time_corr_inv_minus_ref_s=float(np.nanmedian(corr_inv - corr_ref)), time_corr_rms_s=float(np.nanstd(corr_inv - corr_ref)),
                  frac_quality_ok=float(d.quality_ok.mean()), n_in_transit_ok=int(np.sum(y > 0.5)),
                  depth_ppt_modelB=float(c[1] * 1e3), depth_e_s_modelB=float(c[1] * 935), depth_ppt_modelA=float(cA[1] * 1e3),
                  baseline=float(c[0]), near_transit_mean_minus1_ppt=float((np.mean(f[near]) - 1) * 1e3), near_err_ppt=float(np.std(f[near]) / np.sqrt(near.sum()) * 1e3),
                  oot_rstd_ppt=float(rstd(f[oot]) * 1e3), err_median_ppt=float(np.median(d.flux_err[ok]) * 1e3))
    print(f"S{s}: time corr inv-ref {per[s]['time_corr_inv_minus_ref_s']:+.2f} s; quality_ok {per[s]['frac_quality_ok']:.2f}; in-transit ok {per[s]['n_in_transit_ok']}; depth B {c[1]*1e3:.2f} ppt ({c[1]*935:.1f} e/s), A {cA[1]*1e3:.2f} ppt; "
          f"near-transit offset {per[s]['near_transit_mean_minus1_ppt']:+.2f}+-{per[s]['near_err_ppt']:.2f} ppt; OOT rstd {per[s]['oot_rstd_ppt']:.1f} vs flux_err median {per[s]['err_median_ppt']:.1f} ppt")
out["per_sector"] = per
save_json(out, f"{R}/out/csvcheck.json")
