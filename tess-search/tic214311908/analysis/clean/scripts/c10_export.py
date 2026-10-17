"""c10: final adopted cleaned light curves (CSV) of both pipelines.

lc_M3V_INV_clean.csv  (INV = phot model-B method, field variables excluded from every field-wide estimator; all six
                       sectors, every finite cadence)
lc_M3V_REF_clean.csv  (REF = referee local-window method, cleaned auxiliary regressors; one strict window per transit)
lc_M3V_REF_tiled.csv  (REF method run as contiguous tiles, for out-of-transit work; see README caveats)
Column definitions are written to out/lc_columns.json and the README.
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
import pandas as pd

EPH = json.load(open(f"{CL}/out/ephemeris_final.json"))["combined"]
P, T0 = EPH["P_d"], EPH["T0_epoch0_BTJD_TDB"]
NOISE = json.load(open(f"{CL}/out/noise_calib.json"))


def intr(t, half=0.1):
    return np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < half


def epoch_of(t):
    return np.round((t - T0) / P).astype(int)


def inv():
    rows = []
    for s in SECTORS:
        z = np.load(f"{CL}/data/lcinv_all_s{s}.npz")
        t = z["t_tdb"]; fin = np.isfinite(t)
        det = z["det_M3V_B"]; ok = z["quality_ok"].astype(bool) & np.isfinite(det)
        sd = robust_std(det[ok & ~z["intransit_window"].astype(bool)])
        hi = ok & (det > 1 + 4 * sd); fl = hi | np.roll(hi, 1) | np.roll(hi, -1)
        # residual RRab model (3 harmonics at the Gaia period, fitted out of transit per sector) for optional subtraction
        m = ok & ~fl & ~intr(t)
        ph = 2 * np.pi * (t - TMAX_RR) / P_RR
        X = np.column_stack([np.ones(len(t))] + [f(h * ph) for h in (1, 2, 3) for f in (np.cos, np.sin)])
        c, *_ = np.linalg.lstsq(X[m], det[m] - 1, rcond=None)
        rr = X[:, 1:] @ c[1:]
        cad = int(round(float(z["exptime"]) * 86400))
        beta = NOISE["inv_raw"][str(s)]["beta"]
        rows.append(pd.DataFrame(dict(
            time_btjd_tdb=t[fin], time_bjd_tdb=t[fin] + 2457000.0, sector=s, cadence_s=cad,
            quality_ok=(ok & ~fl)[fin].astype(int), flux_norm=det[fin], flux_err=z["err_M3V_B"][fin],
            beta_transit=beta, flux_raw_norm=z["raw_M3V_B"][fin], trend=z["trend_M3V_B"][fin], dflux_e_s=z["dF_M3V_B"][fin],
            flare_flag=fl[fin].astype(int), rrab_resid_model=rr[fin], flux_norm_modelA=z["det_M3V_A"][fin],
            flux_err_modelA=z["err_M3V_A"][fin], epoch=epoch_of(t[fin]), in_transit_window=intr(t[fin]).astype(int),
            tess_quality=z["q"][fin], bkg_e_s_px=z["bkg"][fin], jitter_dx_px=z["dx"][fin], jitter_dy_px=z["dy"][fin],
            time_btjd_tesscut=z["t_tess"][fin])))
    d = pd.concat(rows, ignore_index=True)
    d.to_csv(f"{CL}/lc_M3V_INV_clean.csv", index=False, float_format="%.7f")
    return len(d)


def ref():
    rows = []
    wid = 0
    for s in SECTORS:
        z = np.load(f"{CL}/data/refwin_all_s{s}.npz")
        ex = int(round(float(z["exptime"]) * 86400))
        beta = NOISE["ref"][str(s)]["beta"]
        for j, (tc, ep) in enumerate(zip(z["strict_tc"], z["strict_epoch"])):
            m = z["strict_win"] == j
            t = z["strict_t"][m]; a = z["strict_aM"][m]; oot = z["strict_oot"][m].astype(bool)
            sg = robust_std(a[oot])
            aT = z["strict_aT"][m]
            rows.append(pd.DataFrame(dict(time_btjd_tdb=t, time_bjd_tdb=t + 2457000.0, sector=s, epoch=int(ep), window_id=wid,
                                          cadence_s=ex, dflux_e_s=a, flux_norm=1 + a / F_M3V, flux_err=sg / F_M3V, beta_transit=beta,
                                          oot=oot.astype(int), T_M3V_e_s=aT[:, 0], T_pair_e_s=aT[:, 1], T_Mira_e_s=aT[:, 2],
                                          bkg_e_s_px=z["strict_B"][m])))
            wid += 1
    d = pd.concat(rows, ignore_index=True)
    d.to_csv(f"{CL}/lc_M3V_REF_clean.csv", index=False, float_format="%.7f")
    rows = []
    for s in SECTORS:
        z = np.load(f"{CL}/data/reftile_all_s{s}.npz")
        rows.append(pd.DataFrame(dict(time_btjd_tdb=z["t"], time_bjd_tdb=z["t"] + 2457000.0, sector=s, tile_id=z["win"],
                                      dflux_e_s=z["aM"], flux_norm=1 + z["aM"] / F_M3V, flux_err=z["sig"] / F_M3V,
                                      in_transit_window=intr(z["t"]).astype(int))))
    pd.concat(rows, ignore_index=True).to_csv(f"{CL}/lc_M3V_REF_tiled.csv", index=False, float_format="%.7f")
    return len(d)


if __name__ == "__main__":
    n1 = inv(); n2 = ref()
    print("rows", n1, n2)
