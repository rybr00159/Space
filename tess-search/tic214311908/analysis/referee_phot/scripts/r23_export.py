"""Referee r23: export the referee local-window M3V photometry (r05, strict selection, one window per transit) as a CSV
for downstream cross-checks of the investigator's lc_M3V_all.csv.
Columns: time_btjd_tdb (BJD_TDB - 2457000, barycentred for TIC 214311908 with astropy, r02), sector, epoch, window_id,
cadence_s, dflux_e_s (M3V-only PRF fit, e-/s, per-pixel baseline already removed so out-of-transit ~0),
flux_norm (= 1 + dflux_e_s / 979, 979 e-/s = M3V flux in PRF units from Gaia XP, r09/r17), flux_err (robust OOT
scatter of the window / 979), oot (1 = outside +-0.07 d of the window centre), T_M3V_e_s, T_pair_e_s, T_Mira_e_s
(three-star fit), bkg_e_s_px (referee faint-pixel background). Windows are +-0.35 d (S39 +-0.4, S12 +-0.6) around
T0 3840.91414 + k 5.6915087 (investigator ephemeris, used only to place windows).
Output: lc_M3V_referee_windows.csv"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
F = 979.0
rows = []
for s in SECTORS:
    z = np.load(f"{R}/data/win_s{s}.npz", allow_pickle=True)
    win = z["strict_win"]; ep = z["strict_epoch"]; t = z["strict_t"]; aM = z["strict_aM"]; aT = z["strict_aT"]; oot = z["strict_oot"]; B = z["strict_B"]
    sg = np.zeros(len(t))
    for i in np.unique(win):
        m = win == i; sg[m] = rstd(aM[m & oot])
    rows.append(pd.DataFrame(dict(time_btjd_tdb=t, sector=s, epoch=ep[win], window_id=win, cadence_s=np.round(float(z["exptime"]) * 86400),
                                  dflux_e_s=aM, flux_norm=1 + aM / F, flux_err=sg / F, oot=oot.astype(int),
                                  T_M3V_e_s=aT[:, 0], T_pair_e_s=aT[:, 1], T_Mira_e_s=aT[:, 2], bkg_e_s_px=B)))
D = pd.concat(rows, ignore_index=True)
D.to_csv(f"{R}/lc_M3V_referee_windows.csv", index=False, float_format="%.7f")
print(len(D), "rows;", D.groupby("sector").epoch.unique().to_dict())
