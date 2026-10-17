"""Light curves from the deblended photometry (p05): target-barycentred times, quality flags, normalisation,
transit-masked detrending; written per sector (npz, all stars and models) and as combined CSVs.

Models (p05): B = Mira + M3V + plane (pair not fitted; primary for the M3V); Bp = Mira + pair + plane (pair light
curve; the M3V is not fitted); A = Mira + M3V + pair + plane (all free; noisier, templates correlated).
Normalisation (reference fluxes, e-/s): M3V 940 (TIC 8.2 Tmag 13.0135); pair A+B 894 (TIC Tmag 13.3749 for A plus
Gaia-based T 14.56 for B); Mira = its mean aperture flux in the regression set (the Mira template is normalised to the
r = 2.5 px aperture) minus the M dwarfs' share of that aperture. Relative depths scale with these references.
Quality: quality_ok = finite & QUALITY in {0, 2048 (straylight)} & background excess < 60 e/s/px over the sector
median & pointing offset within 6 robust sigma & per-cadence reduced chi2 of model B below 3x its median.
Detrending: wotan biweight (window 0.75 d, break tolerance 0.5 d) on the normalised flux with the transit windows of
the reference ephemeris (+-0.1 d = +-1.5 T14) masked; the trend is interpolated across the windows.
Flare flag: detrended flux > 1 + 4 robust sigma (and its two neighbours).
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from wotan import flatten

FREF = dict(M3V=tmag_flux(13.0135), pair=tmag_flux(13.3749) + tmag_flux(14.56))
WIN = 0.75
BKG_EXC_MAX = 60.0     # e/s/px above the sector median background (noise and offsets grow beyond, see e09_bkgcut.py)


def detrend(t, f, ok, mask_tr, window=WIN):
    m = ok & np.isfinite(f) & ~mask_tr
    tr = np.full(len(t), np.nan)
    if m.sum() < 20: return tr
    trm = flatten(t[m], f[m], method="biweight", window_length=window, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
    k = np.isfinite(trm)
    fin = np.isfinite(t)
    tr[fin] = np.interp(t[fin], t[m][k], trm[k])
    # do not extrapolate across gaps longer than 0.5 d from the nearest trend point
    tk = t[m][k]
    idx = np.clip(np.searchsorted(tk, t[fin]), 1, len(tk) - 1)
    dist = np.minimum(np.abs(t[fin] - tk[idx - 1]), np.abs(t[fin] - tk[idx]))
    tr_f = tr[fin]; tr_f[dist > 0.5] = np.nan; tr[fin] = tr_f
    return tr


def in_transit_windows(t, P=P_EPH, T0=T0_EPH, half=0.1):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    return np.abs(ph * P) < half


def process(sector, z=None, J=None, tag="", verbose=True):
    if z is None:
        z = dict(np.load(f"{ROOT}/data/phot_s{sector}{tag}.npz"))
        J = json.load(open(f"{ROOT}/out/phot_s{sector}{tag}.json"))
    TC = json.load(open(f"{ROOT}/out/timecorr.json"))
    t = z["t"]; fin = np.isfinite(t)
    add = np.zeros(len(t))
    tc_ = z["timecorr"]; m = fin & np.isfinite(tc_)
    add[m] = roemer_target(sector, t[m], tc_[m])
    t_tdb = t + add / 86400.0
    chiB = z["chi2_B"]
    bexc = z["bkg"] - np.median(z["bkg"][z["q0"]])
    base_ok = z["qok"] & (bexc < BKG_EXC_MAX) & z["jit_ok"] & np.isfinite(chiB)
    base_ok &= chiB < 3 * np.nanmedian(chiB[z["q0"]])
    mask_tr = in_transit_windows(t_tdb)
    # Mira reference flux: mean aperture flux minus the M dwarfs' expected share of the r=2.5 px aperture
    fref = dict(FREF, Mira=float(J["mu"]["F"] - 0.9 * (FREF["M3V"] + FREF["pair"])))
    out = dict(t_tdb=t_tdb, t_tess=t, q=z["q"], cadn=z["cadn"], quality_ok=base_ok, bkg=z["bkg"], dx=z["dx"], dy=z["dy"], fin=fin,
               intransit_window=mask_tr)
    stats = {}
    for mk, stars in (("B", ("M3V", "Mira")), ("Bp", ("pair",)), ("A", ("M3V", "pair", "Mira"))):
        names = list(z[f"names_{mk}"])
        for nm in stars:
            dF = z[f"F_{mk}"][:, names.index(nm)]
            raw = 1.0 + dF / fref[nm]
            tr = detrend(t_tdb, raw, base_ok, mask_tr)
            det = raw / tr
            ok = base_ok & np.isfinite(det)
            sd = robust_std(det[ok & ~mask_tr])
            flare = np.zeros(len(t), bool)
            hi = ok & (det > 1 + 4 * sd)
            for sh in (-1, 0, 1): flare |= np.roll(hi, sh)
            key = f"{nm}_{mk}"
            out[f"dF_{key}"] = dF; out[f"raw_{key}"] = raw; out[f"trend_{key}"] = tr; out[f"det_{key}"] = det
            out[f"flare_{key}"] = flare
            # per-cadence error: sector robust sigma scaled by the local chi2 of the fit
            chi = z[f"chi2_{mk}"]; scale = np.sqrt(np.clip(chi / np.nanmedian(chi[ok]), 1.0, None))
            out[f"err_{key}"] = sd * scale
            stats[key] = dict(rstd=float(sd), rstd_eps=float(sd * fref[nm]), n_ok=int(ok.sum()), n_flare=int(flare.sum()), fref=float(fref[nm]))
    np.savez_compressed(f"{ROOT}/data/lc_s{sector}{tag}.npz", sector=sector, exptime_s=J["exptime"] * 86400, **out)
    if verbose:
        print(f"S{sector}: n_ok {int(base_ok.sum())}/{int(fin.sum())}; time corr {np.nanmean(add):+.1f} s; " +
              "; ".join(f"{k} rstd {v['rstd']*1e3:.1f} ppt ({v['rstd_eps']:.1f} e/s), flares {v['n_flare']}" for k, v in stats.items()), flush=True)
    return out, stats, fref


def write_csv(all_out):
    rows = []
    for s, (o, st, fref) in all_out.items():
        cad = round(float(np.nanmedian(np.diff(o["t_tess"][o["fin"]]))) * 86400)
        sel = o["fin"]
        rows.append(pd.DataFrame(dict(
            time_btjd_tdb=o["t_tdb"][sel], sector=s, cadence_s=cad, quality_ok=o["quality_ok"][sel].astype(int),
            flux_norm=o["det_M3V_B"][sel], flux_err=o["err_M3V_B"][sel], flux_raw_norm=o["raw_M3V_B"][sel], trend=o["trend_M3V_B"][sel],
            dflux_e_s=o["dF_M3V_B"][sel], flare_flag=o["flare_M3V_B"][sel].astype(int),
            flux_norm_modelA=o["det_M3V_A"][sel], flux_err_modelA=o["err_M3V_A"][sel], flux_raw_norm_modelA=o["raw_M3V_A"][sel],
            in_transit_window=o["intransit_window"][sel].astype(int), tess_quality=o["q"][sel], bkg_e_s_px=o["bkg"][sel],
            jitter_dx_px=o["dx"][sel], jitter_dy_px=o["dy"][sel], time_btjd_tesscut=o["t_tess"][sel])))
    d = pd.concat(rows, ignore_index=True)
    d.to_csv(f"{ROOT}/lc_M3V_all.csv", index=False, float_format="%.7f")
    for star, models in (("pair", ("Bp", "A")), ("Mira", ("B", "A"))):
        rows = []
        for s, (o, st, fref) in all_out.items():
            cad = round(float(np.nanmedian(np.diff(o["t_tess"][o["fin"]]))) * 86400); sel = o["fin"]
            m1, m2 = models
            rows.append(pd.DataFrame(dict(
                time_btjd_tdb=o["t_tdb"][sel], sector=s, cadence_s=cad, quality_ok=o["quality_ok"][sel].astype(int),
                flux_norm=o[f"det_{star}_{m1}"][sel], flux_err=o[f"err_{star}_{m1}"][sel], flux_raw_norm=o[f"raw_{star}_{m1}"][sel],
                dflux_e_s=o[f"dF_{star}_{m1}"][sel], flux_norm_modelA=o[f"det_{star}_{m2}"][sel], flux_raw_norm_modelA=o[f"raw_{star}_{m2}"][sel],
                dflux_e_s_modelA=o[f"dF_{star}_{m2}"][sel], in_transit_window=o["intransit_window"][sel].astype(int))))
        pd.concat(rows, ignore_index=True).to_csv(f"{ROOT}/lc_{star}_all.csv", index=False, float_format="%.7f")
    save_json({str(s): dict(stats=st, fref=fref) for s, (o, st, fref) in all_out.items()}, f"{ROOT}/out/lc_stats.json")


if __name__ == "__main__":
    allo = {}
    for s in SECTORS:
        allo[s] = process(s)
    write_csv(allo)
