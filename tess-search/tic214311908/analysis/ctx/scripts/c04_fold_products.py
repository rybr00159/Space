"""Fold every downloaded light-curve product (QLP, TGLC, GSFC-ELEANOR-LITE, TARS, TEQUILA) on the ephemeris and
measure: transit depth (trapezoid template, T14 1.63 h, integrated over the exposure), formal error from the robust
out-of-transit scatter, an EMPIRICAL error from depths at 400 random phases, odd/even depths, phase-0.5 depth and the
absolute flux lost (e-/s) when the product is in e-/s.

Output: out/fold_products.csv (one row per product/star/sector/column plus 'ALL' rows combining sectors),
        out/fold_events_<product>_<id>.csv  per-transit depths for the main columns,
        out/fold_arrays.npz (detrended folded data for plotting).
Run:  cd ctx && OMP_NUM_THREADS=2 python scripts/c04_fold_products.py
"""
import os, sys, glob, re, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
from astropy.io import fits
import common as C

LC = "lc"
EXPT = {12: 1800.0, 39: 600.0}
def exptime(sec):
    return EXPT.get(sec, 200.0) / 86400.0

NAME = {"214311906": "Mira", "214311908": "M3V", "214311910": "pairA", "1515709367": "pairA_dup", "1515709322": "pairB",
        "5925373268528684288": "Mira", "5925373272850784768": "M3V", "5925373272850782208": "pairA",
        "5925373272827253632": "pairB", "5925373204131295872": "ref1", "5925373204131294336": "ref2",
        "s0012-cam1-ccd2-x1515-y0269": "Mira_pixel"}


def load(path):
    """Return list of dicts: product, id, sector, column, t, f(raw, e-/s or normalised), quality-ok mask, unit."""
    b = os.path.basename(path); out = []
    h = fits.open(path); d = h[1].data
    if b.startswith("hlsp_qlp"):
        sec = int(b.split("_s")[1][:4]); tic = str(int(b.split("-")[1][:16]))
        t = d["TIME"] + C.QLP_SHIFT.get(sec, 0.0); ok = (d["QUALITY"] == 0)
        det = "DET_FLUX" if "DET_FLUX" in d.columns.names else "KSPSAP_FLUX"
        for col in ("SAP_FLUX", det):
            out.append(dict(product="QLP", id=tic, sector=sec, column="SAP_FLUX" if col == "SAP_FLUX" else "DET_FLUX", t=t,
                            f=np.array(d[col], float), ok=ok, unit="norm"))   # QLP SAP_FLUX is normalised (median ~1)
    elif b.startswith("hlsp_tglc"):
        sec = int(b.split("-s")[1][:4]); gid = b.split("gaiaid-")[1].split("-")[0]
        t = d["time"]; ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0)
        for col in ("psf_flux", "aperture_flux"):
            out.append(dict(product="TGLC", id=gid, sector=sec, column=col, t=t, f=np.array(d[col], float), ok=ok, unit="e/s"))
    elif b.startswith("hlsp_gsfc-eleanor-lite") and "214311908" not in b:
        return out   # all five ELEANOR-LITE files are byte-identical blends (same aperture L_180); keep one copy
    elif b.startswith("hlsp_gsfc-eleanor-lite"):
        m_ = re.search(r"_s(\d{4})-(\d{16})_", b); sec = int(m_.group(1)); tic = str(int(m_.group(2)))
        t = d["TIME"]; ok = d["QUALITY"] == 0
        for col in ("RAW_FLUX", "CORR_FLUX", "PCA_FLUX"):
            out.append(dict(product="ELEANOR-LITE", id=tic, sector=sec, column=col, t=t, f=np.array(d[col], float), ok=ok, unit="e/s"))
    elif b.startswith("hlsp_tars"):
        sec = int(b.split("_s")[1][:4]); tic = str(int(b.split("-")[1][:16]))
        out.append(dict(product="TARS", id=tic, sector=sec, column="FLUX", t=d["TIME"], f=np.array(d["FLUX"], float),
                        ok=np.ones(len(d), bool), unit="norm"))
    elif b.startswith("hlsp_tequila"):
        out.append(dict(product="TEQUILA", id="s0012-cam1-ccd2-x1515-y0269", sector=12, column="FLUX",
                        t=d["TIME"] - 2457000.0, f=np.array(d["FLUX"], float), ok=d["SAP_QUALITY"] == 0, unit="e/s"))
    return out


def measure(t, f, s, ex, tag):
    tm0 = C.template(t, ex); tm5 = C.template(t, ex, phase0=0.5)
    d0, e0 = C.fit_depth(f, s, tm0); d5, e5 = C.fit_depth(f, s, tm5)
    ep = C.epoch(t)
    res = dict(n_in=int(np.sum(tm0 > 0.5)), depth=d0, err=e0, half=d5, half_err=e5)
    for par, lab in ((0, "even"), (1, "odd")):
        m = (ep % 2 == par)
        res[lab], res[lab + "_err"] = C.fit_depth(f[m], s[m], tm0[m]) if np.sum(tm0[m] > 0.5) > 0 else (np.nan, np.nan)
    rd = C.random_phase_depths(t, f, s, ex, n=400, rng=np.random.default_rng(abs(hash(tag)) % 2**32))
    res["err_emp"] = 1.4826 * np.median(np.abs(rd - np.median(rd))) if len(rd) > 20 else np.nan
    res["err_emp_std"] = np.std(rd) if len(rd) > 20 else np.nan
    res["rand_mean"] = np.mean(rd) if len(rd) > 20 else np.nan
    res["p_rand_ge"] = (np.sum(rd >= d0) + 1) / (len(rd) + 1) if len(rd) > 20 else np.nan
    res["snr"] = d0 / e0; res["snr_emp"] = d0 / max(res["err_emp_std"], e0) if np.isfinite(res["err_emp_std"]) else np.nan
    return res, rd


def main():
    files = sorted(glob.glob(f"{LC}/*.fits") + glob.glob(f"{LC}/tglc/*.fits"))
    rows, ev_rows, arrays = [], [], {}
    store = {}   # (product, id, column) -> list of (t,f,s,ex)
    for fn in files:
        for L in load(fn):
            t, f, ok = L["t"], L["f"], L["ok"] & np.isfinite(L["t"]) & np.isfinite(L["f"])
            t, f = t[ok], f[ok]          # NB: no f>0 cut (TGLC fluxes of faint stars legitimately cross zero)
            if len(t) < 50:
                continue
            med = np.nanmedian(f)
            mask = C.in_window(t, half=1.0 * C.T14) | C.in_window(t, phase0=0.5, half=1.0 * C.T14)
            trend = C.detrend_additive(t, f, mask, window=1.0)
            r = f - trend                                   # residual in the product's own units
            sig_abs = C.robust_std(r[~mask])
            good = ~((~mask) & (np.abs(r) > 5 * sig_abs)) & np.isfinite(r)
            t, r, trend = t[good], r[good], trend[good]
            # depth relative to the product's median level (for TGLC that level is Gaia-predicted, not measured)
            fd = 1 + r / med; sig = sig_abs / abs(med)
            trend_med = trend
            s = np.full(len(t), sig)
            ex = exptime(L["sector"])
            tag = f"{L['product']}_{L['id']}_{L['sector']}_{L['column']}"
            r_, rd = measure(t, fd, s, ex, tag); r = r_
            intr = C.template(t, ex) > 0.5
            lvl = med if L["unit"] == "e/s" else np.nan
            rows.append(dict(product=L["product"], id=L["id"], star=NAME.get(L["id"], "?"), sector=L["sector"], column=L["column"],
                             median_flux=med if L["unit"] == "e/s" else np.nan, level_in_transit=np.median(trend[intr]) if intr.any() else np.nan,
                             rms_ppm=sig * 1e6,
                             dflux_es=r["depth"] * lvl if np.isfinite(lvl) else np.nan,
                             dflux_es_err=r["err"] * lvl if np.isfinite(lvl) else np.nan,
                             dflux_es_err_emp=r["err_emp_std"] * lvl if np.isfinite(lvl) else np.nan, **r))
            store.setdefault((L["product"], L["id"], L["column"]), []).append((t, fd, s, np.full(len(t), ex), np.full(len(t), L["sector"])))
            arrays[tag] = np.vstack([t, fd])
            # per-event depths
            ep = C.epoch(t); tm = C.template(t, ex)
            for e in np.unique(ep[tm > 0.5]):
                m = (ep == e) & (np.abs(t - (C.T0 + e * C.P)) < 3 * C.T14)
                if np.sum(tm[m] > 0.5) >= 1:
                    d, de = C.fit_depth(fd[m], s[m], tm[m])
                    ev_rows.append(dict(product=L["product"], id=L["id"], star=NAME.get(L["id"], "?"), sector=L["sector"],
                                        column=L["column"], epoch=int(e), n_in=int(np.sum(tm[m] > 0.5)), depth=d, err=de))
            print(f"{tag:55s} n_in={r['n_in']:4d} depth={r['depth']*1e6:8.0f} +- {r['err']*1e6:6.0f} (emp {r['err_emp_std']*1e6:6.0f}) "
                  f"S/N {r['snr']:5.1f}/{r['snr_emp']:5.1f}  half {r['half']*1e6:7.0f}+-{r['half_err']*1e6:5.0f}  "
                  f"odd {r['odd']*1e6:7.0f}+-{r['odd_err']*1e6:5.0f} even {r['even']*1e6:7.0f}+-{r['even_err']*1e6:5.0f}  dF={rows[-1]['dflux_es']:.1f} e/s", flush=True)
    # combined over sectors (each sector's residual scatter as its weight)
    for (prod, idd, col), lst in store.items():
        if len(lst) < 2:
            continue
        t = np.concatenate([x[0] for x in lst]); f = np.concatenate([x[1] for x in lst]); s = np.concatenate([x[2] for x in lst])
        exa = np.concatenate([x[3] for x in lst])
        tm0 = np.concatenate([C.template(x[0], x[3][0]) for x in lst]); tm5 = np.concatenate([C.template(x[0], x[3][0], phase0=0.5) for x in lst])
        d0, e0 = C.fit_depth(f, s, tm0); d5, e5 = C.fit_depth(f, s, tm5)
        ep = C.epoch(t); r = dict(n_in=int(np.sum(tm0 > 0.5)), depth=d0, err=e0, half=d5, half_err=e5)
        for par, lab in ((0, "even"), (1, "odd")):
            m = ep % 2 == par; r[lab], r[lab + "_err"] = C.fit_depth(f[m], s[m], tm0[m])
        rng = np.random.default_rng(7); rd = []
        while len(rd) < 400:
            ph0 = rng.uniform(0, 1)
            if min(abs(ph0), abs(ph0 - 0.5), abs(1 - ph0)) < 2.5 * C.T14 / C.P:
                continue
            tm = np.concatenate([C.template(x[0], x[3][0], phase0=ph0) for x in lst]); rd.append(C.fit_depth(f, s, tm)[0])
        rd = np.array(rd)
        r["err_emp"] = 1.4826 * np.median(np.abs(rd - np.median(rd))); r["err_emp_std"] = np.std(rd); r["rand_mean"] = rd.mean()
        r["p_rand_ge"] = (np.sum(rd >= d0) + 1) / (len(rd) + 1); r["snr"] = d0 / e0; r["snr_emp"] = d0 / max(r["err_emp_std"], e0)
        rows.append(dict(product=prod, id=idd, star=NAME.get(idd, "?"), sector="ALL", column=col, median_flux=np.nan,
                         level_in_transit=np.nan, rms_ppm=np.nan, dflux_es=np.nan, dflux_es_err=np.nan, **r))
        print(f"ALL {prod} {idd} {col}: depth {d0*1e6:.0f} +- {e0*1e6:.0f} (emp {r['err_emp_std']*1e6:.0f}) S/N {r['snr']:.1f}/{r['snr_emp']:.1f}; "
              f"half {d5*1e6:.0f}+-{e5*1e6:.0f}; odd {r['odd']*1e6:.0f}+-{r['odd_err']*1e6:.0f} even {r['even']*1e6:.0f}+-{r['even_err']*1e6:.0f}", flush=True)
    pd.DataFrame(rows).to_csv("out/fold_products.csv", index=False)
    pd.DataFrame(ev_rows).to_csv("out/fold_events.csv", index=False)
    np.savez_compressed("out/fold_arrays.npz", **arrays)


if __name__ == "__main__":
    main()
