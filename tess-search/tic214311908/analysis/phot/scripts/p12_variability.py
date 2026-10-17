"""Out-of-transit variability of the M3V (rotation, flares) with control stars processed identically.

Per sector: the p05 pipeline is rerun with up to four control stars (Gaia field stars T 12-14.2, > 4 px from the Mira,
isolated, two of them as red as possible) photometered by the same difference-image PRF fit. For the M3V (model B),
the pair (model Bp) and the controls, the undetrended relative flux (transit windows removed) is regressed against
background excess and pointing offsets (and, for the M3V and pair, the Mira's flux change), then a Lomb-Scargle
periodogram (0.2-13 d) is computed. A period is credible only if it recurs across sectors and is absent from the
controls. The TARS catalogue periods for TIC 214311908 (header PER1..5) are compared.
Flares: detrended M3V flux, events with >= 2 consecutive cadences > 3 sigma or one > 5 sigma; per-sector detection
threshold (5 sigma per cadence) reported.
Outputs: out/variability.json, figs/variability.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p05_phot import setup, photometry
from p06_lc import FREF, BKG_EXC_MAX
from astropy.timeseries import LombScargle
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0 = TR["P"], TR["T0"]


def pick_controls(C, W, g_all):
    g = g_all.copy(); g["x"] += W["global_shift"][0]; g["y"] += W["global_shift"][1]
    mira = g[g.source_id == IDS["Mira"]].iloc[0]
    cand = []
    for r in g.itertuples():
        if not (12.0 <= r.T <= 14.2 and np.isfinite(r.bp_rp)): continue
        if not (2.5 < r.x < C.nx - 3.5 and 2.5 < r.y < C.ny - 3.5): continue
        if np.hypot(r.x - mira.x, r.y - mira.y) < 4.5: continue
        d = np.hypot(g.x - r.x, g.y - r.y)
        if np.any((d > 0) & (d < 1.5) & (g["T"] < r.T + 2.0)): continue
        cand.append((r.source_id, r.x, r.y, r.T, r.bp_rp))
    cand = sorted(cand, key=lambda c: -c[4])
    sel = cand[:2] + [c for c in sorted(cand[2:], key=lambda c: abs(c[3] - 13.0))][:2]
    return sel


def ls(t, f, pmin=0.2, pmax=13.0):
    k = np.isfinite(f) & np.isfinite(t)
    L = LombScargle(t[k], f[k])
    fr, pw = L.autopower(minimum_frequency=1 / pmax, maximum_frequency=1 / pmin, samples_per_peak=10)
    j = int(np.argmax(pw))
    amp = np.sqrt(4 * pw[j] * np.var(f[k]))   # approx semi-amplitude for a sinusoid
    return float(1 / fr[j]), float(pw[j]), float(L.false_alarm_probability(pw[j])), float(amp), fr, pw


def main():
    out = {}
    fig, axs = plt.subplots(len(SECTORS), 2, figsize=(16, 3.2 * len(SECTORS)))
    for row, s in enumerate(SECTORS):
        S = setup(s); C, W, J, g, prf, rBA = S
        ctr = pick_controls(C, W, C.gaia(18.5))
        R = photometry(s, S=S, verbose=False, controls={str(c[0]): (c[1], c[2]) for c in ctr})
        t = R["t"]; bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
        chiB = R["FIT"]["B"]["chi2"]
        ok = R["qok"] & (bexc < BKG_EXC_MAX) & R["jit_ok"] & np.isfinite(chiB) & (chiB < 3 * np.nanmedian(chiB[R["q0"]]))
        ok &= ~(np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1)
        nB = R["FIT"]["B"]["names"]; nBp = R["FIT"]["Bp"]["names"]
        mira = R["FIT"]["B"]["F"][:, nB.index("Mira")]
        series = {"M3V": R["FIT"]["B"]["F"][:, nB.index("M3V")] / FREF["M3V"], "pair": R["FIT"]["Bp"]["F"][:, nBp.index("pair")] / FREF["pair"]}
        for c in ctr: series[f"ctrl_{c[0]}_T{c[3]:.1f}_BR{c[4]:.1f}"] = R["CTRL"][str(c[0])] / tmag_flux(c[3])
        res = {}
        # bin to 1 h for speed
        for nm, f in series.items():
            m = ok & np.isfinite(f)
            regs = [np.ones(m.sum()), bexc[m], R["dx"][m] - np.nanmean(R["dx"][m]), R["dy"][m] - np.nanmean(R["dy"][m])]
            if nm in ("M3V", "pair"): regs.append(mira[m] - np.nanmean(mira[m]))
            A = np.column_stack(regs); c_, *_ = np.linalg.lstsq(A, f[m], rcond=None)
            r = f[m] - A[:, 1:] @ c_[1:]
            tb = np.floor((t[m] - t[m][0]) * 24); u, inv = np.unique(tb, return_inverse=True)
            rb = np.bincount(inv, r) / np.bincount(inv); tbm = np.bincount(inv, t[m]) / np.bincount(inv)
            per, pw, fap, amp, fr, pws = ls(tbm, rb)
            res[nm] = dict(period_d=per, power=pw, fap=fap, semi_amp_ppt=amp * 1e3, rms_1h_ppt=float(np.std(rb) * 1e3))
            if nm == "M3V" or nm.startswith("ctrl") and len([k for k in res if k.startswith("ctrl")]) <= 2:
                axs[row, 1].plot(1 / fr, pws, lw=0.8, label=nm[:22])
            if nm in ("M3V",) or nm.startswith("ctrl"):
                axs[row, 0].plot(tbm, rb * 1e3 + (0 if nm == "M3V" else 30 * (1 + list(series).index(nm))), ".", ms=2, label=nm[:22])
        axs[row, 1].set_xscale("log"); axs[row, 1].set_xlabel("period (d)"); axs[row, 1].legend(fontsize=7); axs[row, 0].set_title(f"S{s} 1-h binned, systematics-regressed (ppt; offsets)")
        out[s] = res
        print(f"S{s}: " + "; ".join(f"{k[:24]} P {v['period_d']:.2f} d pow {v['power']:.3f} A {v['semi_amp_ppt']:.1f} ppt rms {v['rms_1h_ppt']:.1f}" for k, v in res.items()), flush=True)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/variability.png", dpi=70); plt.close(fig)
    # flares
    M = pd.read_csv(f"{ROOT}/lc_M3V_all.csv"); fl = []; thr = {}
    for s in SECTORS:
        d = M[(M.sector == s) & (M.quality_ok == 1)].reset_index(drop=True)
        f = d.flux_norm.values - 1; sd = robust_std(f); thr[s] = dict(sigma_ppt=sd * 1e3, five_sigma_ppt=5 * sd * 1e3, days=float(len(d) * d.cadence_s[0] / 86400))
        hi3 = f > 3 * sd; hi5 = f > 5 * sd
        i = 0
        while i < len(f):
            if hi5[i] or (i + 1 < len(f) and hi3[i] and hi3[i + 1] and d.time_btjd_tdb[i + 1] - d.time_btjd_tdb[i] < 2 * d.cadence_s[0] / 86400):
                j = i
                while j + 1 < len(f) and hi3[j + 1]: j += 1
                fl.append(dict(sector=s, t_start=float(d.time_btjd_tdb[i]), n_cad=int(j - i + 1), peak_ppt=float(f[i:j + 1].max() * 1e3), peak_sigma=float(f[i:j + 1].max() / sd)))
                i = j + 1
            else: i += 1
    out["flares"] = fl; out["flare_thresholds"] = thr
    tot_days = sum(v["days"] for v in thr.values())
    # expected false positives from Gaussian noise
    print(f"flare candidates: {len(fl)} in {tot_days:.1f} d of good data: " + ", ".join(f"S{f_['sector']} {f_['t_start']:.3f} ({f_['n_cad']} cad, {f_['peak_ppt']:.0f} ppt, {f_['peak_sigma']:.1f} sig)" for f_ in fl))
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/variability.json")


if __name__ == "__main__":
    main()
