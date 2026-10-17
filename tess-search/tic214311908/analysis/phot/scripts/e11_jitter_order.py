import sys; sys.path.insert(0, ".")
from p05_phot import *
from p06_lc import detrend, FREF, BKG_EXC_MAX
from tlib import ephem_frac, fit_single
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"]
for s in (103, 66, 104):
    S = setup(s)
    for jo in (1, 2):
        R = photometry(s, S=S, verbose=False, jitter_order=jo)
        t = R["t"]; chiB = R["FIT"]["B"]["chi2"]; bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
        ok = R["qok"] & (bexc < BKG_EXC_MAX) & R["jit_ok"] & np.isfinite(chiB) & (chiB < 3 * np.nanmedian(chiB[R["q0"]]))
        mask = np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1
        msg = []
        for mk, nm in (("B", "M3V"), ("A", "M3V"), ("A", "pair")):
            names = R["FIT"][mk]["names"]; raw = 1 + R["FIT"][mk]["F"][:, names.index(nm)] / FREF[nm]
            d = raw / detrend(t, raw, ok, mask); m = ok & np.isfinite(d)
            sd = robust_std(d[m & ~mask])
            eps = []
            for k in np.unique(np.round((t[m] - T0) / P)):
                tc0 = T0 + k * P; w = m & (np.abs(t - tc0) < 0.3)
                if np.sum(m & (np.abs(t - tc0) < T14 / 2)) < 10: continue
                r = fit_single(t[w], d[w], np.full(w.sum(), sd), R["exptime"], tc0, T14, FL)
                if r: eps.append(f"{int(k)}:{r['depth']*1e3:.1f}")
            msg.append(f"{mk}:{nm} rstd {sd*1e3:.1f} [{' '.join(eps)}]")
        print(f"S{s} jitter_order {jo}: " + " | ".join(msg), flush=True)
