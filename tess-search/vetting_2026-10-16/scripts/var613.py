"""Shape and period of the TIC 61332759 variability: periodogram and Fourier fits at P and 2P per sector (QLP SAP)."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "61332759"]
import numpy as np
from astropy.timeseries import LombScargle
from wotan import flatten
import stage2
from lcan import read, bin10
P0 = 0.740299
p = stage2.products(61332759); out = {}
for _, r in p[p.prov.isin(["QLP", "TESS-SPOC"])].iterrows():
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), r.prov)
    c = "SAP_FLUX"
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[c]) & (cols[c] > 0)
    tt, ff = bin10(t[ok], cols[c][ok] / np.nanmedian(cols[c][ok]), r.exptime)
    fl = flatten(tt, ff, method="biweight", window_length=4.0, break_tolerance=0.4)
    m = np.isfinite(fl); tt, fl = tt[m], fl[m]
    ls = LombScargle(tt, fl); fr, pw = ls.autopower(minimum_frequency=1 / 10, maximum_frequency=1 / 0.1, samples_per_peak=20)
    top = np.argsort(pw)[::-1]; peaks = []
    for i in top:
        if all(abs(fr[i] - fr[j]) > 0.02 for j in peaks): peaks.append(i)
        if len(peaks) == 3: break
    # Fourier fit at 2P with 4 harmonics: odd harmonics measure differences between alternate cycles
    P2 = 2 * P0; ph = 2 * np.pi * (tt - 3148.7857) / P2
    A = np.column_stack([np.ones(len(tt))] + [f(k * ph) for k in range(1, 5) for f in (np.cos, np.sin)])
    cf, *_ = np.linalg.lstsq(A, fl - 1, rcond=None); res = fl - 1 - A @ cf
    cov = np.linalg.pinv(A.T @ A) * np.var(res); e = np.sqrt(np.diag(cov))
    amp = [np.hypot(cf[2 * k - 1], cf[2 * k]) * 1e6 for k in range(1, 5)]
    err = [np.hypot(e[2 * k - 1], e[2 * k]) * 1e6 for k in range(1, 5)]
    key = f"S{int(r.sector)} {r.prov}"
    out[key] = dict(peaks=[(1 / fr[i], float(pw[i])) for i in peaks], amp_2P_harm=amp, err=err)
    print(f"{key:15s} LS peaks (d): " + ", ".join(f"{1/fr[i]:.4f} ({pw[i]:.3f})" for i in peaks) +
          "  | amplitudes at 2P harmonics k=1..4 (ppm): " + ", ".join(f"{a:.0f}±{x:.0f}" for a, x in zip(amp, err)))
json.dump(out, open("out/61332759/variability.json", "w"), indent=1)
