"""R04b: which QLP depth estimator is unbiased? Inject a 15000-ppm trapezoid (T14 1.63 h, flat 70%) into the M3V QLP
SAP/DET light curves (S39-S104, real noise) at phase offsets that avoid the real transit and phase 0.5, then recover it
with (a) the referee's local-quadratic joint fit (rlib/r04) and (b) a re-implementation of ctx's estimator (1-d
biweight trend with +-T14 masks at phase 0 and 0.5 around the TRUE ephemeris, then a global template fit) applied
exactly as ctx did, i.e. WITHOUT masking the injected (unknown) transit... and (c) ctx's estimator with the injected
transit masked (as it would be for the real transit). Output: out/r04b_depth_injection.json
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from wotan import flatten
import rlib as R
import r04_qlp_fold as Q

def ctx_est(t, f, ex, ph_inj, mask_inj):
    ph = ((t - R.T0) / R.P - ph_inj + 0.5) % 1 - 0.5
    m0 = np.abs((((t - R.T0) / R.P + 0.5) % 1 - 0.5) * R.P) < R.T14
    m5 = np.abs((((t - R.T0) / R.P) % 1 - 0.5) * R.P) < R.T14
    mi = np.abs(ph * R.P) < R.T14
    mask = m0 | m5 | (mi if mask_inj else False)
    tr = flatten(t[~mask], f[~mask], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    k = np.isfinite(tr); trend = np.interp(t, t[~mask][k], tr[k]); r = f / trend
    s = 1.4826 * np.median(np.abs(r[~mask] - np.median(r[~mask])))
    good = ~((~mask) & (np.abs(r - 1) > 5 * s))
    tm = R.shape(ph[good] * R.P, 0.0, ex)
    return np.sum(tm * (1 - r[good])) / np.sum(tm * tm), s

res = {}
for col in ("DET_FLUX", "SAP_FLUX"):
    data = []
    for sec in Q.SECT:
        L = Q.load(Q.TIC["M3V"], sec, col); data.append((L[0], L[1], R.exptime(sec)))
    rec = {"referee": [], "ctx_unmasked": [], "ctx_masked": []}
    for phi in np.arange(0.12, 0.88, 0.04):
        if abs(phi - 0.5) < 0.08: continue
        dinj = []
        for (t, f, ex) in data:
            ph = ((t - R.T0) / R.P - phi + 0.5) % 1 - 0.5
            dinj.append((t, f * (1 - 0.015 * R.shape(ph * R.P, 0.0, ex)), ex))
        # referee estimator at the injected phase (offset = phi*P) minus the same without injection
        d_inj = Q.measure(dinj, offset=phi * R.P)[0]; d_0 = Q.measure(data, offset=phi * R.P)[0]
        rec["referee"].append(d_inj - d_0)
        for key, mk in (("ctx_unmasked", False), ("ctx_masked", True)):
            num = den = 0.0; num0 = 0.0
            for (t, f, ex), (_, fi, _) in zip(data, dinj):
                a, s = ctx_est(t, fi, ex, phi, mk); b, _ = ctx_est(t, f, ex, phi, mk)
                w = 1 / s ** 2; num += w * a; num0 += w * b; den += w
            rec[key].append((num - num0) / den)
    res[col] = {k: dict(mean_recovered_ppm=float(np.mean(v) * 1e6), sd_ppm=float(np.std(v) * 1e6), n=len(v),
                        bias_pct=float((np.mean(v) / 0.015 - 1) * 100)) for k, v in rec.items()}
    print(col, res[col], flush=True)
json.dump(res, open("out/r04b_depth_injection.json", "w"), indent=1)
