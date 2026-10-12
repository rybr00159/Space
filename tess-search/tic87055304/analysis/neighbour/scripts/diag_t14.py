"""Diagnostic: maximum-likelihood transit shape (k, b, rho, T14) per light curve (pair PSF vs QLP SAP) and per sector
subset, to understand the T14 difference between the pair MCMC (4.22 h) and the first-round QLP SAP MCMC (4.81 h).
D fixed: pair 5.32, QLP SAP 4.39 (empirical S13-S93, results/fluxratio.json) / 6.0 for S104. Noise per sector = OOT MAD.
Output: results/diag_t14.json"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbmodel as tm
ROOT = tm.NBROOT
out = {}
for lcname, Ddef in (("pair", 5.32), ("qlp", 4.39)):
    for subset in ([13, 39, 66, 93, 104], [13, 39, 66, 93], [13], [39], [66], [93], [104]):
        data = tm.load(os.path.join(ROOT, f"data/lc_{lcname}.npy"), sectors=subset)
        D = Ddef if not (lcname == "qlp" and subset == [104]) else 6.0
        def chi2(x):
            T0, k, b, rho = x
            if not (0.05 < k < 0.3 and 0 <= b < 1 + k and 0.02 < rho < 3): return 1e12
            mm = tm.phys_model([T0, 4.1362, k, b, rho, 0.306, 0.272], data, D)
            c = 0
            for s in subset:
                d = data[s]; r = d["f"] - mm[s]; r -= np.median(r[np.abs(d["t"] - T0 - np.round((d["t"] - T0) / 4.1362) * 4.1362) > 0.15])
                c += np.sum(r ** 2) / d["sigma0"] ** 2
            return c
        best = None
        for b0 in (0.1, 0.4, 0.6, 0.8):
            for r0 in (0.15, 0.3, 0.5):
                x = minimize(chi2, [3271.5326, 0.14, b0, r0], method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-6, fatol=1e-3)).x
                v = chi2(x)
                if best is None or v < best[0]: best = (v, x)
        T0, k, b, rho = best[1]
        a = tm.a_from_rho(rho, 4.1362)
        der = tm.derived(np.array([[T0, 4.1362, k, b, rho, 0.3, 0.27]]))
        key = f"{lcname}_{'-'.join(map(str, subset))}"
        out[key] = dict(chi2=float(best[0]), k=float(k), b=float(b), rho=float(rho), a=float(a), T14_h=float(der["T14_h"][0]), T23_h=float(der["T23_h"][0]))
        print(key, {kk: round(v, 4) for kk, v in out[key].items()}, flush=True)
json.dump(out, open(os.path.join(ROOT, "results/diag_t14.json"), "w"), indent=1)
