"""Diagnostic 2: same as diag_t14.py but on the UNDETRENDED normalised flux (f_raw) with a per-transit quadratic
baseline solved linearly (variable projection) inside +-0.45 d windows, i.e. no global biweight. If the pair/QLP T14
difference survives, it is in the photometry, not the detrending. Also profile chi2 vs T14 (b fixed on a grid).
Output: results/diag_t14b.json"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbmodel as tm
ROOT = tm.NBROOT
P = 4.1362006
out = {}


def prep(lcname, subset):
    lc = np.load(os.path.join(ROOT, f"data/lc_{lcname}.npy"), allow_pickle=True).item()
    blocks = []
    for s in subset:
        v = lc[s]; t = v["t"]; y = v["f_raw"]
        n = np.round((t - 3271.5326) / P); tau = t - 3271.5326 - n * P
        for e in np.unique(n[np.abs(tau) < 0.45]):
            m = (n == e) & (np.abs(tau) < 0.45)
            if np.sum(m & (np.abs(tau) < 0.08)) < 3 or np.sum(m & (np.abs(tau) > 0.15)) < 6: continue
            blocks.append(dict(s=s, t=t[m], y=y[m] / np.median(y[m]), sig=v["sigma"], E=tm.EXPT[s]))
    return blocks


def chi2_blocks(x, blocks, D):
    T0, k, b, rho = x
    if not (0.05 < k < 0.3 and 0 <= b < 1 + k and 0.02 < rho < 3): return 1e12
    data = {i: dict(t=bl["t"], exptime=bl["E"]) for i, bl in enumerate(blocks)}
    mm = tm.phys_model([T0, P, k, b, rho, 0.306, 0.272], data, D)
    c = 0.0
    for i, bl in enumerate(blocks):
        x_ = (bl["t"] - bl["t"].mean())
        A = np.vstack([mm[i], mm[i] * x_, mm[i] * x_ ** 2]).T
        coef, *_ = np.linalg.lstsq(A, bl["y"], rcond=None)
        c += np.sum((bl["y"] - A @ coef) ** 2) / bl["sig"] ** 2
    return c


for lcname, D in (("pair", 5.32), ("qlp", 4.39)):
    for subset in ([13, 39, 66, 93, 104], [13, 39, 66, 93], [66, 93, 104], [39], [66], [93], [104]):
        Dd = 6.0 if (lcname == "qlp" and subset == [104]) else D
        bl = prep(lcname, subset)
        best = None
        for b0 in (0.1, 0.5, 0.75):
            for r0 in (0.15, 0.4):
                r = minimize(chi2_blocks, [3271.5326, 0.14, b0, r0], args=(bl, Dd), method="Nelder-Mead",
                             options=dict(maxiter=3000, xatol=1e-6, fatol=1e-3))
                if best is None or r.fun < best[0]: best = (r.fun, r.x)
        T0, k, b, rho = best[1]
        der = tm.derived(np.array([[T0, P, k, b, rho, 0.3, 0.27]]))
        key = f"{lcname}_{'-'.join(map(str, subset))}"
        out[key] = dict(n_tr=len(bl), chi2=float(best[0]), k=float(k), b=float(b), rho=float(rho), T14_h=float(der["T14_h"][0]),
                        T23_h=float(der["T23_h"][0]))
        print(key, {kk: round(v, 4) for kk, v in out[key].items()}, flush=True)
json.dump(out, open(os.path.join(ROOT, "results/diag_t14b.json"), "w"), indent=1)
