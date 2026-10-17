"""R02c: which directions of the 6-star dF space carry ctx's chi2? Eigen-decomposition of the ctx-style covariance B,
residuals of the M3V and best free-position models projected on each eigenvector.
Output: out/r02c_eigen.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import r02_leak_gls as A
names = A.names
w, V = np.linalg.eigh(A.CB)
Ci = np.linalg.inv(A.CB)
out = {"eigenvalues": w.tolist(), "eigenvectors_columns": np.round(V, 3).tolist(), "names": names}
for h, m in (("M3V", A.HYPV["M3V"]), ("G17.7_5.5", A.HYPV["G17.7_5.5"])):
    d, e, c2 = A.gls(A.o_ctx, Ci, m); r = A.o_ctx - d * m
    proj = V.T @ r; contrib = proj ** 2 / w
    out[h] = dict(dF=d, chi2=c2, resid=r.round(3).tolist(), chi2_per_eigen=contrib.round(2).tolist(),
                  sigma_per_eigen=np.sqrt(w).round(3).tolist(), resid_proj=proj.round(3).tolist())
    print(h, "chi2", round(c2, 2), "per-eigen", contrib.round(2), "sigma", np.sqrt(w).round(3))
for i in range(len(w)):
    print(f"eig {i}: sigma {np.sqrt(w[i]):6.3f}  vec", dict(zip(names, V[:, i].round(2))))
# fractional coefficient precision needed: sigma_min / (dF * |V_min . c|) etc.
json.dump(out, open("out/r02c_eigen.json", "w"), indent=1, default=float)
