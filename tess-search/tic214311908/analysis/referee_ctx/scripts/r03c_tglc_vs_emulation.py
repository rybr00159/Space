"""R03c: is the real TGLC S39 dF vector consistent with a linear emulation of TGLC applied to the same pixels (r03)?
Difference vector (TGLC - emulation), both from the same transit windows and the same local-baseline estimator; its
size is compared with ctx's covariance B (r02) and projected on B's eigenvectors.
Output: out/r03c_tglc_vs_emulation.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import r02_leak_gls as A
r3 = json.load(open("out/r03_pixels_s39.json"))
names = A.names
emu = np.array([r3["emulated_tglc_psf_dF"][k] for k in names])
tg = np.array([r3["real_tglc_psf_dF_referee"][k] for k in names])
d = tg - emu
w, V = np.linalg.eigh(A.CB)
proj = V.T @ d
out = dict(tglc=tg.tolist(), emulated=emu.round(2).tolist(), diff=d.round(2).tolist(),
           chi2_diff_covB=float(d @ np.linalg.inv(A.CB) @ d), per_eigen_sigma=np.sqrt(w).round(3).tolist(),
           per_eigen_proj=proj.round(3).tolist(), per_eigen_chi2=(proj ** 2 / w).round(2).tolist(),
           pairA_minus_pairB=dict(tglc=float(tg[1] - tg[2]), emulated=float(emu[1] - emu[2])),
           M3V_minus_pairA=dict(tglc=float(tg[0] - tg[1]), emulated=float(emu[0] - emu[1])))
print(json.dumps(out, indent=1))
json.dump(out, open("out/r03c_tglc_vs_emulation.json", "w"), indent=1)
