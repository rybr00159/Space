import time, sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
s = int(sys.argv[1])
t0 = time.time()
S = setup(s, "none", jitter_from="phot"); print("setup", time.time() - t0)
R = photometry(s, S); print("phot", time.time() - t0)
z = np.load(f"{PH}/data/phot_s{s}.npz")
a = R["FIT"]["B"]["F"][:, 1]; b = z["F_B"][:, 1]
k = np.isfinite(a) & np.isfinite(b)
print("max |diff| M3V_B vs phot:", np.max(np.abs(a[k] - b[k])), "rms", np.std(b[k]))
t1 = time.time(); J = jitter_clean(s, "none", S["C"], save=False); print("jitter", time.time() - t1)
Jp = np.load(f"{PH}/data/jitter_s{s}.npz")
for kk in ("dx", "dy", "bkg"):
    k = np.isfinite(J[kk]) & np.isfinite(Jp[kk]); print(kk, np.max(np.abs(J[kk][k] - Jp[kk][k])))
