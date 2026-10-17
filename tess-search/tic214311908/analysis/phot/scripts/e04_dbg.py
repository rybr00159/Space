import sys; sys.path.insert(0, ".")
from photlib import *
from wotan import flatten
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
s = int(sys.argv[1])
z = np.load(f"{ROOT}/data/phot_s{s}.npz"); names = list(z["names"])
t = z["t"]; ok = z["qok"]
fig, ax = plt.subplots(5, 1, figsize=(14, 14), sharex=True)
for k, nm in enumerate(("M3V", "pair", "Mira")):
    f = z["F"][:, names.index(nm)]
    ax[k].plot(t[ok], f[ok], ",", alpha=0.5); ax[k].set_ylabel(nm)
    m = ok & np.isfinite(f)
    tr = flatten(t[m], f[m], method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.3, edge_cutoff=0)[1]
    ax[k].plot(t[m], tr, "r-")
    lo, hi = np.nanpercentile(f[ok], [1, 99]); ax[k].set_ylim(lo - 50, hi + 50)
ax[3].plot(t[ok], z["bkg"][ok], ","); ax[3].set_ylabel("bkg")
ax[4].plot(t[ok], z["dx"][ok], ",", label="dx"); ax[4].plot(t[ok], z["dy"][ok], ",", label="dy"); ax[4].legend()
for tc in T0_EPH + P_EPH * np.arange(-300, 100):
    if t[ok].min() < tc < t[ok].max():
        for a in ax: a.axvline(tc, color="g", lw=0.5)
fig.tight_layout(); fig.savefig(f"{ROOT}/explore/e04_dbg_s{s}.png", dpi=60)
