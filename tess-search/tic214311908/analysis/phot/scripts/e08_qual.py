import sys; sys.path.insert(0, ".")
from photlib import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(6, 1, figsize=(16, 22))
for ax, s in zip(axs, SECTORS):
    z = np.load(f"{ROOT}/data/phot_s{s}.npz"); J = json.load(open(f"{ROOT}/out/phot_s{s}.json")); L = np.load(f"{ROOT}/data/lc_s{s}.npz")
    t = z["t"]; f = z["F_B"][:, 1]; qok = z["qok"]
    chi = z["chi2_B"]; ch = chi / np.nanmedian(chi[z["q0"]])
    # simple running-median detrend over 0.75 d for all qok cadences
    from scipy.ndimage import median_filter
    o = np.argsort(t); ff = f.copy()
    ax.plot(t[qok], f[qok], ",", color="k", alpha=0.5)
    bad = qok & ~L["quality_ok"]
    ax.plot(t[bad], f[bad], ".", ms=2, color="r", label="rejected (current flags)")
    ax2 = ax.twinx(); ax2.plot(t[qok], z["bkg"][qok], "g-", lw=0.6); ax2.axhline(J["bkg_lim"], color="g", ls=":"); ax2.set_ylabel("bkg")
    ax3 = ax.twinx(); ax3.plot(t[qok], ch[qok], "m-", lw=0.4, alpha=0.6); ax3.set_ylim(0, 20); ax3.axhline(4, color="m", ls=":")
    ax.set_ylim(np.nanpercentile(f[qok], 0.5) - 30, np.nanpercentile(f[qok], 99.5) + 30)
    for tc in T0_EPH + P_EPH * np.arange(-300, 100):
        if t[qok].min() < tc < t[qok].max(): ax.axvline(tc, color="b", lw=0.6)
    ax.set_title(f"S{s} M3V model B dF (e/s), bkg (green), chi2r/median (magenta)")
fig.tight_layout(); fig.savefig(f"{ROOT}/explore/e08_qual.png", dpi=55)
