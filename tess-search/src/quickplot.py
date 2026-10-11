"""Quick look: per-sector phase folds (full phase and zoom) plus a Lomb-Scargle periodogram."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

warnings.filterwarnings("ignore")
import stage2


def binned(x, y, edges):
    idx = np.digitize(x, edges); c, v = [], []
    for k in range(1, len(edges)):
        m = idx == k
        if m.sum() >= 3: c.append(0.5 * (edges[k - 1] + edges[k])); v.append(np.mean(y[m]))
    return np.array(c), np.array(v)


def quick(tic, P, t0, dur, out, secs=None, title=""):
    secs = secs or stage2.load_all(tic)
    n = len(secs)
    fig, ax = plt.subplots(3, n, figsize=(3.0 * n, 7.5), squeeze=False)
    for i, s in enumerate(secs):
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5
        c, v = binned(ph, (s["f"] - 1) * 1e6, np.linspace(-0.5, 0.5, 101))
        ax[0, i].plot(ph, (s["f"] - 1) * 1e6, ",", color="0.75"); ax[0, i].plot(c, v, "o-", ms=2, lw=0.6)
        lim = max(4 * np.std(v) if len(v) else 1000, 300); ax[0, i].set_ylim(-lim, lim)
        ax[0, i].set_title(f"S{s['sector']} {s['prov']}", fontsize=8)
        h = ph * P * 24; m = np.abs(h) < 5 * dur * 24
        c, v = binned(h[m], (s["f"][m] - 1) * 1e6, np.linspace(-5 * dur * 24, 5 * dur * 24, 41))
        ax[1, i].plot(h[m], (s["f"][m] - 1) * 1e6, ".", ms=1, color="0.75"); ax[1, i].plot(c, v, "o-", ms=3, lw=0.8)
        ax[1, i].set_ylim(-lim, lim); ax[1, i].set_xlabel("hours", fontsize=7)
        fr, pw = LombScargle(s["t"], s["f"]).autopower(minimum_frequency=0.05, maximum_frequency=48, samples_per_peak=5)
        ax[2, i].plot(fr, pw, lw=0.5); ax[2, i].axvline(1 / P, color="r", lw=0.5, ls=":")
        ax[2, i].set_xlabel("cycles/day", fontsize=7)
        top = fr[np.argmax(pw)]; ax[2, i].set_title(f"LS peak {top:.3f} c/d", fontsize=7)
    fig.suptitle(title or f"TIC {tic} P={P:.5f} d", fontsize=10); fig.tight_layout(); fig.savefig(out, dpi=100); plt.close(fig)
    return secs


if __name__ == "__main__":
    tic, P, t0, dur, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    quick(tic, P, t0, dur, out)
