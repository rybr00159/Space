"""Data preparation for the transit-model analysis of TIC 87055304.

Builds data/transit_lc.npz-like dict (saved in this folder as lc_prepared.npy) from the QLP SAP_FLUX:
  * QUALITY == 0 only, finite flux;
  * S104 scattered-light event BTJD 4178.80-4179.20 removed;
  * momentum dumps (QUALITY bit 32, value 32): detrending segments are split at each dump and cadences within
    DUMP_PAD of a dump are removed (documented in README);
  * per-sector detrending: wotan biweight, window WIN d, transit-masked (mask full width MASKW d around the
    ephemeris), segments split at gaps > 0.5 d and at dumps;
  * also stores the raw normalised SAP flux for local-baseline fits (timing, systematics checks).
Also writes diagnostic plots of raw flux, background and centroids with transit and dump times marked.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wotan import flatten, transit_mask

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = "/home/user/work/tic87/data"
P0, T00, DUR0 = 4.136195, 1658.41325, 0.1875
EXPT = {13: 1800.0, 39: 600.0, 66: 200.0, 93: 200.0, 104: 200.0}
SCAT = (4178.80, 4179.20)
DUMP_PAD = 1.0 / 24  # remove +-1 h around momentum dumps (see README)
GAP = 0.5


def segments(t, breaks):
    """Integer segment labels: new segment after gaps > GAP or at any break time."""
    lab = np.zeros(len(t), int)
    k = 0
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > GAP or np.any((breaks > t[i - 1]) & (breaks <= t[i])):
            k += 1
        lab[i] = k
    return lab


def build(win=1.0, maskw=0.30, dump_pad=DUMP_PAD, method="biweight", P=P0, T0=T00):
    raw = np.load(f"{DATA}/qlp_raw_columns.npy", allow_pickle=True).item()
    out = {}
    for s in sorted(raw):
        v = raw[s]
        t = np.asarray(v["TIME"], float); f = np.asarray(v["SAP_FLUX"], float); q = np.asarray(v["QUALITY"], np.int64)
        bkg = np.asarray(v["SAP_BKG"], float); x = np.asarray(v["SAP_X"], float); y = np.asarray(v["SAP_Y"], float)
        dumps = np.unique(np.round(t[(q & 32) > 0], 3))
        ok = (q == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0)
        ok &= ~((t > SCAT[0]) & (t < SCAT[1]))
        if dump_pad > 0 and len(dumps):
            near = np.min(np.abs(t[:, None] - dumps[None, :]), axis=1) < dump_pad
            ok &= ~near
        t, f, bkg, x, y = t[ok], f[ok], bkg[ok], x[ok], y[ok]
        o = np.argsort(t); t, f, bkg, x, y = t[o], f[o], bkg[o], x[o], y[o]
        f = f / np.nanmedian(f)
        lab = segments(t, dumps)
        mask = transit_mask(time=t, period=P, duration=maskw, T0=T0)
        trend = np.full_like(f, np.nan)
        for k in np.unique(lab):
            m = lab == k
            if m.sum() < 5:
                continue
            ff = f[m].copy(); ff[mask[m]] = np.nan
            if np.isfinite(ff).sum() < 5:
                continue
            tr = flatten(t[m], ff, method=method, window_length=win, return_trend=True, edge_cutoff=0)[1]
            # interpolate trend across masked transit cadences
            g = np.isfinite(tr)
            if g.sum() >= 2:
                tr = np.interp(t[m], t[m][g], tr[g])
            trend[m] = tr
        good = np.isfinite(trend)
        fd = f / trend
        # robust per-sector white-noise estimate from out-of-transit flattened flux
        oot = good & ~mask
        mad = 1.4826 * np.nanmedian(np.abs(fd[oot] - np.nanmedian(fd[oot])))
        # 5-sigma clip of out-of-transit outliers only (in-transit points never clipped)
        clip = oot & (np.abs(fd - 1) > 5 * mad)
        keep = good & ~clip
        out[s] = dict(t=t[keep], f=fd[keep], f_raw=f[keep], trend=trend[keep], bkg=bkg[keep], x=x[keep], y=y[keep],
                      seg=lab[keep], exptime=EXPT[s], sigma=mad, dumps=dumps, nclip=int(clip.sum()))
    return out


def diag_plot(lc, png):
    fig, axes = plt.subplots(len(lc), 1, figsize=(14, 2.6 * len(lc)))
    for ax, s in zip(axes, sorted(lc)):
        v = lc[s]
        ax.plot(v["t"], v["f_raw"], ",", color="0.5")
        ax.plot(v["t"], v["trend"], "-", color="C1", lw=0.8)
        n = np.arange(np.ceil((v["t"].min() - T00) / P0), np.floor((v["t"].max() - T00) / P0) + 1)
        for tt in T00 + n * P0:
            ax.axvspan(tt - DUR0 / 2, tt + DUR0 / 2, color="C0", alpha=0.25)
        for d in v["dumps"]:
            ax.axvline(d, color="C3", lw=0.6, ls=":")
        ax2 = ax.twinx(); ax2.plot(v["t"], v["bkg"], ",", color="C2"); ax2.set_ylabel("bkg", color="C2", fontsize=7)
        ax.set_ylabel(f"S{s} SAP"); lo, hi = np.percentile(v["f_raw"], [0.5, 99.5]); ax.set_ylim(lo - 0.004, hi + 0.003)
    axes[-1].set_xlabel("BTJD")
    fig.tight_layout(); fig.savefig(png, dpi=110); plt.close(fig)


if __name__ == "__main__":
    lc = build()
    np.save(os.path.join(HERE, "lc_prepared.npy"), lc, allow_pickle=True)
    for s in sorted(lc):
        v = lc[s]
        print(s, len(v["t"]), "sigma %.0f ppm" % (v["sigma"] * 1e6), "clipped", v["nclip"], "dumps", v["dumps"])
    diag_plot(lc, os.path.join(HERE, "diag_raw_sap.png"))
