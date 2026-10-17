"""Figure: TIC 214311908 transit folded in independent light-curve products (detrended, normalised; from c04's
out/fold_arrays.npz). Grey points = individual cadences (shown only for 200-s/600-s data); blue = 15-min bins (30-min for Sector 12 panels).
Output: fig/fold_products.png"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import common as C
a = np.load("out/fold_arrays.npz")
panels = [("QLP M3V (DET_FLUX), S39-S104", [f"QLP_214311908_{s}_DET_FLUX" for s in (39, 66, 93, 103, 104)]),
          ("QLP Mira TIC 214311906 (DET_FLUX), S39-S104", [f"QLP_214311906_{s}_DET_FLUX" for s in (39, 66, 93, 103, 104)]),
          ("QLP pair TIC 214311910 (DET_FLUX), S39-S104", [f"QLP_214311910_{s}_DET_FLUX" for s in (39, 66, 93, 103, 104)]),
          ("TGLC M3V psf_flux, S39 (Gaia-deblended)", ["TGLC_5925373272850784768_39_psf_flux"]),
          ("TGLC Mira psf_flux, S39", ["TGLC_5925373268528684288_39_psf_flux"]),
          ("TARS M3V (1-pixel aperture), S12-S93", [f"TARS_214311908_{s}_FLUX" for s in (12, 39, 66, 93)]),
          ("GSFC-ELEANOR-LITE S12 (PCA_FLUX; blend)", ["ELEANOR-LITE_214311908_12_PCA_FLUX"]),
          ("TGLC M3V psf_flux, S12 (30-min)", ["TGLC_5925373272850784768_12_psf_flux"])]
fig, axs = plt.subplots(4, 2, figsize=(10, 11), sharex=True)
for ax, (title, keys) in zip(axs.ravel(), panels):
    t = np.concatenate([a[k][0] for k in keys]); f = np.concatenate([a[k][1] for k in keys])
    ph = ((t - C.T0) / C.P + 0.5) % 1 - 0.5; dt = ph * C.P * 24
    m = np.abs(dt) < 4
    ax.plot(dt[m], (f[m] - 1) * 100, ".", ms=1.5, color="#b8b7b0", zorder=1)
    bw = 0.5 if all("_12_" in k for k in keys) else 0.25; bins = np.arange(-4, 4.01, bw); idx = np.digitize(dt[m], bins)
    bc = 0.5 * (bins[1:] + bins[:-1]); bm = np.array([np.nanmean(f[m][idx == i]) if np.sum(idx == i) >= 2 else np.nan for i in range(1, len(bins))])
    be = np.array([np.nanstd(f[m][idx == i]) / np.sqrt(np.sum(idx == i)) if np.sum(idx == i) >= 2 else np.nan for i in range(1, len(bins))])
    ax.errorbar(bc, (bm - 1) * 100, (be) * 100, fmt="o", ms=4, color="#2a78d6", lw=1, zorder=2)
    lo = np.nanpercentile((bm - 1) * 100, 1); hi = np.nanpercentile((bm - 1) * 100, 99); pad = max(0.3, 0.6 * (hi - lo))
    ax.set_ylim(lo - pad, hi + pad); ax.axvspan(-C.T14 * 12, C.T14 * 12, color="#eda100", alpha=0.12, lw=0)
    ax.set_title(title, fontsize=9); ax.grid(alpha=0.25, lw=0.5)
for ax in axs[-1]: ax.set_xlabel("hours from predicted mid-transit")
for ax in axs[:, 0]: ax.set_ylabel("relative flux - 1 (%)")
fig.suptitle("TIC 214311908 candidate (P = 5.6914968 d) in independent TESS light-curve products", fontsize=11)
plt.tight_layout(); plt.savefig("fig/fold_products.png", dpi=90)
print("saved fig/fold_products.png")
