"""Summary figures: per-sector light curves of the three stars, individual M3V transits, per-sector folds."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from tlib import trap_frac
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL, DEP = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"], TR["depth"]


def binlc(t, f, w):
    if len(t) == 0: return np.array([]), np.array([])
    k = np.floor((t - t.min()) / w).astype(int); u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    return (np.bincount(inv, t) / n)[n >= 2], (np.bincount(inv, f) / n)[n >= 2]


# 1. per-sector light curves
fig, axs = plt.subplots(len(SECTORS), 1, figsize=(15, 3.0 * len(SECTORS)))
for ax, s in zip(axs, SECTORS):
    L = np.load(f"{ROOT}/data/lc_s{s}.npz"); t = L["t_tdb"]; ok = L["quality_ok"]
    for k, (key, lab, off) in enumerate((("M3V_B", "M3V (model B)", 0), ("pair_A", "pair (model A)", 60), ("Mira_B", "Mira (model B, x10)", 110))):
        d = (L[f"det_{key}"] - 1) * 1e3 * (10 if key.startswith("Mira") else 1)
        m = ok & np.isfinite(d)
        tb, fb = binlc(t[m], d[m], 1 / 24)
        ax.plot(tb, fb + off, ".", ms=2.5, label=lab)
    for e in np.arange(np.ceil((np.nanmin(t) - T0) / P), np.floor((np.nanmax(t) - T0) / P) + 1):
        ax.axvline(T0 + e * P, color="r", lw=0.6, alpha=0.6)
    ax.set_ylim(-60, 150); ax.set_ylabel("ppt (offset)"); ax.set_title(f"Sector {s}: deblended, detrended, 1-h bins; red = predicted transits")
axs[0].legend(loc="upper right", fontsize=8); axs[-1].set_xlabel("BTJD (TDB, target-barycentred)")
fig.tight_layout(); fig.savefig(f"{ROOT}/figs/lc_sectors.png", dpi=70); plt.close(fig)

# 2. individual transits (strict, model B)
PT = pd.read_csv(f"{ROOT}/out/per_transit_M3V.csv")
st = PT[(PT.model == "B") & (~PT.relaxed)].sort_values("epoch").reset_index(drop=True)
rl = PT[(PT.model == "B") & (PT.relaxed) & (~PT.epoch.isin(st.epoch))].sort_values("epoch").reset_index(drop=True)
M = pd.read_csv(f"{ROOT}/lc_M3V_all.csv")
allev = pd.concat([st.assign(kind="strict"), rl.assign(kind="relaxed")], ignore_index=True)
n = len(allev); nc = 6; nr = int(np.ceil(n / nc))
fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc, 2.6 * nr), sharey=True)
for ax, r in zip(axs.ravel(), allev.itertuples()):
    tc0 = T0 + r.epoch * P
    if r.kind == "strict":
        d = M[(np.abs(M.time_btjd_tdb - tc0) < 0.3) & (M.quality_ok == 1)]; f = d.flux_norm.values
    else:
        bmed = M.groupby("sector").bkg_e_s_px.transform("median")
        d = M[(np.abs(M.time_btjd_tdb - tc0) < 0.3) & ((M.bkg_e_s_px - bmed) < 800)]; f = d.flux_raw_norm.values
        f = f / np.nanmedian(f[np.abs(d.time_btjd_tdb.values - tc0) > T14])
    x = (d.time_btjd_tdb.values - tc0) * 24
    ax.plot(x, (f - 1) * 1e3, ".", color="0.7", ms=2)
    xb, fb = binlc(x / 24, (f - 1) * 1e3, 20 / 1440); ax.plot(xb * 24, fb, "o", ms=3, color="C0" if r.kind == "strict" else "C1")
    xm = np.linspace(-7, 7, 400) / 24
    ax.plot(xm * 24 , -r.depth * 1e3 * trap_frac(xm + tc0, r.tc, T14, FL, 0.0014), "r-", lw=1)
    ax.set_title(f"S{r.sector} E{r.epoch} {r.kind}: {r.depth*1e3:.0f}+-{max(r.depth_err, r.depth_err_pb)*1e3:.0f} ppt, O-C {r.oc_min:+.0f} min", fontsize=7)
    ax.set_ylim(-60, 40); ax.axvline(0, color="k", lw=0.4)
for ax in axs.ravel()[n:]: ax.axis("off")
fig.tight_layout(); fig.savefig(f"{ROOT}/figs/transits_M3V.png", dpi=70); plt.close(fig)

# 3. per-sector folds for every light curve
fig, axs = plt.subplots(2, len(SECTORS), figsize=(3.2 * len(SECTORS), 6.5), sharey="row")
for j, s in enumerate(SECTORS):
    L = np.load(f"{ROOT}/data/lc_s{s}.npz"); t = L["t_tdb"]; ok = L["quality_ok"]
    ph = (((t - T0) / P + 0.5) % 1 - 0.5) * P * 24
    stt = json.load(open(f"{ROOT}/out/lc_stats.json"))[str(s)]["fref"]
    for i, keys in enumerate((("M3V_B", "M3V_A"), ("pair_A", "pair_Bp", "Mira_B"))):
        for key in keys:
            d = (L[f"det_{key}"] - 1) * stt[key.split("_")[0]]     # e/s
            m = ok & np.isfinite(d) & (np.abs(ph) < 6)
            b = np.arange(-6, 6.01, 0.4); idx = np.digitize(ph[m], b)
            yb = np.array([np.nanmean(d[m][idx == k]) if np.sum(idx == k) > (0 if s == 12 else 3) else np.nan for k in range(1, len(b))])
            axs[i, j].plot(b[:-1] + 0.2, yb, ".-", ms=3, label=key)
        axs[i, j].axhline(0, color="k", lw=0.4); axs[i, j].set_title(f"S{s}", fontsize=9)
        xm = np.linspace(-6, 6, 300) / 24
        axs[i, j].plot(xm * 24, -DEP * 940 * trap_frac(xm, 0.0, T14, FL, 0.0014), "k:", lw=0.8)
    axs[0, j].set_ylim(-35, 20); axs[1, j].set_ylim(-60, 60)
axs[0, 0].set_ylabel("e/s (M3V)"); axs[1, 0].set_ylabel("e/s (pair, Mira)"); axs[0, 0].legend(fontsize=7); axs[1, 0].legend(fontsize=7)
for a in axs[1]: a.set_xlabel("h from mid-transit")
fig.tight_layout(); fig.savefig(f"{ROOT}/figs/fold_by_sector.png", dpi=75); plt.close(fig)
print("figures written")
