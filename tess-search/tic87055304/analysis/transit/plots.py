"""Figures: per-sector folds with the free-density MCMC model (fold_sectors.png), posterior corner plots
(corner_free.png and an overlay of the three density priors corner_priors.png), and residual diagnostics."""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import tmodel as tm

HERE = os.path.dirname(os.path.abspath(__file__))
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
ch = np.load(os.path.join(HERE, "mcmc_free.npz"))
chain = ch["chain"]; names = list(ch["names"])
med = np.array([mc["params"][n]["med"] for n in names])
data = tm.load()

# ---- per-sector folds
fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharey=True)
axes = axes.ravel()
T0, P = med[0], med[1]
rng = np.random.default_rng(1)
draws = chain[rng.choice(len(chain), 100, replace=False)]
binw = {13: 30, 39: 20, 66: 20, 93: 20, 104: 20}
for i, s in enumerate(tm.SECTORS + ["all"]):
    ax = axes[i]
    secs = tm.SECTORS if s == "all" else [s]
    for ss in secs:
        d = data[ss]
        j = tm.SECTORS.index(ss)
        c = med[7 + j]
        tau = (d["t"] - T0 - np.round((d["t"] - T0) / P) * P) * 24
        f = d["f"] / (1 + c)
        if s != "all":
            ax.plot(tau, (f - 1) * 1e6, ".", color="0.75", ms=1.5, zorder=0)
    # binned
    tt = np.concatenate([data[ss]["t"] for ss in secs])
    ff = np.concatenate([data[ss]["f"] / (1 + med[7 + tm.SECTORS.index(ss)]) for ss in secs])
    ww = np.concatenate([np.full(len(data[ss]["t"]), np.exp(-2 * med[12 + tm.SECTORS.index(ss)])) for ss in secs])
    tau = (tt - T0 - np.round((tt - T0) / P) * P) * 24
    bw = (binw[s] if s != "all" else 15) / 60
    be = np.arange(-0.45 * 24, 0.45 * 24 + 1e-6, bw)
    idx = np.digitize(tau, be)
    bc = []; bv = []; bs = []
    for k in range(1, len(be)):
        m = idx == k
        if m.sum() < 3: continue
        bc.append(0.5 * (be[k - 1] + be[k])); bv.append(np.sum(ff[m] * ww[m]) / np.sum(ww[m])); bs.append(1 / np.sqrt(np.sum(ww[m])))
    ax.errorbar(bc, (np.array(bv) - 1) * 1e6, np.array(bs) * 1e6, fmt="o", color="k", ms=3, lw=0.8, zorder=3)
    # model (exposure of that sector; 'all' shows 200 s)
    E = tm.EXPT[s] if s != "all" else 200.0
    tg = np.linspace(-0.45, 0.45, 1500)
    dm = {"x": dict(t=tg + T0, exptime=E)}
    for dr in draws[:60]:
        ax.plot(tg * 24, (tm.phys_model(dr[:7], dm)["x"] - 1) * 1e6, color="C0", alpha=0.05, lw=0.8, zorder=1)
    ax.plot(tg * 24, (tm.phys_model(med[:7], dm)["x"] - 1) * 1e6, color="C3", lw=1.3, zorder=2)
    nt = len(np.unique(np.round((tt - T0) / P)))
    lab = f"S{s} ({tm.EXPT[s]:.0f} s, {nt} windows, {binw[s]}-min bins)" if s != "all" else "all sectors (15-min bins; model at 200 s)"
    ax.set_title(lab, fontsize=9)
    ax.set_xlim(-10.8, 10.8); ax.set_ylim(-7000, 3500)
    ax.axhline(0, color="k", lw=0.4)
    if i >= 3: ax.set_xlabel("hours from mid-transit")
    if i % 3 == 0: ax.set_ylabel("relative flux (ppm)")
p = mc["params"]
fig.suptitle(f"TIC 87055304 SAP (biweight-detrended), free-density emcee fit: k = {p['k']['med']:.4f}, b = {p['b']['med']:.3f}, "
             f"rho* = {p['rho']['med']:.3f} rho_sun, depth = {p['depth_ppm']['med']:.0f} ppm, T14 = {p['T14_h']['med']:.2f} h", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fold_sectors.png"), dpi=130); plt.close(fig)

# ---- corner (free)
der = tm.derived(chain)
keys = ["k", "b", "rho", "u1", "u2", "P", "T0"]
X = np.vstack([chain[:, names.index(k)] for k in keys] + [der["T14_h"], der["inc"]]).T
X[:, 5] = (X[:, 5] - np.median(X[:, 5])) * 86400; X[:, 6] = (X[:, 6] - np.median(X[:, 6])) * 1440
labels = ["k=Rp/R*", "b", "rho* [rho_sun]", "u1", "u2", "P - med [s]", "T0 - med [min]", "T14 [h]", "inc [deg]"]
fig = corner.corner(X, labels=labels, quantiles=[0.1587, 0.5, 0.8413], show_titles=True, title_fmt=".3f", title_kwargs=dict(fontsize=8), label_kwargs=dict(fontsize=8))
fig.savefig(os.path.join(HERE, "corner_free.png"), dpi=110); plt.close(fig)

# ---- overlay of density priors (k, b, rho, T14)
fig = None
cols = {"free": "k", "dwarf": "C0", "subgiant": "C3", "dil051": "C2"}
for nm, c in cols.items():
    path = os.path.join(HERE, f"mcmc_{nm}.npz")
    if not os.path.exists(path):
        continue
    cc = np.load(path)["chain"]
    dd = tm.derived(cc)
    Y = np.vstack([cc[:, 2], cc[:, 3], np.log10(cc[:, 4]), dd["T14_h"]]).T
    fig = corner.corner(Y, fig=fig, color=c, labels=["k", "b", "log10 rho*", "T14 [h]"], plot_datapoints=False, plot_density=False,
                        levels=(0.393, 0.865), smooth=1.0, hist_kwargs=dict(density=True), range=[(0.04, 0.12), (0, 1.1), (-1.8, 0.0), (3.8, 5.2)])
if fig is not None:
    fig.legend(handles=[matplotlib.lines.Line2D([], [], color=c, label=n) for n, c in cols.items() if os.path.exists(os.path.join(HERE, f"mcmc_{n}.npz"))],
               loc="upper right", fontsize=9)
    fig.savefig(os.path.join(HERE, "corner_priors.png"), dpi=110)
print("ok")
