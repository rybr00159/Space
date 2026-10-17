"""s20: aggregation - empirical false-alarm thresholds of the BLS searches, significance of the real peaks,
injection-recovery completeness map, summary figures.

Null families (s02, all on the RRab-cleaned light curve, identical search):
  M3V-based  : 'roll' (rolled within sectors, not inverted), 'inv_roll' (inverted + rolled), and the 13 inverted+rolled
               global maxima of an earlier interrupted run (logs/s02_inv_roll_1_killed.log, global maxima only);
  other LCs  : deblended pair (model A) and Mira (model B), both signs; control stars and phantom positions (s03; these
               contain leaked real variables, e.g. the 1.595-d Gaia EB, so they are conservative).
Threshold: Gumbel fit (method of moments) to the M3V-based global maxima -> Delta-chi2 at global FAP 1 % and 0.1 %;
empirical p-value of each real peak = fraction of null maxima >= it. Per-sector thresholds from the same realisations.
Completeness (s09): recovered if Delta-chi2 at the injected period/phase >= threshold (FAP 1 %).
Outputs: out/summary.json, figs/bls_clean.png, figs/null_dist.png, figs/completeness.png
"""
import sys, glob, re; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def gumbel_fit(x):
    x = np.asarray(x); beta = np.std(x) * np.sqrt(6) / np.pi; mu = np.mean(x) - 0.5772 * beta
    return mu, beta


def gumbel_q(mu, beta, fap):
    return float(mu - beta * np.log(-np.log(1 - fap)))


def gumbel_p(mu, beta, x):
    return float(1 - np.exp(-np.exp(-(x - mu) / beta)))


def load_nulls():
    fam = {}
    for f in glob.glob(f"{ROOT}/out/null_*_clean_*.json"):
        d = json.load(open(f)); fam.setdefault(d["kind"], []).extend(d["res"])
    killed = []
    fn = f"{ROOT}/logs/s02_inv_roll_1_killed.log"
    if os.path.exists(fn):
        for line in open(fn):
            m = re.match(r"^\d+ ([\d.]+) P ([\d.]+)", line)
            if m: killed.append(float(m.group(1)))
    return fam, killed


def main():
    real = json.load(open(f"{ROOT}/out/bls_real_clean.json"))
    fam, killed = load_nulls()
    m3v_glob = [r["max"] for k in ("roll", "inv_roll") for r in fam.get(k, [])] + killed
    out = dict(n_null_m3v=len(m3v_glob), null_families={k: len(v) for k, v in fam.items()}, n_killed_inv=len(killed))
    mu, beta = gumbel_fit(m3v_glob)
    thr1, thr01 = gumbel_q(mu, beta, 0.01), gumbel_q(mu, beta, 0.001)
    out["global"] = dict(null_max_median=float(np.median(m3v_glob)), null_max_p90=float(np.percentile(m3v_glob, 90)),
                         null_max_max=float(np.max(m3v_glob)), gumbel_mu=mu, gumbel_beta=beta, thr_fap1=thr1, thr_fap01=thr01)
    for k, v in fam.items():
        out["global"][f"null_{k}_maxima"] = [round(r["max"], 1) for r in v]
    # real peaks
    pk = []
    for p in real["all_clean"]["peaks"][:8]:
        pk.append(dict(P=p["P"], dchi2=p["dchi2"], depth_ppt=p["depth_ppt"], n_events=p["n_events"], max_event_frac=p["max_event_frac"],
                       p_emp=float(np.mean(np.array(m3v_glob) >= p["dchi2"])), p_gumbel=gumbel_p(mu, beta, p["dchi2"])))
    out["real_peaks"] = pk
    # per sector
    sec = {}
    for g in ("S12", "S39", "S66", "S93", "S103", "S104", "S103_104"):
        vals = [r["sectors"][g]["max"] for k in ("roll", "inv_roll") for r in fam.get(k, []) if g in r.get("sectors", {})]
        if len(vals) < 5: continue
        m_, b_ = gumbel_fit(vals)
        rp = real[g]["peaks"][0]
        sec[g] = dict(n=len(vals), null_median=float(np.median(vals)), thr_fap1=gumbel_q(m_, b_, 0.01), real_max=rp["dchi2"], real_P=rp["P"],
                      real_n_events=rp["n_events"], real_max_event_frac=rp["max_event_frac"],
                      p_emp=float(np.mean(np.array(vals) >= rp["dchi2"])), p_gumbel=gumbel_p(m_, b_, rp["dchi2"]))
    out["per_sector"] = sec
    # TLS
    tls = json.load(open(f"{ROOT}/out/tls.json"))
    tn = json.load(open(f"{ROOT}/out/tls_null_S103_104.json")) if os.path.exists(f"{ROOT}/out/tls_null_S103_104.json") else []
    out["tls"] = dict(per_sector={k: dict(SDE=v["SDE"], P=v["P"], n_transits=v["n_transits"]) for k, v in tls.items()},
                      S103_104_null_SDE_roll=[r["SDE"] for r in tn if r["sign"] == 1], S103_104_null_SDE_invroll=[r["SDE"] for r in tn if r["sign"] == -1])
    if tn:
        sd = np.array([r["SDE"] for r in tn])
        out["tls"]["S103_104_p_emp_all"] = float(np.mean(sd >= tls["S103_104"]["SDE"]))
        out["tls"]["S103_104_p_emp_roll"] = float(np.mean(np.array(out["tls"]["S103_104_null_SDE_roll"]) >= tls["S103_104"]["SDE"]))
    # completeness
    files = sorted(glob.glob(f"{ROOT}/out/inject_[0-9]*.csv"))
    if files:
        I = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        I["det"] = I.dchi2_match >= thr1
        Pb = np.array([0.5, 1, 2, 4, 8, 15, 30]); Rb = np.array([0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])
        grid = np.full((len(Rb) - 1, len(Pb) - 1), np.nan); ngrid = np.zeros_like(grid, dtype=int)
        for i in range(len(Rb) - 1):
            for j in range(len(Pb) - 1):
                k = (I.rp >= Rb[i]) & (I.rp < Rb[i + 1]) & (I.P >= Pb[j]) & (I.P < Pb[j + 1])
                ngrid[i, j] = k.sum()
                if k.sum(): grid[i, j] = I.det[k].mean()
        # radius of 50 % / 90 % completeness per period bin (logistic fit in log R)
        from scipy.optimize import minimize
        r50 = {}
        for j in range(len(Pb) - 1):
            k = (I.P >= Pb[j]) & (I.P < Pb[j + 1])
            x = np.log(I.rp[k].values); yv = I.det[k].values.astype(float)
            def nll(p):
                pr = 1 / (1 + np.exp(-(x - p[0]) / max(p[1], 1e-3))); pr = np.clip(pr, 1e-6, 1 - 1e-6)
                return -np.sum(yv * np.log(pr) + (1 - yv) * np.log(1 - pr))
            res = minimize(nll, [np.log(2.5), 0.15], method="Nelder-Mead")
            r50[f"{Pb[j]}-{Pb[j + 1]}"] = dict(R50=float(np.exp(res.x[0])), R90=float(np.exp(res.x[0] + 2.197 * abs(res.x[1]))), n=int(k.sum()))
        big = I.rp >= 4.5
        out["completeness"] = dict(n_injections=len(I), threshold=thr1, P_bins=Pb.tolist(), R_bins=Rb.tolist(), frac=grid.tolist(), n=ngrid.tolist(),
                                   R50_R90_by_P=r50, frac_rp_ge_4p5=float(I.det[big].mean()),
                                   frac_rp_ge_4p5_P_lt_15=float(I.det[big & (I.P < 15)].mean()),
                                   frac_best_peak_is_match_when_det=float(I.best_is_match[I.det].mean()) if I.det.any() else np.nan)
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", ["#f2f6fc", "#9cc0ee", C1, "#0d3f7a"])
        im = ax.pcolormesh(np.arange(len(Pb)), np.arange(len(Rb)), grid, cmap=cmap, vmin=0, vmax=1, edgecolors="#fcfcfb", linewidth=2)
        for i in range(len(Rb) - 1):
            for j in range(len(Pb) - 1):
                if np.isfinite(grid[i, j]):
                    ax.text(j + 0.5, i + 0.5, f"{grid[i, j]:.0%}\n(n={ngrid[i, j]})", ha="center", va="center", fontsize=7,
                            color="#fcfcfb" if grid[i, j] > 0.55 else INK)
        ax.set_xticks(np.arange(len(Pb))); ax.set_xticklabels([f"{p:g}" for p in Pb]); ax.set_yticks(np.arange(len(Rb))); ax.set_yticklabels([f"{r:g}" for r in Rb])
        ax.set_xlabel("injected period (d)"); ax.set_ylabel("injected radius (R_earth), R* = 0.381 R_sun")
        ax.set_title(f"Injection-recovery: coherent BLS, all 6 sectors, detection at global FAP 1% (dchi2 >= {thr1:.0f})", fontsize=9, color=INK)
        cb = plt.colorbar(im, ax=ax, fraction=0.04); cb.set_label("fraction recovered")
        fig.tight_layout(); fig.savefig(f"{ROOT}/figs/completeness.png", dpi=110); plt.close(fig)
    # periodogram figure
    z = np.load(f"{ROOT}/data/bls_real_all_clean.npz")
    P = 1 / z["freq"]; d2 = z["dchi2"]
    edges = np.logspace(np.log10(0.3), np.log10(30), 1500)

    def binmax(PP, yy):
        idx = np.digitize(PP, edges); o = np.full(len(edges) - 1, np.nan)
        s_ = pd.Series(yy).groupby(idx).max()
        for i, v in s_.items():
            if 1 <= i < len(edges): o[i - 1] = v
        return o
    pm = binmax(P, d2)
    zb = np.load(f"{ROOT}/data/bls_real_all_B.npz"); pmb = binmax(1 / zb["freq"], zb["dchi2"])
    fig, ax = plt.subplots(2, 1, figsize=(8, 5.6), sharex=True)
    ec = np.sqrt(edges[1:] * edges[:-1])
    ax[0].plot(ec, pmb, lw=0.8, color=C2, label="phot model-B flux as delivered")
    ax[0].axvline(0.648528, color=INK2, lw=0.6, ls=":"); ax[0].text(0.66, np.nanmax(pmb) * 0.95, "Gaia DR3 RRab 200\" away, P = 0.648528 d (and harmonics)", fontsize=7, color=INK2)
    ax[0].set_ylabel("BLS dchi2"); ax[0].legend(loc="upper right", frameon=False, fontsize=7)
    ax[1].plot(ec, pm, lw=0.8, color=C1, label="RRab leak removed (searched light curve)")
    ax[1].axhline(thr1, color=INK, lw=0.8, ls="--"); ax[1].text(12, thr1 * 1.03, f"global FAP 1% ({thr1:.0f})", fontsize=7, color=INK); ax[1].set_ylim(0, thr1 * 1.2)
    ax[1].set_xscale("log"); ax[1].set_xlabel("period (d)"); ax[1].set_ylabel("BLS dchi2"); ax[1].legend(loc="lower left", frameon=False, fontsize=7)
    for a in ax: a.grid(axis="y", color=GRID, lw=0.5)
    fig.suptitle("Coherent BLS of the deblended M3V light curve, 5.69-d transits masked (S12-S104, 2019-2026)", fontsize=9, color=INK)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/bls_clean.png", dpi=110); plt.close(fig)
    # null distribution figure
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    bins = np.logspace(np.log10(60), np.log10(max(400, max([r["max"] for v in fam.values() for r in v] + [0]) * 1.1)), 40)
    ax.hist(m3v_glob, bins=bins, color=C1, alpha=0.85, label=f"M3V rolled / inverted+rolled (n={len(m3v_glob)})", edgecolor="#fcfcfb")
    oth = [r["max"] for k in fam if k not in ("roll", "inv_roll") for r in fam[k]]
    if oth: ax.hist(oth, bins=bins, color=C2, alpha=0.7, label=f"pair, Mira, controls, phantoms (n={len(oth)})", edgecolor="#fcfcfb")
    ax.axvline(pk[0]["dchi2"], color=INK, lw=1.2); ax.text(pk[0]["dchi2"] + 3, ax.get_ylim()[1] * 0.85, f"real max {pk[0]['dchi2']:.0f}\n(P = {pk[0]['P']:.2f} d)", fontsize=7)
    ax.axvline(thr1, color=INK2, lw=0.8, ls="--"); ax.text(thr1 + 3, ax.get_ylim()[1] * 0.55, "FAP 1%", fontsize=7, color=INK2)
    ax.set_xscale("log"); ax.set_xlabel("global max BLS dchi2 (0.3-30 d), log scale"); ax.set_ylabel("realisations"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/null_dist.png", dpi=110); plt.close(fig)
    save_json(out, f"{ROOT}/out/summary.json")
    print(json.dumps({k: v for k, v in out.items() if k not in ("completeness",)}, indent=0, default=str)[:6000])
    if "completeness" in out: print(json.dumps(out["completeness"]["R50_R90_by_P"]), out["completeness"]["frac_rp_ge_4p5"])


if __name__ == "__main__":
    main()
