"""Per-sector, odd/even and per-transit depths of TIC 87055304 with a shared (fixed) transit shape.

Template: free-density MCMC posterior medians (exposure-integrated per sector).
Model per group g: f = c_s * (1 - A_g * s(t)), s = 1 - template; A_g scales the depth; c_s a per-sector offset.
Solved as weighted linear least squares; uncertainties from the white noise (MCMC per-sector sigma) inflated by the
per-sector red-noise factor beta (vshape.json).  Depth_g = A_g * template depth (mid-transit, unsmeared).
Per-transit depths: from transit_fits_raw.json (raw SAP, local linear baseline, Tc fixed on the ephemeris), errors
scaled by the per-sector injection factor; and from the detrended flux here (constant offset per transit, beta-scaled).
Consistency: chi^2 about the weighted mean, p-values.
Also odd/even with different timing (half-period secondary-free check) and alternative: 2P fold depth.
Writes depths.json and depths.png.
"""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tmodel as tm
from transit_fits import Template

HERE = os.path.dirname(os.path.abspath(__file__))
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th)
DEPTH0 = tmpl.depth * 1e6
T0, P = th[0], th[1]
SIG = {s: mc["params"][f"sigma{s}_ppm"]["med"] * 1e-6 for s in tm.SECTORS}
vs = json.load(open(os.path.join(HERE, "vshape.json")))
BETA = {int(s): float(b) for s, b in vs["beta"].items()}
data = tm.load()


def design(groups_fn):
    """Return X, y, w for groups defined by groups_fn(sector, epoch) -> label."""
    rows_y, rows_w, cols = [], [], []
    labels = []
    blocks = []
    for s in tm.SECTORS:
        d = data[s]
        ep = np.round((d["t"] - T0) / P).astype(int)
        sx = 1 - tmpl(d["t"] - T0 - ep * P, d["exptime"])
        lab = np.array([groups_fn(s, e) for e in ep], dtype=object)
        blocks.append((s, d["f"], sx, lab))
        for l in sorted(set(lab.tolist()) - {None}, key=str):
            if l not in labels:
                labels.append(l)
    n = sum(len(b[1]) for b in blocks)
    X = np.zeros((n, len(tm.SECTORS) + len(labels)))
    y = np.zeros(n); w = np.zeros(n); i0 = 0
    for j, (s, f, sx, lab) in enumerate(blocks):
        k = slice(i0, i0 + len(f))
        X[k, j] = 1.0
        for li, l in enumerate(labels):
            m = lab == l
            col = np.zeros(len(f)); col[m] = -sx[m]
            X[k, len(tm.SECTORS) + li] = col
        y[k] = f; w[k] = 1 / (SIG[s] * BETA[s]) ** 2; i0 += len(f)
    return X, y, w, labels


def wls(X, y, w):
    Xw = X * w[:, None]
    cov = np.linalg.inv(X.T @ Xw)
    p = cov @ (Xw.T @ y)
    return p, cov


def group_fit(fn, name):
    X, y, w, labels = design(fn)
    p, cov = wls(X, y, w)
    ns = len(tm.SECTORS)
    # p[ns+i] = c * A ~ A (c ~ 1)
    A = p[ns:]; eA = np.sqrt(np.diag(cov)[ns:])
    out = {str(l): dict(depth_ppm=float(a * DEPTH0), err_ppm=float(e * DEPTH0), A=float(a), A_err=float(e)) for l, a, e in zip(labels, A, eA)}
    return out, A, eA, cov[ns:, ns:]


def consistency(v, e):
    w = 1 / e ** 2
    mu = np.sum(v * w) / np.sum(w)
    chi2 = float(np.sum(((v - mu) / e) ** 2)); dof = len(v) - 1
    return dict(mean=float(mu), mean_err=float(1 / np.sqrt(np.sum(w))), chi2=chi2, dof=dof, p=float(stats.chi2.sf(chi2, dof)))


res = dict(template_depth_ppm=DEPTH0, beta=BETA, sigma_ppm={s: SIG[s] * 1e6 for s in tm.SECTORS})
# all
allf, A, eA, _ = group_fit(lambda s, e: "all", "all")
res["combined"] = allf["all"]
# per sector
sec, A, eA, _ = group_fit(lambda s, e: f"S{s}", "sector")
res["per_sector"] = sec
res["per_sector_consistency"] = consistency(A * DEPTH0, eA * DEPTH0)
# S104 versus the rest
rest = [k for k in sec if k != "S104"]
vr = np.array([sec[k]["depth_ppm"] for k in rest]); er = np.array([sec[k]["err_ppm"] for k in rest])
c_rest = consistency(vr, er)
res["S104_vs_rest"] = dict(rest_mean=c_rest["mean"], rest_err=c_rest["mean_err"], rest_chi2=c_rest["chi2"], rest_p=c_rest["p"],
                           S104=sec["S104"]["depth_ppm"], S104_err=sec["S104"]["err_ppm"],
                           diff_sigma=float((sec["S104"]["depth_ppm"] - c_rest["mean"]) / np.hypot(sec["S104"]["err_ppm"], c_rest["mean_err"])),
                           ratio=float(sec["S104"]["depth_ppm"] / c_rest["mean"]))
# odd / even (absolute epoch parity relative to E_ref convention: epoch number from T0_INIT)
E_REF = json.load(open(os.path.join(HERE, "ml_start.json")))["E_ref"]
oe, A, eA, C = group_fit(lambda s, e: "odd" if (e + E_REF) % 2 else "even", "oddeven")
d = oe["odd"]["depth_ppm"] - oe["even"]["depth_ppm"]
ed = np.sqrt(C[0, 0] + C[1, 1] - 2 * C[0, 1]) * DEPTH0
res["odd_even"] = dict(odd=oe["odd"], even=oe["even"], diff_ppm=float(d), diff_err_ppm=float(ed), diff_sigma=float(d / ed))
# odd/even per sector (to check whether any one sector drives it)
oes = {}
for s in tm.SECTORS:
    try:
        r, A2, e2, C2 = group_fit(lambda ss, e, s=s: (("odd" if (e + E_REF) % 2 else "even") if ss == s else None), "x")
        if "odd" in r and "even" in r:
            dd = r["odd"]["depth_ppm"] - r["even"]["depth_ppm"]
            oes[f"S{s}"] = dict(odd=r["odd"]["depth_ppm"], even=r["even"]["depth_ppm"], diff_sigma=float(dd / (np.sqrt(C2[0, 0] + C2[1, 1] - 2 * C2[0, 1]) * DEPTH0)))
    except np.linalg.LinAlgError:
        pass
res["odd_even_per_sector"] = oes
# per transit (detrended, beta-scaled white noise)
pt, A, eA, _ = group_fit(lambda s, e: f"{s}:{e + E_REF}", "transit")
# keep only transits with >= 50 % in-transit coverage
keep = {}
for key, v in pt.items():
    s, e = key.split(":"); s = int(s); e = int(e)
    d_ = data[s]; tc = T0 + (e - E_REF) * P
    cov = min(1.0, np.sum(np.abs(d_["t"] - tc) < tmpl.T14 / 2) * d_["exptime"] / 86400 / tmpl.T14)
    if cov >= 0.5:
        keep[key] = dict(v, coverage=float(cov))
res["per_transit_detrended"] = keep
vv = np.array([v["depth_ppm"] for v in keep.values()]); ee = np.array([v["err_ppm"] for v in keep.values()])
res["per_transit_detrended_consistency"] = consistency(vv, ee)
res["per_transit_detrended_consistency"]["max_abs_dev_sigma"] = float(np.max(np.abs(vv - res["per_transit_detrended_consistency"]["mean"]) / ee))
# per transit from the raw-flux local fits with injection-calibrated errors
tf_path = os.path.join(HERE, "transit_fits_raw.json")
if os.path.exists(tf_path):
    tf = json.load(open(tf_path))
    inj = tf["injection"]
    rows = [r for r in tf["transits"] if r.get("used") and "depth_A" in r]
    v = np.array([r["depth_A"] for r in rows]) * tf["template_depth_ppm"]
    e = np.array([r["depth_A_err"] * max(1.0, inj[str(r["sector"])]["depth_scale_sd"]) for r in rows]) * tf["template_depth_ppm"]
    c = consistency(v, e)
    c["max_abs_dev_sigma"] = float(np.max(np.abs(v - c["mean"]) / e))
    c["n"] = len(v)
    res["per_transit_raw_injcal"] = dict(consistency=c, transits=[dict(sector=r["sector"], epoch=r["epoch"] + E_REF, depth_ppm=float(a), err_ppm=float(b))
                                                                   for r, a, b in zip(rows, v, e)])
    # odd/even from per-transit raw depths
    par = np.array([(r["epoch"] + E_REF) % 2 for r in rows])
    co = consistency(v[par == 1], e[par == 1]); ce = consistency(v[par == 0], e[par == 0])
    res["odd_even_raw_injcal"] = dict(odd=co["mean"], odd_err=co["mean_err"], even=ce["mean"], even_err=ce["mean_err"],
                                      diff_sigma=float((co["mean"] - ce["mean"]) / np.hypot(co["mean_err"], ce["mean_err"])))
    # sector means from raw per-transit
    sm = {}
    for s in tm.SECTORS:
        k = np.array([r["sector"] == s for r in rows])
        if k.sum():
            cs = consistency(v[k], e[k]); sm[f"S{s}"] = dict(depth_ppm=cs["mean"], err_ppm=cs["mean_err"], n=int(k.sum()), chi2=cs["chi2"], dof=cs["dof"])
    res["per_sector_raw_injcal"] = sm
    vv2 = np.array([x["depth_ppm"] for x in sm.values()]); ee2 = np.array([x["err_ppm"] for x in sm.values()])
    res["per_sector_raw_injcal_consistency"] = consistency(vv2, ee2)

json.dump(res, open(os.path.join(HERE, "depths.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("per_transit_detrended", "per_transit_raw_injcal")}, indent=1))

# figure
fig, ax = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 2.2]))
labs = [f"S{s}" for s in tm.SECTORS]
ax[0].errorbar(range(5), [sec[l]["depth_ppm"] for l in labs], [sec[l]["err_ppm"] for l in labs], fmt="o", color="k", label="template scaling (detrended, beta)")
if "per_sector_raw_injcal" in res:
    ax[0].errorbar(np.arange(5) + 0.15, [res["per_sector_raw_injcal"][l]["depth_ppm"] for l in labs],
                   [res["per_sector_raw_injcal"][l]["err_ppm"] for l in labs], fmt="s", color="C1", label="per-transit raw, injection errors")
ax[0].axhline(res["combined"]["depth_ppm"], color="C0", lw=0.8)
ax[0].errorbar([5.2], [oe["odd"]["depth_ppm"]], [oe["odd"]["err_ppm"]], fmt="v", color="C2", label="odd")
ax[0].errorbar([5.5], [oe["even"]["depth_ppm"]], [oe["even"]["err_ppm"]], fmt="^", color="C3", label="even")
ax[0].set_xticks(list(range(5)) + [5.35]); ax[0].set_xticklabels(labs + ["odd/even"])
ax[0].set_ylabel("depth (ppm)"); ax[0].legend(fontsize=7)
ax[0].set_title(f"per-sector chi2 = {res['per_sector_consistency']['chi2']:.1f}/4 (p = {res['per_sector_consistency']['p']:.3f})", fontsize=9)
ks = list(keep.keys())
tcs = [T0 + (int(k.split(':')[1]) - E_REF) * P for k in ks]
colors = {13: "C0", 39: "C1", 66: "C2", 93: "C3", 104: "C4"}
for i, k in enumerate(ks):
    s = int(k.split(":")[0])
    ax[1].errorbar(i, keep[k]["depth_ppm"], keep[k]["err_ppm"], fmt="o", color=colors[s])
ax[1].axhline(res["per_transit_detrended_consistency"]["mean"], color="k", lw=0.7)
ax[1].set_xticks(range(len(ks))); ax[1].set_xticklabels([k.replace(":", "\n") for k in ks], fontsize=6)
ax[1].set_title(f"per-transit depths (detrended, beta-scaled): chi2 = {res['per_transit_detrended_consistency']['chi2']:.1f}/{res['per_transit_detrended_consistency']['dof']}"
                f" (p = {res['per_transit_detrended_consistency']['p']:.3f})", fontsize=9)
ax[1].set_ylabel("depth (ppm)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "depths.png"), dpi=130)
