"""s05: pixel-level localisation of a periodic signal (default: the 0.6485-d modulation found by s01), with
injection calibration.

Per sector, from the cached difference-image crops (s04): cadences quality_ok and outside the 5.69-d windows, binned to
30 min. Per pixel: least squares of D_j(t) on a cubic B-spline in time (knots every 1.5 d; absorbs the Mira's slow
variation, scattered light and any slow trend) plus sin/cos of the first NH harmonics of the trial period.
-> 2*NH amplitude images. Per-pixel errors: the same amplitude images at 60 fake periods (0.55-0.80 d, avoiding the
trial period and its 1/2, 2x aliases by > 0.01 d) -> empirical std per pixel and component.
Spatial fit (pixels within 2.6 px of the M3V): every component image = flux_c x template + const_c, jointly
(one template, free flux per component). Templates: PRF at the M3V, pair (A+B flux-weighted), bg3 (3.1" NE),
bg5 (5.5" NNE), the Mira's empirical image, and a PRF at a free position (grid +-1.2 px, step 0.04 px).
Injections: a signal with the measured temporal pattern (amplitudes of the M3V fit) at a fake period, placed on
the M3V / pair / bg3 PRF, added to D; the recovered free position and the hypothesis chi2 differences give
the calibrated errors and the discrimination power.
Usage: s05_pix_periodic.py [period_d] [tag] [n_inject]
Outputs: out/pix_periodic_<tag>.json, figs/pix_periodic_<tag>.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from scipy.interpolate import BSpline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NH = 2
RNG = np.random.default_rng(11)


def bspline_basis(t, spacing=1.5, k=3):
    t0, t1 = t.min() - 1e-3, t.max() + 1e-3
    n = max(2, int(np.ceil((t1 - t0) / spacing)))
    inner = np.linspace(t0, t1, n + 1)
    knots = np.concatenate([[t0] * k, inner, [t1] * k])
    nb_ = len(knots) - k - 1
    B = np.zeros((len(t), nb_))
    for i in range(nb_):
        c = np.zeros(nb_); c[i] = 1
        B[:, i] = BSpline(knots, c, k, extrapolate=False)(t)
    B = np.nan_to_num(B)
    return B[:, B.sum(0) > 1e-6]


def load_sector(s, binmin=30.0):
    z = np.load(f"{ROOT}/data/pix_s{s}.npz")
    t = z["t_tdb"]; ok = z["ok"].astype(bool) & ~z["intr"].astype(bool) & np.isfinite(t)
    D = z["D"].astype(float); ok &= np.all(np.isfinite(D.reshape(len(t), -1)), 1)
    t, D = t[ok], D[ok]
    w = max(binmin * 60, float(z["exptime"]) * 86400) / 86400.
    # bin within contiguous segments
    k = np.floor((t - t[0]) / w).astype(int)
    u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    tb = np.bincount(inv, t) / n
    Db = np.zeros((len(u),) + D.shape[1:])
    np.add.at(Db, inv, D); Db /= n[:, None, None]
    keep = n >= max(1, 0.5 * np.median(n))
    meta = dict(names=list(z["pos_names"]), xy=z["pos_xy"], prf_names=list(z["prf_names"]), prfs=z["prf_imgs"], b=z["b"], jac=z["jac"])
    return tb[keep], Db[keep], meta


def amp_images(t, D, P, Bs):
    """Per-pixel LSQ: D = spline + sum_h (a_h sin + b_h cos). Returns (2NH, ny, nx) amplitude images."""
    cols = []
    for h in range(1, NH + 1):
        cols += [np.sin(2 * np.pi * h * t / P), np.cos(2 * np.pi * h * t / P)]
    X = np.column_stack([Bs] + cols)
    Y = D.reshape(len(t), -1)
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return c[Bs.shape[1]:].reshape((2 * NH,) + D.shape[1:])


def fit_template(A, E, T, m):
    """Joint fit: A_c = F_c T + k_c over mask m; returns chi2, F (2NH)."""
    chi2 = 0.0; F = []
    for c in range(A.shape[0]):
        w = 1 / E[c][m]
        X = np.column_stack([T[m], np.ones(m.sum())]) * w[:, None]
        y = A[c][m] * w
        cc, *_ = np.linalg.lstsq(X, y, rcond=None)
        chi2 += float(np.sum((X @ cc - y) ** 2)); F.append(cc[0])
    return chi2, np.array(F)


class PRFGrid:
    """PRF images at arbitrary crop positions via the phot PRF model (built lazily per sector)."""
    def __init__(self, s, x0, y0):
        sys.path.insert(0, f"{PHOT}/scripts")
        from photlib import Cutout, PRFModel
        C = Cutout(s); self.prf = PRFModel(C)
        W = json.load(open(f"{PHOT}/out/wcs_s{s}.json")); self.prf.blur = W["blur"]
        self.x0, self.y0 = x0, y0
        self.cache = {}

    def img(self, x, y, shape):
        key = (round(x, 3), round(y, 3))
        if key not in self.cache:
            full = self.prf.image(x + self.x0, y + self.y0)
            self.cache[key] = full[self.y0:self.y0 + shape[0], self.x0:self.x0 + shape[1]]
        return self.cache[key]


def localise(A, E, meta, grid, m, step=0.04, half=1.2):
    res = {}
    for nm in ("M3V", "pair", "bg3", "bg5"):
        T = meta["prfs"][meta["prf_names"].index(nm)]
        res[nm] = fit_template(A, E, T, m)
    b = meta["b"] / meta["b"].max()
    res["Mira"] = fit_template(A, E, b, m)
    xm, ym = meta["xy"][meta["names"].index("M3V")]
    best = (np.inf, None, None)
    xs = np.arange(xm - half, xm + half + 1e-9, step); ys = np.arange(ym - half, ym + half + 1e-9, step)
    chi = np.zeros((len(ys), len(xs)))
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            c2, F = fit_template(A, E, grid.img(x, y, A.shape[1:]), m)
            chi[i, j] = c2
            if c2 < best[0]: best = (c2, (x, y), F)
    res["free"] = (best[0], best[2]); res["free_xy"] = best[1]
    # offsets in arcsec (E, N) relative to the M3V
    dxy = np.array(best[1]) - np.array([xm, ym])
    res["free_dEN"] = (meta["jac"] @ dxy).tolist()
    return res, chi, xs, ys


def fake_periods(P, n):
    out = []
    while len(out) < n:
        p = RNG.uniform(0.55, 0.80)
        if all(abs(p - a) > 0.012 for a in (P, P / 2, 2 * P, P * 1.0)): out.append(p)
    return out


def main():
    P = float(sys.argv[1]) if len(sys.argv) > 1 else 0.64853
    tag = sys.argv[2] if len(sys.argv) > 2 else "p0p6485"
    ninj = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    out = {}
    fig, axs = plt.subplots(len(SECTORS), 4, figsize=(14, 3.2 * len(SECTORS)))
    for row, s in enumerate(SECTORS):
        t, D, meta = load_sector(s)
        z = np.load(f"{ROOT}/data/pix_s{s}.npz"); x0, y0 = int(z["x0"]), int(z["y0"])
        grid = PRFGrid(s, x0, y0)
        Bs = bspline_basis(t)
        A = amp_images(t, D, P, Bs)
        fakes = fake_periods(P, 60)
        AF = np.array([amp_images(t, D, pf, Bs) for pf in fakes])
        E = AF.std(0)
        ny, nx = D.shape[1:]
        yy, xx = np.mgrid[0:ny, 0:nx]
        xm, ym = meta["xy"][meta["names"].index("M3V")]
        m = np.hypot(xx - xm, yy - ym) <= 2.6
        res, chi, xs, ys = localise(A, E, meta, grid, m)
        F3 = res["M3V"][1]
        amp_e = float(np.sqrt(np.sum(F3[:2] ** 2)))      # first-harmonic semi-amplitude (e-/s) if on the M3V
        # null distribution of the M3V-fit amplitude at fake periods
        amp_fake = np.array([np.sqrt(np.sum(fit_template(a, E, meta["prfs"][meta["prf_names"].index("M3V")], m)[1][:2] ** 2)) for a in AF])
        sec = dict(n=len(t), chi2={k: float(v[0]) for k, v in res.items() if k not in ("free_xy", "free_dEN")},
                   flux_M3V_template=F3.tolist(), amp1_M3V_e_s=amp_e, amp1_fake_median=float(np.median(amp_fake)),
                   amp1_fake_p99=float(np.percentile(amp_fake, 99)), free_xy=res["free_xy"], free_dEN_arcsec=res["free_dEN"],
                   ndof=int(m.sum() * 2 * NH - 2 * NH * 2))
        # injections: same temporal pattern as measured (M3V template fluxes), fake periods, on M3V / pair / bg3
        inj = {}
        for host in ("M3V", "pair", "bg3"):
            T = meta["prfs"][meta["prf_names"].index(host)]
            rec = []
            for k in range(ninj):
                pf = fakes[k % len(fakes)] + RNG.uniform(-0.003, 0.003)
                ph = RNG.uniform(0, 2 * np.pi)
                sig = np.zeros(len(t))
                for h in range(1, NH + 1):
                    a, b_ = F3[2 * (h - 1)], F3[2 * (h - 1) + 1]
                    amp = np.hypot(a, b_); phi = np.arctan2(b_, a)
                    sig += amp * np.sin(2 * np.pi * h * t / pf + phi + h * ph)
                D2 = D + sig[:, None, None] * T[None]
                A2 = amp_images(t, D2, pf, Bs)
                r2, _, _, _ = localise(A2, E, meta, grid, m)
                rec.append(dict(dEN=r2["free_dEN"], chi2={kk: float(v[0]) for kk, v in r2.items() if kk not in ("free_xy", "free_dEN")}))
            dEN = np.array([r["dEN"] for r in rec])
            inj[host] = dict(dEN_mean=dEN.mean(0).tolist(), dEN_std=dEN.std(0).tolist(),
                             dchi2_M3V_minus_free=[r["chi2"]["M3V"] - r["chi2"]["free"] for r in rec],
                             dchi2_pair_minus_M3V=[r["chi2"]["pair"] - r["chi2"]["M3V"] for r in rec],
                             dchi2_bg3_minus_M3V=[r["chi2"]["bg3"] - r["chi2"]["M3V"] for r in rec],
                             dchi2_Mira_minus_M3V=[r["chi2"]["Mira"] - r["chi2"]["M3V"] for r in rec],
                             dEN=dEN.tolist())
        sec["inject"] = inj
        out[s] = sec
        print(f"S{s}: amp1(M3V) {amp_e:.2f} e/s (fake median {np.median(amp_fake):.2f}, p99 {np.percentile(amp_fake, 99):.2f}); chi2 " +
              ", ".join(f"{k} {v:.1f}" for k, v in sec["chi2"].items()) + f"; free dEN {res['free_dEN'][0]:+.2f},{res['free_dEN'][1]:+.2f}\"; "
              f"inj M3V std {np.round(inj['M3V']['dEN_std'], 2)} mean {np.round(inj['M3V']['dEN_mean'], 2)}", flush=True)
        # figure: first-harmonic amplitude image (projected on the measured phase), error, M3V model residual, chi2 map
        ph = np.arctan2(F3[1], F3[0])
        Aproj = A[0] * np.cos(ph) + A[1] * np.sin(ph)
        ax = axs[row]
        ax[0].imshow(Aproj, origin="lower", cmap="RdBu_r"); ax[0].set_title(f"S{s} amp image (e/s), 1st harmonic")
        for nm, (x, y) in zip(meta["names"], meta["xy"]):
            ax[0].plot(x, y, "k+" if nm != "M3V" else "g+", ms=8); ax[0].text(x + 0.1, y + 0.1, nm, fontsize=6)
        ax[0].plot(*res["free_xy"], "rx", ms=8)
        ax[1].imshow(E[0], origin="lower"); ax[1].set_title("error (fake periods)")
        T3 = meta["prfs"][meta["prf_names"].index("M3V")]
        ax[2].imshow(Aproj - np.hypot(F3[0], F3[1]) * T3, origin="lower", cmap="RdBu_r"); ax[2].set_title("resid (M3V PRF)")
        ax[3].contourf(xs, ys, chi - chi.min(), levels=[0, 2.3, 6.2, 11.8, 25, 50, 100], cmap="viridis_r")
        for nm, (x, y) in zip(meta["names"], meta["xy"]):
            if abs(x - xm) < 1.3 and abs(y - ym) < 1.3: ax[3].plot(x, y, "r+"); ax[3].text(x, y, nm, fontsize=6)
        ax[3].set_title("free-position dchi2")
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/pix_periodic_{tag}.png", dpi=70); plt.close(fig)
    # combination over sectors
    comb = {}
    use = [s for s in SECTORS if out[s]["amp1_M3V_e_s"] > out[s]["amp1_fake_p99"]]
    for k in ("M3V", "pair", "bg3", "bg5", "Mira", "free"):
        comb[f"chi2_{k}"] = float(sum(out[s]["chi2"][k] for s in SECTORS))
    dEN = np.array([out[s]["free_dEN_arcsec"] for s in SECTORS]); sd = np.array([out[s]["inject"]["M3V"]["dEN_std"] for s in SECTORS])
    bias = np.array([out[s]["inject"]["M3V"]["dEN_mean"] for s in SECTORS])
    wts = 1 / sd ** 2
    comb["free_dEN_weighted_mean_arcsec"] = ((dEN * wts).sum(0) / wts.sum(0)).tolist()
    comb["free_dEN_err_arcsec"] = (1 / np.sqrt(wts.sum(0))).tolist()
    comb["inj_M3V_bias_weighted_mean"] = ((bias * wts).sum(0) / wts.sum(0)).tolist()
    comb["sectors_signal_detected_p99"] = use
    out["combined"] = comb
    out["period_d"] = P; out["NH"] = NH
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/pix_periodic_{tag}.json")
    print("combined:", json.dumps(comb), flush=True)


if __name__ == "__main__":
    main()
