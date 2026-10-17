"""s05b: full-cutout (21x21 px) search for the source of a periodic signal (default 0.64853 d) in the difference
images of every sector (s04b cache).

Per pixel: cubic B-spline (knots every 1.5 d) + NH = 3 harmonics of the trial period, on 30-min bins of quality_ok
cadences outside the 5.69-d windows. Per-pixel, per-component errors from 40 fake periods (0.55-0.80 d).
A single PRF (phot SPOC PRF model, field-star blur) + a constant per component is fitted jointly to the 2*NH component
images over all pixels, on a coarse grid of positions over the whole cutout (0.2 px), then refined (0.03 px).
The best position is converted to arcsec offsets from the M3V and compared with every Gaia star (T < 19).
Also reports, at the best position, the source's flux semi-amplitudes (e-/s) and, for each Gaia star within 1.5 px of
the best position, the fractional amplitude it would need.
Injections: the measured pattern placed at the M3V and at the best-fit position (fake periods), 10 each per sector,
to calibrate the position error.
Usage: s05b_pixfull_periodic.py [period] [tag]
Outputs: out/pixfull_periodic_<tag>.json, figs/pixfull_periodic_<tag>.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s05_pix_periodic import bspline_basis
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NH = 3
RNG = np.random.default_rng(21)


def load(s, binmin=30.0):
    z = np.load(f"{ROOT}/data/pixfull_s{s}.npz", allow_pickle=True)
    t = z["t_tdb"]; ok = z["ok"].astype(bool) & ~z["intr"].astype(bool) & np.isfinite(t)
    D = z["D"].astype(float); ok &= np.all(np.isfinite(D.reshape(len(t), -1)), 1)
    t, D = t[ok], D[ok]
    w = max(binmin * 60, float(z["exptime"]) * 86400) / 86400.
    k = np.floor((t - t[0]) / w).astype(int)
    u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    tb = np.bincount(inv, t) / n
    Db = np.zeros((len(u),) + D.shape[1:]); np.add.at(Db, inv, D); Db /= n[:, None, None]
    keep = n >= max(1, 0.5 * np.median(n))
    return tb[keep], Db[keep], z


def amp_images(t, D, P, Bs):
    cols = []
    for h in range(1, NH + 1):
        cols += [np.sin(2 * np.pi * h * t / P), np.cos(2 * np.pi * h * t / P)]
    X = np.column_stack([Bs] + cols)
    c, *_ = np.linalg.lstsq(X, D.reshape(len(t), -1), rcond=None)
    return c[Bs.shape[1]:].reshape((2 * NH,) + D.shape[1:])


class PRFs:
    def __init__(self, s):
        sys.path.insert(0, f"{PHOT}/scripts")
        from photlib import Cutout, PRFModel
        C = Cutout(s); self.prf = PRFModel(C)
        W = json.load(open(f"{PHOT}/out/wcs_s{s}.json")); self.prf.blur = W["blur"]

    def img(self, x, y):
        return self.prf.image(x, y)


def fit_all(A, E, T, m):
    """Joint: A_c = F_c T + k_c (weighted), returns chi2 and F."""
    chi2 = 0.0; F = []
    for c in range(A.shape[0]):
        w = 1 / E[c][m]
        X = np.column_stack([T[m], np.ones(m.sum())]) * w[:, None]
        cc, *_ = np.linalg.lstsq(X, A[c][m] * w, rcond=None)
        chi2 += float(np.sum((X @ cc - A[c][m] * w) ** 2)); F.append(cc[0])
    return chi2, np.array(F)


def search(A, E, prfs, m, xr, yr, step):
    best = (np.inf, None, None); grid = []
    for y in np.arange(yr[0], yr[1] + 1e-9, step):
        for x in np.arange(xr[0], xr[1] + 1e-9, step):
            c2, F = fit_all(A, E, prfs.img(x, y), m)
            grid.append((x, y, c2))
            if c2 < best[0]: best = (c2, (float(x), float(y)), F)
    return best, np.array(grid)


def locate(A, E, prfs, m, ny, nx):
    b0, grid = search(A, E, prfs, m, (0.0, nx - 1.0), (0.0, ny - 1.0), 0.25)
    x, y = b0[1]
    b1, g1 = search(A, E, prfs, m, (x - 0.3, x + 0.3), (y - 0.3, y + 0.3), 0.03)
    return b1, grid


def main():
    P = float(sys.argv[1]) if len(sys.argv) > 1 else 0.64853
    tag = sys.argv[2] if len(sys.argv) > 2 else "p0p6485"
    ninj = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    out = {}
    fig, axs = plt.subplots(2, 3, figsize=(15, 9)); axs = axs.ravel()
    for i, s in enumerate(SECTORS):
        t, D, z = load(s)
        ny, nx = D.shape[1:]
        prfs = PRFs(s)
        Bs = bspline_basis(t)
        A = amp_images(t, D, P, Bs)
        fakes = [p for p in RNG.uniform(0.55, 0.80, 60) if abs(p - P) > 0.012 and abs(p - 2 * P) > 0.012][:40]
        AF = np.array([amp_images(t, D, pf, Bs) for pf in fakes])
        E = AF.std(0) + 1e-6
        m = np.isfinite(E[0]) & (E[0] > 0)
        best, grid = locate(A, E, prfs, m, ny, nx)
        x, y = best[1]; F = best[2]
        m3 = z["m3v_xy"]; jac = z["jac"]
        dEN = (jac @ (np.array([x, y]) - m3)).tolist()
        # chi2 at catalogue positions
        sid = z["star_id"]; sx = z["star_x"]; sy = z["star_y"]; sT = z["star_T"]; sn = z["star_name"]
        dist = np.hypot(sx - x, sy - y)
        near = np.argsort(dist)[:6]
        cat = []
        for k in near:
            c2, Fk = fit_all(A, E, prfs.img(sx[k], sy[k]), m)
            amp1 = float(np.hypot(Fk[0], Fk[1]))
            cat.append(dict(source_id=int(sid[k]), name=str(sn[k]), T=float(sT[k]), dist_px=float(dist[k]),
                            dEN_from_M3V=(jac @ (np.array([sx[k], sy[k]]) - m3)).tolist(), chi2=c2, dchi2=c2 - best[0],
                            amp1_e_s=amp1, amp1_frac=amp1 / (15000 * 1.107 * 10 ** (-0.4 * (sT[k] - 10)))))
        c2m3, Fm3 = fit_all(A, E, prfs.img(*m3), m)
        # injections at the M3V (position error)
        amp_pattern = F
        inj = []
        for k in range(ninj):
            pf = fakes[k % len(fakes)] + 0.002 * RNG.standard_normal()
            sig = np.zeros(len(t)); ph0 = RNG.uniform(0, 2 * np.pi)
            for h in range(1, NH + 1):
                a, b_ = amp_pattern[2 * (h - 1)], amp_pattern[2 * (h - 1) + 1]
                sig += np.hypot(a, b_) * np.sin(2 * np.pi * h * t / pf + np.arctan2(b_, a) + h * ph0)
            D2 = D + sig[:, None, None] * prfs.img(*m3)[None]
            A2 = amp_images(t, D2, pf, Bs)
            b2, _ = search(A2, E, prfs, m, (m3[0] - 1.5, m3[0] + 1.5), (m3[1] - 1.5, m3[1] + 1.5), 0.05)
            inj.append((jac @ (np.array(b2[1]) - m3)).tolist())
        inj = np.array(inj)
        out[s] = dict(best_xy=[x, y], best_dEN_from_M3V_arcsec=dEN, chi2_best=best[0], npix=int(m.sum()), ncomp=2 * NH,
                      amp_e_s=[float(np.hypot(F[2 * h], F[2 * h + 1])) for h in range(NH)],
                      chi2_M3V=c2m3, dchi2_M3V=c2m3 - best[0], amp1_M3V_e_s=float(np.hypot(Fm3[0], Fm3[1])),
                      nearest_stars=cat, inj_M3V_dEN_mean=inj.mean(0).tolist(), inj_M3V_dEN_std=inj.std(0).tolist())
        print(f"S{s}: best ({x:.2f},{y:.2f}) = dEN {dEN[0]:+.1f},{dEN[1]:+.1f}\" from M3V; amp1 {np.hypot(F[0], F[1]):.2f} e/s; dchi2(M3V) {c2m3 - best[0]:.1f}; "
              f"inj@M3V scatter {np.round(inj.std(0), 1)}\"; nearest: " +
              "; ".join(f"{c['source_id']} {c['name']} T{c['T']:.1f} d{c['dist_px']:.2f}px dchi2 {c['dchi2']:.1f} frac {c['amp1_frac']:.2f}" for c in cat[:4]), flush=True)
        # figure: projected first-harmonic image
        ph = np.arctan2(F[1], F[0]); Ap = A[0] * np.cos(ph) + A[1] * np.sin(ph)
        ax = axs[i]; im = ax.imshow(Ap / E[0], origin="lower", cmap="RdBu_r", vmin=-6, vmax=6)
        ax.plot(*m3, "g+", ms=12, label="M3V"); ax.plot(*z["mira_xy"], "kx", ms=10, label="Mira"); ax.plot(*z["pair_xy"], "m+", ms=10, label="pair")
        ax.plot(x, y, "ro", mfc="none", ms=12, label="best")
        bright = sT < 16
        ax.plot(sx[bright], sy[bright], ".", color="0.3", ms=3)
        ax.set_title(f"S{s} 1st-harmonic amplitude / error (P = {P} d)"); ax.legend(fontsize=6)
        plt.colorbar(im, ax=ax, fraction=0.04)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/pixfull_periodic_{tag}.png", dpi=75); plt.close(fig)
    out["period_d"] = P
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/pixfull_periodic_{tag}.json")


if __name__ == "__main__":
    main()
