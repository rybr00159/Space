"""s05c: full-cutout (21x21 px) check of where the 5.69-d transit signal lives, using the same machinery as s05b.

Motivation: s05b showed that the phot pipeline couples a variable star 200" away (Gaia DR3 RRab, P = 0.648528 d) into
the M3V light curve (via the cutout-wide background/regression terms). This tests whether the 5.69-d dip could
likewise come from a source elsewhere in the cutout.
Per pixel: cubic B-spline (knots every 1.5 d) + depth x exposure-integrated trapezoid (T14 1.6 h, flat 70 %) at the
referee ephemeris, on 10-min bins (30-min for S12) of quality_ok cadences. Errors: the same fit at 60 fake epochs
(ephemeris shifted by random phases 0.1-0.9, so that no fake window overlaps a real transit).
A single PRF + constant is fitted over all pixels (coarse 0.25 px grid over the whole cutout, refined 0.03 px).
Output: out/pixfull_transit.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s05_pix_periodic import bspline_basis
from s05b_pixfull_periodic import PRFs, search

RNG = np.random.default_rng(33)
T14, FLAT = 1.6 / 24, 0.7


def trap(t, tc_list, cad_d):
    """Exposure-integrated trapezoid (depth 1) for transit centres tc_list."""
    ss = np.linspace(-0.5, 0.5, 7) * cad_d
    out = np.zeros(len(t))
    for s_ in ss:
        tt = t + s_
        ph = ((tt - tc_list[0]) / P_B + 0.5) % 1 - 0.5
        x = np.abs(ph * P_B)
        ing = (T14 - FLAT * T14) / 2
        v = np.where(x <= FLAT * T14 / 2, 1.0, np.where(x < T14 / 2, (T14 / 2 - x) / ing, 0.0))
        out += v
    return out / len(ss)


def load(s, binmin=10.0):
    z = np.load(f"{ROOT}/data/pixfull_s{s}.npz", allow_pickle=True)
    t = z["t_tdb"]; ok = z["ok"].astype(bool) & np.isfinite(t)
    D = z["D"].astype(float); ok &= np.all(np.isfinite(D.reshape(len(t), -1)), 1)
    t, D = t[ok], D[ok]
    w = max(binmin * 60, float(z["exptime"]) * 86400) / 86400.
    k = np.floor((t - t[0]) / w).astype(int)
    u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    tb = np.bincount(inv, t) / n
    Db = np.zeros((len(u),) + D.shape[1:]); np.add.at(Db, inv, D); Db /= n[:, None, None]
    keep = n >= max(1, 0.5 * np.median(n))
    return tb[keep], Db[keep], z, w


def depth_image(t, D, Bs, tmpl):
    X = np.column_stack([Bs, -tmpl])
    c, *_ = np.linalg.lstsq(X, D.reshape(len(t), -1), rcond=None)
    return c[-1].reshape(D.shape[1:])


def fit1(A, E, T, m):
    w = 1 / E[m]
    X = np.column_stack([T[m], np.ones(m.sum())]) * w[:, None]
    c, *_ = np.linalg.lstsq(X, A[m] * w, rcond=None)
    return float(np.sum((X @ c - A[m] * w) ** 2)), c[0]


def main():
    out = {}
    for s in SECTORS:
        t, D, z, cad = load(s)
        prfs = PRFs(s)
        Bs = bspline_basis(t)
        tm = trap(t, [T0_B], cad)
        A = depth_image(t, D, Bs, tm)
        AF = []
        for k in range(60):
            ph = RNG.uniform(0.1, 0.9)
            tf = trap(t, [T0_B + ph * P_B], cad)
            AF.append(depth_image(t, D, Bs, tf))
        E = np.std(AF, 0) + 1e-6
        m = np.ones(A.shape, bool)
        # whole-cutout scan with one PRF (single 'component')
        A3, E3 = A[None], E[None]
        b0, grid = search(A3, E3, prfs, m, (0.0, A.shape[1] - 1.0), (0.0, A.shape[0] - 1.0), 0.25)
        x, y = b0[1]
        b1, _ = search(A3, E3, prfs, m, (x - 0.3, x + 0.3), (y - 0.3, y + 0.3), 0.03)
        x, y = b1[1]; F = b1[2][0]
        m3 = z["m3v_xy"]; jac = z["jac"]
        c2m3, F3 = fit1(A, E, prfs.img(*m3), m)
        # second-best local minimum farther than 2 px from the best (any other source in the cutout?)
        far = np.hypot(grid[:, 0] - x, grid[:, 1] - y) > 2.0
        j2 = np.argmin(np.where(far, grid[:, 2], np.inf))
        out[s] = dict(best_xy=[x, y], best_dEN_from_M3V_arcsec=(jac @ (np.array([x, y]) - m3)).tolist(), flux_loss_e_s=float(F),
                      chi2_best=b1[0], chi2_M3V=c2m3, flux_loss_M3V_e_s=float(F3), npix=int(m.sum()),
                      second_min_xy=[float(grid[j2, 0]), float(grid[j2, 1])], second_min_dchi2=float(grid[j2, 2] - b1[0]),
                      chi2_noPRF=float(np.sum(((A - np.sum(A / E ** 2) / np.sum(1 / E ** 2)) / E) ** 2)))
        print(f"S{s}: best dEN {out[s]['best_dEN_from_M3V_arcsec'][0]:+.1f},{out[s]['best_dEN_from_M3V_arcsec'][1]:+.1f}\" flux loss {F:.1f} e/s; "
              f"chi2 best {b1[0]:.1f} M3V {c2m3:.1f} ({F3:.1f} e/s) noPRF {out[s]['chi2_noPRF']:.1f} / {m.sum()} px; "
              f"2nd min (>2 px away) dchi2 {out[s]['second_min_dchi2']:.1f}", flush=True)
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/pixfull_transit.json")


if __name__ == "__main__":
    main()
