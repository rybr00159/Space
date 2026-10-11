"""Attribute a transit's flux loss between the target and a known neighbour by fitting the difference image
with two PSFs at fixed catalogue positions (Gaussian PSF whose width is fitted on the direct image)."""
import sys, json, warnings
import numpy as np
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
from centroid import cutout, diff_image


def gauss(nx, ny, x0, y0, s):
    yy, xx = np.mgrid[0:ny, 0:nx]
    return np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2)) / (2 * np.pi * s ** 2)


def fit_sector(ra, dec, nra, ndec, sector, P, t0, dur, nboot=200):
    h = cutout(ra, dec, sector)
    r = diff_image(h, P, t0, dur)
    if r is None: return None
    diffs, outs, w = r
    tx, ty = w.world_to_pixel_values(ra, dec); nx_, ny_ = w.world_to_pixel_values(nra, dec if ndec is None else ndec)
    D = diffs.mean(0); O = outs.mean(0); ny, nx = D.shape
    # PSF width from the direct image around the target (target dominates the flux)
    def rd(p):
        s, a, c, dx, dy = p
        return (O - (a * gauss(nx, ny, tx + dx, ty + dy, s) + c)).ravel()[np.hypot(*np.meshgrid(np.arange(nx) - tx, np.arange(ny) - ty)).ravel() < 4]
    fs = least_squares(rd, [0.9, O.max() * 5, np.median(O), 0, 0], bounds=([0.3, 0, -np.inf, -1, -1], [3, np.inf, np.inf, 1, 1]))
    sig, _, _, dx, dy = fs.x
    # two-source fit to the difference image (positions corrected by the direct-image offset dx, dy)
    mask = np.hypot(*np.meshgrid(np.arange(nx) - tx, np.arange(ny) - ty)) < 4
    A = np.vstack([gauss(nx, ny, tx + dx, ty + dy, sig)[mask], gauss(nx, ny, nx_ + dx, ny_ + dy, sig)[mask], np.ones(mask.sum())]).T
    coef, *_ = np.linalg.lstsq(A, D[mask], rcond=None)
    rng = np.random.default_rng(1); boot = []
    for _ in range(nboot):
        k = rng.integers(0, len(diffs), len(diffs))
        cb, *_ = np.linalg.lstsq(A, diffs[k].mean(0)[mask], rcond=None); boot.append(cb)
    boot = np.array(boot)
    # single-source fits: residual chi2 when all the flux loss is put at the target vs at the neighbour
    res = {}
    for name, col in (("target_only", 0), ("neighbour_only", 1)):
        A1 = A[:, [col, 2]]; c1, *_ = np.linalg.lstsq(A1, D[mask], rcond=None)
        res[name] = float(np.sum((D[mask] - A1 @ c1) ** 2))
    noise = np.std(boot @ A.T - (boot @ A.T).mean(0), axis=0)
    return dict(sector=sector, n_transits=len(diffs), psf_sigma=float(sig), sep_px=float(np.hypot(nx_ - tx, ny_ - ty)),
                flux_target=float(coef[0]), flux_target_err=float(boot[:, 0].std()),
                flux_neighbour=float(coef[1]), flux_neighbour_err=float(boot[:, 1].std()),
                chi2_target_only=res["target_only"] / np.median(noise) ** 2,
                chi2_neighbour_only=res["neighbour_only"] / np.median(noise) ** 2)


if __name__ == "__main__":
    ra, dec, nra, ndec, P, t0, dur = map(float, sys.argv[1:8]); sectors = [int(s) for s in sys.argv[8].split(",")]
    out = [fit_sector(ra, dec, nra, ndec, s, P, t0, dur) for s in sectors]
    print(json.dumps(out, indent=1))
