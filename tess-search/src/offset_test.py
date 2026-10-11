"""Is the transit on the target or on a specific neighbour?

For each sector: positions propagated to the observation epoch with Gaia proper motions; difference-image centroid
(same estimator as centroid.py) converted to an on-sky offset; compared with the offset the same estimator returns for a
synthetic difference image made by a PSF at the neighbour's position (PSF width fitted on the direct image)."""
import sys, json, warnings
import numpy as np
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
from centroid import cutout, diff_image, centroid


def gauss(nx, ny, x0, y0, s):
    yy, xx = np.mgrid[0:ny, 0:nx]
    return np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2))


def sky_offset(w, x0, y0, x1, y1):
    a = w.pixel_to_world_values(x0, y0); b = w.pixel_to_world_values(x1, y1)
    return ((b[0] - a[0]) * np.cos(np.radians(a[1])) * 3600, (b[1] - a[1]) * 3600)


def test(tgt, nbr, sector, P, t0, dur, epoch_year, nboot=300):
    (ra, dec, pmra, pmde), (nra, nde, npmra, npmde) = tgt, nbr
    dt = epoch_year - 2016.0
    ra_t = ra + pmra * dt / 3.6e6 / np.cos(np.radians(dec)); de_t = dec + pmde * dt / 3.6e6
    ra_n = nra + npmra * dt / 3.6e6 / np.cos(np.radians(nde)); de_n = nde + npmde * dt / 3.6e6
    h = cutout(ra, dec, sector)
    r = diff_image(h, P, t0, dur)
    if r is None: return None
    diffs, outs, w = r
    tx, ty = w.world_to_pixel_values(ra_t, de_t); nx_, ny_ = w.world_to_pixel_values(ra_n, de_n)
    D = diffs.mean(0); O = outs.mean(0); ny, nx = D.shape
    cx, cy = centroid(D, tx, ty)
    rng = np.random.default_rng(2); bs = []
    for _ in range(nboot):
        k = rng.integers(0, len(diffs), len(diffs)); bs.append(centroid(diffs[k].mean(0), tx, ty))
    bs = np.array(bs)
    # PSF width from the direct image
    m = np.hypot(*np.meshgrid(np.arange(nx) - tx, np.arange(ny) - ty)) < 3
    def res(p):
        s, a, c = p
        return (O - (a * gauss(nx, ny, tx, ty, s) + c))[m]
    s = least_squares(res, [1.0, O.max(), np.median(O)], bounds=([0.4, 0, -np.inf], [3, np.inf, np.inf])).x[0]
    ex, ey = centroid(gauss(nx, ny, nx_, ny_, s), tx, ty)   # expected centroid if the loss were on the neighbour
    meas = sky_offset(w, tx, ty, cx, cy); exp_n = sky_offset(w, tx, ty, ex, ey)
    cov_px = np.cov(bs.T) + np.eye(2) * 0.15 ** 2
    d_t = np.array([cx - tx, cy - ty]); d_n = np.array([cx - ex, cy - ey])
    chi_t = float(d_t @ np.linalg.solve(cov_px, d_t)); chi_n = float(d_n @ np.linalg.solve(cov_px, d_n))
    return dict(sector=sector, n_transits=len(diffs), psf_sigma_px=float(s), neighbour_sep_px=float(np.hypot(nx_ - tx, ny_ - ty)),
                measured_offset_arcsec=[round(v, 2) for v in meas], expected_if_neighbour_arcsec=[round(v, 2) for v in exp_n],
                err_px=[float(np.sqrt(cov_px[0, 0])), float(np.sqrt(cov_px[1, 1]))], chi2_target=chi_t, chi2_neighbour=chi_n)


if __name__ == "__main__":
    cfg = json.loads(sys.argv[1])
    out = []
    for s, yr in cfg["sectors"]:
        out.append(test(cfg["target"], cfg["neighbour"], s, cfg["P"], cfg["t0"], cfg["dur"], yr))
    print(json.dumps(out, indent=1))
    tot_t = sum(o["chi2_target"] for o in out if o); tot_n = sum(o["chi2_neighbour"] for o in out if o)
    print(f"total chi2 (2 dof per sector): target {tot_t:.1f}, neighbour {tot_n:.1f}")
