"""Matched-filter amplitude images and TESS-PRF localisation of a periodic signal in TESScut cutouts.

For every pixel, flux(t) = per-segment smooth baseline (cubic B-spline, knot spacing KNOT d) + jitter/background
regressors + alpha * template(t). alpha is the signal amplitude in that pixel (positive = the pixel dims during
the event). The alpha image is fitted with the TESS PRF at a free position; the fitted position is compared with
Gaia DR3 stars propagated to the sector epoch, with WCS errors removed by fitting the direct image of the target.
Conventions: pixel coords 0-based, images [row, col]; sky offsets (dE, dN) arcsec relative to the target.
"""
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import BSpline
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
CUT = "/home/user/work/tesscut5"
PRFDIR = "/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/prf_files"


def tmag_from_gaia(g, c):
    c = np.where(np.isfinite(c), c, 1.0)
    return g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473


def circle_mask(shape, xc, yc, rad):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(xx - xc, yy - yc) <= rad


class Cut:
    def __init__(self, tic, sector, cfg):
        f = glob.glob(f"{CUT}/{tic}/tess-s{sector:04d}-*_astrocut.fits")[0]
        with fits.open(f) as h:
            self.hdr0 = h[0].header; h1 = h[1].header; self.wcs = WCS(h[2].header)
            d = h[1].data
            self.t = np.array(d["TIME"], float); self.flux = np.array(d["FLUX"], float); self.q = np.array(d["QUALITY"], int)
        self.sector = sector; self.camera = int(self.hdr0["CAMERA"]); self.ccd = int(self.hdr0["CCD"])
        self.col0 = int(h1["1CRV4P"]); self.row0 = int(h1["2CRV4P"])
        self.ny, self.nx = self.flux.shape[1:]
        fin = np.all(np.isfinite(self.flux.reshape(len(self.t), -1)), axis=1) & np.isfinite(self.t)
        self.good = (self.q == 0) & fin
        for a, b in cfg.get("bad", {}).get(str(sector), []): self.good &= ~((self.t > a) & (self.t < b))
        self.exptime = float(np.nanmedian(np.diff(self.t[self.good])))
        ep = 2000.0 + (np.nanmedian(self.t[self.good]) + 2457000.0 - 2451545.0) / 365.25
        g = pd.read_csv(cfg["gaia"])
        dt = ep - 2016.0
        g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
        g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(g.ra_ep.values, g.dec_ep.values, 0)
        g["x"], g["y"] = x, y
        if "T" not in g: g["T"] = tmag_from_gaia(g.G, g.bp_rp)
        it = int(np.where(g.source_id == cfg["gaia_id"])[0][0])
        self.it = it; self.g = g
        self.xt, self.yt = float(x[it]), float(y[it])
        ra_t, de_t = g.ra_ep[it], g.dec_ep[it]
        g["dE"] = (g.ra_ep - ra_t) * np.cos(np.radians(de_t)) * 3600; g["dN"] = (g.dec_ep - de_t) * 3600
        # pixel -> sky Jacobian at the target
        r0 = np.array(self.wcs.all_pix2world([[self.xt, self.yt]], 0))[0]
        def off(r): return np.array([(r[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
        rx = np.array(self.wcs.all_pix2world([[self.xt + 0.01, self.yt]], 0))[0]
        ry = np.array(self.wcs.all_pix2world([[self.xt, self.yt + 0.01]], 0))[0]
        self.J = np.column_stack([off(rx) / 0.01, off(ry) / 0.01])

    def prf(self):
        import PRF
        sub = "s0001" if self.sector < 4 else "s0004"
        p = PRF.TESS_PRF(self.camera, self.ccd, self.sector, self.col0 + self.xt, self.row0 + self.yt, localdatadir=f"{PRFDIR}/{sub}/")
        shape = (self.ny, self.nx)
        def f(x, y):
            if x < -6 or y < -6 or x > shape[1] + 5 or y > shape[0] + 5: return np.zeros(shape)
            return p.locate(x, y, shape)
        return f


def bspline_basis(t, knot, deg=3):
    """Cubic B-spline basis over [t.min, t.max] with given knot spacing (columns sum to 1)."""
    lo, hi = t.min() - 1e-6, t.max() + 1e-6
    n = max(1, int(np.ceil((hi - lo) / knot)))
    inner = np.linspace(lo, hi, n + 1)
    k = np.concatenate([[lo] * deg, inner, [hi] * deg])
    nb = len(k) - deg - 1
    B = np.zeros((len(t), nb))
    for i in range(nb):
        c = np.zeros(nb); c[i] = 1
        B[:, i] = BSpline(k, c, deg, extrapolate=False)(t)
    return np.nan_to_num(B)


def template(t, cfg):
    P, t0 = cfg["P"], cfg["t0"]
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    if cfg["mode"] == "eclipse":
        a = np.abs(ph * P); T14 = cfg["T14"]; T23 = T14 * cfg.get("flat", 0.0)
        y = np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-6)))
        if cfg.get("sec_frac"):   # optional secondary eclipse at phase 0.5 with relative depth
            a2 = np.abs((((t - t0) / P) % 1 - 0.5) * P)
            y = y + cfg["sec_frac"] * np.where(a2 <= T23 / 2, 1.0, np.where(a2 >= T14 / 2, 0.0, (T14 / 2 - a2) / max((T14 - T23) / 2, 1e-6)))
        return y
    # Fourier template: dimming = -(sum of harmonics), normalised to unit peak-to-peak
    c = cfg["fourier"]; y = np.zeros_like(t)
    for k, (a, b) in enumerate(c, start=1):
        y += a * np.cos(2 * np.pi * k * (t - t0) / P) + b * np.sin(2 * np.pi * k * (t - t0) / P)
    return -y


def amplitude_image(C, cfg, extra_t=None):
    g = C.good
    t = C.t[g]; cube = C.flux[g].reshape(g.sum(), -1)
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
    cols = []
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < 20: continue
        B = bspline_basis(t[m], cfg.get("knot", 1.0))
        full = np.zeros((len(t), B.shape[1])); full[m] = B; cols.append(full)
    X = np.hstack(cols)
    # background and jitter regressors (per-pixel coefficients): faintest-pixel median and target centroid
    med = np.median(cube, 0); faint = np.argsort(med)[: max(10, cube.shape[1] // 8)]
    bkg = np.median(cube[:, faint], 1)
    m = circle_mask((C.ny, C.nx), C.xt, C.yt, 2.0).ravel()
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]; fs = cube[:, m].sum(1)
    cx = (cube[:, m] * xx.ravel()[m]).sum(1) / fs; cy = (cube[:, m] * yy.ravel()[m]).sum(1) / fs
    reg = [bkg - np.median(bkg), cx - np.median(cx), cy - np.median(cy)]
    tpl = template(t, cfg)
    X = np.column_stack([X] + [r for r in reg] + [-tpl])
    coef, *_ = np.linalg.lstsq(X, cube, rcond=None)
    res = cube - X @ coef
    dof = len(t) - X.shape[1]
    rms = np.sqrt((res ** 2).sum(0) / dof)
    XtX_inv = np.linalg.pinv(X.T @ X)
    sig_a = rms * np.sqrt(XtX_inv[-1, -1])
    alpha = coef[-1]
    # red-noise inflation from the residual time series binned at the signal timescale
    return alpha.reshape(C.ny, C.nx), sig_a.reshape(C.ny, C.nx), med.reshape(C.ny, C.nx), dict(n=len(t), ncols=X.shape[1])


def fit_psf(img, sig, psf, x0, y0, rad=3.0, extra=None, box=4.0):
    mask = circle_mask(img.shape, x0, y0, rad)
    w = 1 / sig[mask]; b = img[mask] * w
    def lin(p):
        cols = [psf(p[0], p[1])[mask]] + ([e[mask] for e in extra] if extra else []) + [np.ones(mask.sum())]
        A = np.column_stack(cols) * w[:, None]
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        return A, c
    def res(p):
        A, c = lin(p); return A @ c - b
    best = None
    for gx in np.linspace(-1.0, 1.0, 5):
        for gy in np.linspace(-1.0, 1.0, 5):
            p = np.array([x0 + gx, y0 + gy]); c2 = (res(p) ** 2).sum()
            if best is None or c2 < best[0]: best = (c2, p)
    r = least_squares(res, best[1], bounds=([x0 - box, y0 - box], [x0 + box, y0 + box]), diff_step=1e-4, x_scale=0.05)
    A, c = lin(r.x)
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * max(1.0, (r.fun ** 2).sum() / max(1, mask.sum() - 4 - (len(extra) if extra else 0)))
    return dict(x=r.x[0], y=r.x[1], cov=cov, amp=c[0], chi2=float((r.fun ** 2).sum()), npix=int(mask.sum()))


def localise(tic, sector, cfg, n_inj=6, seed=0):
    C = Cut(tic, sector, cfg)
    psf = C.prf()
    A, S, O, info = amplitude_image(C, cfg)
    snr_img = A / S
    # direct image: target position from a PRF fit with bright neighbours as extra components
    g = C.g
    nb = g[(np.hypot(g.x - C.xt, g.y - C.yt) < 4) & (g.index != C.it) & (g["T"] < g["T"][C.it] + 3)]
    extra = [psf(r.x, r.y) for r in nb.itertuples()]
    Od = fit_psf(O, np.sqrt(np.abs(O)) + 1, psf, C.xt, C.yt, rad=2.5, extra=extra, box=1.0)
    wcs_off = np.array([Od["x"] - C.xt, Od["y"] - C.yt])  # WCS/PRF offset of the target (pixels)
    # start the amplitude fit at the strongest positive S/N pixel and at the target; keep the better fit
    yy, xx = np.unravel_index(np.argmax(np.where(np.isfinite(snr_img), snr_img, -1e9)), snr_img.shape)
    fits_ = [fit_psf(A, S, psf, float(xx), float(yy)), fit_psf(A, S, psf, C.xt + wcs_off[0], C.yt + wcs_off[1])]
    F = min(fits_, key=lambda f: f["chi2"] / f["npix"])
    pos = np.array([F["x"], F["y"]]) - wcs_off            # source position in WCS pixel frame
    dEN = C.J @ (pos - np.array([C.xt, C.yt]))
    covEN = C.J @ F["cov"] @ C.J.T
    # injection test: fake signal of the same amplitude on the target at shifted phase
    rng = np.random.default_rng(seed); inj = []
    for k in range(n_inj):
        cfg2 = dict(cfg); shift = cfg["P"] * (0.25 + 0.5 * rng.random()) * (1 if k % 2 else -1)
        cfg2["t0"] = cfg["t0"] + shift
        C2 = C; fl0 = C.flux.copy()
        tpl = template(C.t, cfg2)
        C2.flux = fl0 - F["amp"] * tpl[:, None, None] * psf(C.xt + wcs_off[0], C.yt + wcs_off[1])[None]
        A2, S2, _, _ = amplitude_image(C2, cfg2)
        F2 = fit_psf(A2, S2, psf, C.xt + wcs_off[0], C.yt + wcs_off[1])
        C.flux = fl0
        d2 = C.J @ (np.array([F2["x"], F2["y"]]) - wcs_off - np.array([C.xt, C.yt]))
        inj.append(d2.tolist())
    inj = np.array(inj)
    # candidate stars: chi2 of the fitted position against each Gaia star within 2'
    covT = covEN + np.cov(inj.T) if len(inj) > 2 else covEN
    cand = []
    for r in g.itertuples():
        if np.hypot(r.dE, r.dN) > 120: continue
        d = dEN - np.array([r.dE, r.dN])
        cand.append(dict(source_id=int(r.source_id), sep=float(np.hypot(r.dE, r.dN)), T=float(r.T),
                         chi2=float(d @ np.linalg.solve(covT, d))))
    cand = sorted(cand, key=lambda c: c["chi2"])
    tgt = next(c for c in cand if c["source_id"] == cfg["gaia_id"])
    out = dict(tic=tic, sector=sector, camera=C.camera, ccd=C.ccd, n_cad=info["n"], amp=F["amp"], amp_snr=float(np.nanmax(snr_img)),
               fit_px=[F["x"], F["y"]], wcs_offset_px=wcs_off.tolist(), dE=float(dEN[0]), dN=float(dEN[1]),
               err_E=float(np.sqrt(covEN[0, 0])), err_N=float(np.sqrt(covEN[1, 1])), cov_EN=covEN.tolist(),
               inj_offsets=inj.tolist(), inj_bias=inj.mean(0).tolist() if len(inj) else None,
               inj_scatter=inj.std(0).tolist() if len(inj) > 1 else None,
               chi2_target=tgt["chi2"], best=cand[:6])
    return out, (A, S, O, C)


def plot(out, imgs, path, cfg):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    A, S, O, C = imgs
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
    ax[0].imshow(O, origin="lower", cmap="gray_r"); ax[0].set_title(f"S{C.sector} direct image")
    ax[1].imshow(A / S, origin="lower", cmap="RdBu_r", vmin=-np.nanmax(np.abs(A / S)), vmax=np.nanmax(np.abs(A / S)))
    ax[1].set_title("signal amplitude S/N (red = dims)")
    g = C.g; m = (g.x > -0.5) & (g.x < C.nx - 0.5) & (g.y > -0.5) & (g.y < C.ny - 0.5) & (g["T"] < g["T"][C.it] + cfg.get("dmag_plot", 6))
    for a in ax:
        a.scatter(g.x[m], g.y[m], s=np.clip(60 - 8 * (g["T"][m] - g["T"][C.it]), 4, 80), facecolors="none", edgecolors="g", lw=0.7)
        a.plot(C.xt, C.yt, "y+", ms=12, mew=2)
        a.plot(out["fit_px"][0] - out["wcs_offset_px"][0], out["fit_px"][1] - out["wcs_offset_px"][1], "mx", ms=10, mew=2)
    fig.suptitle(f"TIC {out['tic']} S{out['sector']}: fitted source offset ({out['dE']:.1f}, {out['dN']:.1f})\" ± ({out['err_E']:.1f}, {out['err_N']:.1f})")
    fig.tight_layout(); fig.savefig(path, dpi=90); plt.close(fig)


if __name__ == "__main__":
    tic, sector = int(sys.argv[1]), int(sys.argv[2])
    cfg = json.load(open(f"{HERE}/loccfg.json"))[str(tic)]
    out, imgs = localise(tic, sector, cfg)
    os.makedirs(f"{HERE}/out/{tic}/loc", exist_ok=True)
    json.dump(out, open(f"{HERE}/out/{tic}/loc/s{sector}.json", "w"), indent=1, default=float)
    plot(out, imgs, f"{HERE}/out/{tic}/loc/s{sector}.png", cfg)
    b = out["best"]
    print(f"TIC {tic} S{sector}: amp S/N {out['amp_snr']:.1f}; offset ({out['dE']:+.1f}, {out['dN']:+.1f})\" ± ({out['err_E']:.1f}, {out['err_N']:.1f}); "
          f"inj bias {np.round(out['inj_bias'], 1).tolist()} scatter {np.round(out['inj_scatter'], 1).tolist()}; chi2 target {out['chi2_target']:.1f}; "
          f"best: " + "; ".join(f"{c['source_id']} ({c['sep']:.0f}\", T {c['T']:.1f}, chi2 {c['chi2']:.1f})" for c in b[:3]))


# ---------------------------------------------------------------------------------------------------------------
# Additions: (1) phase-free Fourier mode (K harmonics, each a separate regressor; images fitted jointly with one
# position), (2) joint chi2 maps on a sky grid, summed over sectors, for weak signals.

def component_images(C, cfg):
    """Per-pixel regression with several signal columns. Returns list of (alpha_img, sigma_img)."""
    g = C.good; t = C.t[g]; cube = C.flux[g].reshape(g.sum(), -1)
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
    cols = []
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < 20: continue
        B = bspline_basis(t[m], cfg.get("knot", 1.0)); full = np.zeros((len(t), B.shape[1])); full[m] = B; cols.append(full)
    med = np.median(cube, 0); faint = np.argsort(med)[: max(10, cube.shape[1] // 8)]
    bkg = np.median(cube[:, faint], 1)
    mm = circle_mask((C.ny, C.nx), C.xt, C.yt, 2.0).ravel(); yy, xx = np.mgrid[0:C.ny, 0:C.nx]; fs = cube[:, mm].sum(1)
    cx = (cube[:, mm] * xx.ravel()[mm]).sum(1) / fs; cy = (cube[:, mm] * yy.ravel()[mm]).sum(1) / fs
    if cfg["mode"] == "fourier_free":
        sig = []
        for k in range(1, cfg.get("K", 2) + 1):
            ph = 2 * np.pi * k * (t - cfg["t0"]) / cfg["P"]; sig += [np.cos(ph), np.sin(ph)]
    else:
        sig = [-template(t, cfg)]
    ns = len(sig)
    X = np.column_stack(cols + [bkg - np.median(bkg), cx - np.median(cx), cy - np.median(cy)] + sig)
    coef, *_ = np.linalg.lstsq(X, cube, rcond=None)
    res = cube - X @ coef; dof = len(t) - X.shape[1]
    rms = np.sqrt((res ** 2).sum(0) / dof); Xi = np.linalg.pinv(X.T @ X)
    out = []
    for j in range(ns):
        idx = X.shape[1] - ns + j
        out.append((coef[idx].reshape(C.ny, C.nx), (rms * np.sqrt(Xi[idx, idx])).reshape(C.ny, C.nx)))
    return out, med.reshape(C.ny, C.nx)


def fit_psf_multi(comps, psf, x0, y0, rad=3.0, box=4.0):
    """Common free position, separate amplitude + constant per component image."""
    mask = circle_mask(comps[0][0].shape, x0, y0, rad)
    def res(p):
        P_ = psf(p[0], p[1])[mask]; r = []
        for A, S in comps:
            w = 1 / S[mask]; M = np.column_stack([P_, np.ones(mask.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(M, A[mask] * w, rcond=None); r.append(M @ c - A[mask] * w)
        return np.concatenate(r)
    best = None
    for gx in np.linspace(-1.5, 1.5, 7):
        for gy in np.linspace(-1.5, 1.5, 7):
            p = np.array([x0 + gx, y0 + gy]); c2 = (res(p) ** 2).sum()
            if best is None or c2 < best[0]: best = (c2, p)
    r = least_squares(res, best[1], bounds=([x0 - box, y0 - box], [x0 + box, y0 + box]), diff_step=1e-4, x_scale=0.05)
    J = r.jac; n = len(r.fun); cov = np.linalg.pinv(J.T @ J) * max(1.0, (r.fun ** 2).sum() / max(1, n - 2 - 2 * len(comps)))
    return dict(x=r.x[0], y=r.x[1], cov=cov, chi2=float((r.fun ** 2).sum()), n=n)


def chi2_map(comps, psf, C, wcs_off, grid, rad=5.0):
    """chi2 of [amp * PRF(position) + const] per component over a fixed mask, on a grid of sky offsets."""
    mask = circle_mask(comps[0][0].shape, C.xt + wcs_off[0], C.yt + wcs_off[1], rad)
    Ji = np.linalg.inv(C.J); out = np.zeros(len(grid))
    for i, (dE, dN) in enumerate(grid):
        px = np.array([C.xt, C.yt]) + Ji @ np.array([dE, dN]) + wcs_off
        P_ = psf(px[0], px[1])[mask]; c2 = 0.0
        for A, S in comps:
            w = 1 / S[mask]; M = np.column_stack([P_, np.ones(mask.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(M, A[mask] * w, rcond=None); c2 += ((M @ c - A[mask] * w) ** 2).sum()
        out[i] = c2
    return out


def direct_offset(C, psf, O):
    g = C.g
    nb = g[(np.hypot(g.x - C.xt, g.y - C.yt) < 4) & (g.index != C.it) & (g["T"] < g["T"][C.it] + 3)]
    extra = [psf(r.x, r.y) for r in nb.itertuples()]
    Od = fit_psf(O, np.sqrt(np.abs(O)) + 1, psf, C.xt, C.yt, rad=2.5, extra=extra, box=1.0)
    return np.array([Od["x"] - C.xt, Od["y"] - C.yt])
