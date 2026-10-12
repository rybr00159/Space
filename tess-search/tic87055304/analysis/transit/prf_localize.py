"""PRF localisation of the transit source from the per-pixel flux-drop (dF) maps, with the PRF/WCS zero-point
calibrated on the out-of-transit image of the same frames (cross-check of the pixel_depth.py finding that the
neighbour's PRF fits the dF maps better than the target's).

Per sector (central 11x11 pixels of the TESScut cutout, background subtracted as in pixel_depth.py):
 1. Out-of-transit image O (median of the out-of-transit cadences in the transit windows). Model:
    sum_i f_i PRF(x_i + dx, y_i + dy) + b over Gaia DR3 stars (G < 19.5) inside the region (+2 px margin),
    free fluxes for the 12 brightest, one common scale for the rest, free constant b; (dx, dy) = PRF/WCS offset
    found by non-linear least squares (uniform weights; the fit is dominated by PRF-model errors).
 2. dF map (as pixel_depth.py: joint fit of all transits with per-transit linear baselines). Model:
    A PRF(xs, ys) (+ optional constant); (xs, ys) free. Source offset = (xs - dx, ys - dy) - target position,
    converted to arcsec East/North with the cutout WCS.
 3. Uncertainty: bootstrap over transits (300 resamples of the transit blocks, refitting the dF map and the
    position), plus the formal chi^2 curvature (errors scaled to chi^2_red = 1).
 4. Delta chi^2 of the dF map between the source fixed at the (offset-corrected) target and at the 9" neighbour.
Sectors are combined by inverse-variance weighting in (East, North).
Writes prf_localize.json and prf_localize.png.
"""
import os, json, glob, warnings
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tmodel as tm
from transit_fits import Template
from aperture_depth import gaia, TARGET_ID, NB_ID, RA_T, DEC_T, raw

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = "/home/user/work/tic87/data"
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th); T0, P = th[0], th[1]
W, HALF, NBOOT = 0.25, 5, 300


def load_sector(fpath):
    s = int(os.path.basename(fpath).split("-")[1][1:])
    h = fits.open(fpath); d = h[1].data; w = WCS(h[2].header)
    cam, ccd = h[0].header["CAMERA"], h[0].header["CCD"]
    col0, row0 = h[1].header["1CRV4P"] - h[1].header["1CRPX4"], h[1].header["2CRV4P"] - h[1].header["2CRPX4"]
    t = d["TIME"]; flux = d["FLUX"].astype(float); q = d["QUALITY"]
    dumps = np.unique(np.round(np.asarray(raw[s]["TIME"])[(np.asarray(raw[s]["QUALITY"]) & 32) > 0], 3))
    ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
    ok &= ~((t > 4178.80) & (t < 4179.20))
    if len(dumps):
        ok &= np.min(np.abs(t[:, None] - dumps[None, :]), axis=1) >= 1 / 24
    t, flux = t[ok], flux[ok]
    tx, ty = w.world_to_pixel_values(RA_T, DEC_T)
    med = np.median(flux, axis=0)
    yy, xx = np.mgrid[0:med.shape[0], 0:med.shape[1]]
    far = np.hypot(xx - tx, yy - ty) >= 5
    sky = far & (med < np.percentile(med[far], 35))
    flux = flux - np.median(flux[:, sky], axis=1)[:, None, None]
    E = tm.EXPT[s]
    lab = np.zeros(len(t), int); k = 0
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > 0.5 or np.any((dumps > t[i - 1]) & (dumps <= t[i])):
            k += 1
        lab[i] = k
    blocks = []
    for e in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1).astype(int):
        Tc = T0 + e * P
        m = np.abs(t - Tc) < W
        inn = m & (np.abs(t - Tc) < tmpl.T14 / 2)
        if inn.sum() == 0: continue
        vals, cnt = np.unique(lab[inn], return_counts=True)
        m &= lab == vals[np.argmax(cnt)]
        cov = min(1, np.sum(m & (np.abs(t - Tc) < tmpl.T14 / 2)) * E / 86400 / tmpl.T14)
        pre = np.sum(m & (t < Tc - tmpl.T14 / 2)) * E / 86400; post = np.sum(m & (t > Tc + tmpl.T14 / 2)) * E / 86400
        if cov >= 0.7 and pre >= 0.5 / 24 and post >= 0.5 / 24:
            blocks.append((e, np.where(m)[0], Tc))
    ix, iy = int(np.round(tx)), int(np.round(ty))
    ys, xs = slice(iy - HALF, iy + HALF + 1), slice(ix - HALF, ix + HALF + 1)
    return dict(s=s, t=t, flux=flux, w=w, cam=cam, ccd=ccd, col0=col0, row0=row0, tx=float(tx), ty=float(ty), ix=ix, iy=iy,
                ys=ys, xs=xs, blocks=blocks, E=E, shape=med.shape)


def dF_map(S, blocks):
    t, E = S["t"], S["E"]
    idx = np.concatenate([b[1] for b in blocks]); nb = len(blocks)
    X = np.zeros((len(idx), 2 * nb + 1)); r0 = 0
    oot = []
    for i, (e, ii, Tc) in enumerate(blocks):
        n = len(ii)
        X[r0:r0 + n, 2 * i] = 1; X[r0:r0 + n, 2 * i + 1] = t[ii] - Tc
        X[r0:r0 + n, -1] = -(1 - tmpl(t[ii] - Tc, E)) / tmpl.depth
        oot.append(ii[np.abs(t[ii] - Tc) > tmpl.T14 / 2 + 0.01])
        r0 += n
    Y = S["flux"][idx][:, S["ys"], S["xs"]].reshape(len(idx), -1)
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    R = Y - X @ coef
    sig = np.sqrt(np.sum(R ** 2, axis=0) / (len(idx) - X.shape[1]))
    XtXi = np.linalg.inv(X.T @ X)
    n = 2 * HALF + 1
    O = np.median(S["flux"][np.concatenate(oot)][:, S["ys"], S["xs"]], axis=0)
    return coef[-1].reshape(n, n), (sig * np.sqrt(XtXi[-1, -1])).reshape(n, n), O


class PRFModel:
    def __init__(self, S):
        import PRF
        self.prf = PRF.TESS_PRF(S["cam"], S["ccd"], S["s"], S["col0"] + S["tx"], S["row0"] + S["ty"])
        self.S = S
        stars = gaia[(gaia.phot_g_mean_mag < 19.5)].copy()
        px, py = S["w"].world_to_pixel_values(stars.ra.values, stars.dec.values)
        stars["px"], stars["py"] = px, py
        lo_x, hi_x = S["ix"] - HALF - 2, S["ix"] + HALF + 2; lo_y, hi_y = S["iy"] - HALF - 2, S["iy"] + HALF + 2
        self.stars = stars[(stars.px > lo_x) & (stars.px < hi_x) & (stars.py > lo_y) & (stars.py < hi_y)].sort_values("Tmag").reset_index(drop=True)

    def img(self, x, y):
        im = self.prf.locate(float(x), float(y), self.S["shape"])
        return im[self.S["ys"], self.S["xs"]]

    def fit_direct(self, O, nfree=12):
        st = self.stars
        def design(dx, dy):
            cols = []
            rest = np.zeros(O.shape)
            for i, r in st.iterrows():
                im = self.img(r.px + dx, r.py + dy)
                if i < nfree:
                    cols.append(im.ravel())
                else:
                    rest += im * 10 ** (-0.4 * (r.Tmag - 10))
            cols.append(rest.ravel()); cols.append(np.ones(O.size))
            return np.array(cols).T
        def chi(p):
            X = design(*p)
            c, *_ = np.linalg.lstsq(X, O.ravel(), rcond=None)
            return np.sum((O.ravel() - X @ c) ** 2)
        best = None
        for dx0 in np.arange(-0.6, 0.61, 0.2):
            for dy0 in np.arange(-0.6, 0.61, 0.2):
                v = chi((dx0, dy0))
                if best is None or v < best[0]:
                    best = (v, dx0, dy0)
        x0 = np.array([best[1], best[2]])
        simplex = np.array([x0, x0 + [0.05, 0.0], x0 + [0.0, 0.05]])
        r = minimize(chi, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3, initial_simplex=simplex))
        X = design(*r.x); c, *_ = np.linalg.lstsq(X, O.ravel(), rcond=None)
        resid = O.ravel() - X @ c
        return dict(dx=float(r.x[0]), dy=float(r.x[1]), rms_frac=float(np.std(resid) / O.max()), flux_target=float(c[0]),
                    nstars=int(len(st)), fluxes=c[:min(nfree, len(st))].tolist(), ids=st.source_id.values[:nfree].tolist(),
                    tmag=st.Tmag.values[:nfree].tolist())

    def fit_source(self, dF, edF, x0, y0, const=False):
        wgt = 1 / edF.ravel() ** 2
        def lin(x, y):
            m = self.img(x, y).ravel()
            X = np.vstack([m] + ([np.ones_like(m)] if const else [])).T
            Xw = X * wgt[:, None]
            c = np.linalg.solve(X.T @ Xw, Xw.T @ dF.ravel())
            return c, np.sum((dF.ravel() - X @ c) ** 2 * wgt)
        f = lambda p: lin(*p)[1]
        r = minimize(f, [x0, y0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-5))
        c, chi2 = lin(*r.x)
        return float(r.x[0]), float(r.x[1]), float(c[0]), float(chi2), lin


def to_EN(w, x, y, x0, y0):
    """pixel offset -> arcsec East/North at (x0, y0)."""
    ra0, dec0 = w.pixel_to_world_values(x0, y0); ra, dec = w.pixel_to_world_values(x, y)
    return float((ra - ra0) * np.cos(np.radians(dec0)) * 3600), float((dec - dec0) * 3600)


def main():
    res = dict(sectors={})
    files = sorted(glob.glob(f"{DATA}/tesscut/*.fits"))
    nbr = gaia[gaia.source_id == NB_ID].iloc[0]
    tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]
    rng = np.random.default_rng(11)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    allE = []
    for fpath in files:
        S = load_sector(fpath)
        M = PRFModel(S)
        dF, edF, O = dF_map(S, S["blocks"])
        dirf = M.fit_direct(O)
        dx, dy = dirf["dx"], dirf["dy"]
        # positions in stamp coordinates (relative to the 11x11 region origin handled inside img via full-frame coords)
        xt, yt = S["w"].world_to_pixel_values(tgt.ra, tgt.dec); xn, yn = S["w"].world_to_pixel_values(nbr.ra, nbr.dec)
        xs, ys_, A, chi2, lin = M.fit_source(dF, edF, float(xt) + dx, float(yt) + dy)
        _, chi_t = lin(float(xt) + dx, float(yt) + dy)
        _, chi_n = lin(float(xn) + dx, float(yn) + dy)
        _, chi_t0 = lin(float(xt), float(yt))
        _, chi_n0 = lin(float(xn), float(yn))
        # with free constant
        xs_c, ys_c, A_c, chi2_c, _ = M.fit_source(dF, edF, float(xt) + dx, float(yt) + dy, const=True)
        red = chi2 / (dF.size - 3)
        # bootstrap
        bx, by = [], []
        nb = len(S["blocks"])
        for b in range(NBOOT):
            pick = [S["blocks"][i] for i in rng.integers(0, nb, nb)]
            # unique epochs required for design-matrix independence: re-label duplicated blocks is fine (separate columns)
            d2, e2, _ = dF_map(S, pick)
            x2, y2, _, _, _ = M.fit_source(d2, e2, xs, ys_)
            bx.append(x2); by.append(y2)
        bx, by = np.array(bx), np.array(by)
        off_px = (xs - dx - float(xt), ys_ - dy - float(yt))
        offE, offN = to_EN(S["w"], xs - dx, ys_ - dy, float(xt), float(yt))
        bE, bN = zip(*[to_EN(S["w"], x - dx, y - dy, float(xt), float(yt)) for x, y in zip(bx, by)])
        bE, bN = np.array(bE), np.array(bN)
        sE, sN = float(np.std(bE)), float(np.std(bN))
        nE, nN = to_EN(S["w"], float(xn), float(yn), float(xt), float(yt))
        sres = dict(n_transits=nb, direct_fit=dirf, source_px=[xs, ys_], offset_px_corrected=list(map(float, off_px)),
                    offset_E_arcsec=offE, offset_N_arcsec=offN, err_E_arcsec_boot=sE, err_N_arcsec_boot=sN,
                    offset_E_arcsec_with_const=to_EN(S["w"], xs_c - dx, ys_c - dy, float(xt), float(yt))[0],
                    offset_N_arcsec_with_const=to_EN(S["w"], xs_c - dx, ys_c - dy, float(xt), float(yt))[1],
                    neighbour_E_N_arcsec=[nE, nN], chi2_free=chi2, chi2_red=red, npix=int(dF.size),
                    chi2_at_target_corr=float(chi_t), chi2_at_neighbour_corr=float(chi_n),
                    chi2_at_target_uncorr=float(chi_t0), chi2_at_neighbour_uncorr=float(chi_n0),
                    dchi2_neighbour_minus_target_corr_scaled=float((chi_n - chi_t) / red),
                    dchi2_neighbour_minus_target_uncorr_scaled=float((chi_n0 - chi_t0) / red),
                    dist_from_target_arcsec=float(np.hypot(offE, offN)),
                    dist_from_neighbour_arcsec=float(np.hypot(offE - nE, offN - nN)))
        res["sectors"][S["s"]] = sres
        allE.append((offE, offN, sE, sN, nE, nN))
        print(S["s"], "PRF/WCS offset (%.3f, %.3f) px" % (dx, dy), "source offset E %.2f +- %.2f N %.2f +- %.2f arcsec" % (offE, sE, offN, sN),
              "neighbour at (%.1f, %.1f)" % (nE, nN), "chi2_red %.2f" % red,
              "dchi2(nb - tgt) corr %.1f uncorr %.1f" % (sres["dchi2_neighbour_minus_target_corr_scaled"], sres["dchi2_neighbour_minus_target_uncorr_scaled"]), flush=True)
        c = {13: "C0", 39: "C1", 66: "C2", 93: "C3", 104: "C4"}[S["s"]]
        axes[0].errorbar(offE, offN, xerr=sE, yerr=sN, fmt="o", color=c, label=f"S{S['s']}")
        axes[0].plot(bE, bN, ",", color=c, alpha=0.3)
    A = np.array(allE)
    wE, wN = 1 / A[:, 2] ** 2, 1 / A[:, 3] ** 2
    mE, mN = np.sum(A[:, 0] * wE) / wE.sum(), np.sum(A[:, 1] * wN) / wN.sum()
    eE, eN = 1 / np.sqrt(wE.sum()), 1 / np.sqrt(wN.sum())
    chiE = float(np.sum((A[:, 0] - mE) ** 2 * wE)); chiN = float(np.sum((A[:, 1] - mN) ** 2 * wN))
    nE, nN = A[0, 4], A[0, 5]
    res["combined"] = dict(E=float(mE), E_err=float(eE), N=float(mN), N_err=float(eN), chi2_scatter_E=chiE, chi2_scatter_N=chiN, dof=len(A) - 1,
                           dist_target_arcsec=float(np.hypot(mE, mN)), sigma_from_target=float(np.hypot(mE / eE, mN / eN)),
                           sigma_from_neighbour=float(np.hypot((mE - nE) / eE, (mN - nN) / eN)), neighbour_E_N=[float(nE), float(nN)])
    # inflate if the sectors scatter more than their errors
    infl = max(1.0, np.sqrt((chiE + chiN) / (2 * (len(A) - 1))))
    res["combined"]["scatter_inflation"] = float(infl)
    res["combined"]["sigma_from_target_infl"] = float(res["combined"]["sigma_from_target"] / infl)
    res["combined"]["sigma_from_neighbour_infl"] = float(res["combined"]["sigma_from_neighbour"] / infl)
    print(res["combined"])
    json.dump(res, open(os.path.join(HERE, "prf_localize.json"), "w"), indent=1)
    ax = axes[0]
    ax.errorbar(mE, mN, xerr=eE * infl, yerr=eN * infl, fmt="s", color="k", ms=8, label="combined")
    ax.plot(0, 0, "k+", ms=15, mew=2, label="target")
    ax.plot(nE, nN, "rx", ms=12, mew=2, label='9" neighbour (G=14.0)')
    for r in gaia[(gaia.sep < 30) & (gaia.phot_g_mean_mag < 19)].itertuples():
        if r.source_id in (TARGET_ID, NB_ID): continue
        e_ = (r.ra - tgt.ra) * np.cos(np.radians(tgt.dec)) * 3600; n_ = (r.dec - tgt.dec) * 3600
        ax.plot(e_, n_, ".", color="0.5", ms=4)
    ax.set_xlim(25, -25); ax.set_ylim(-25, 25); ax.set_aspect("equal"); ax.set_xlabel("East offset (arcsec)"); ax.set_ylabel("North offset (arcsec)")
    ax.legend(fontsize=7); ax.set_title("dF-map PRF source position (PRF/WCS offset from O image removed)", fontsize=8)
    ax = axes[1]
    for s, v in res["sectors"].items():
        ax.bar(str(s), v["dchi2_neighbour_minus_target_corr_scaled"], color="C0")
    ax.axhline(0, color="k"); ax.set_ylabel("chi2(source at neighbour) - chi2(source at target), scaled")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "prf_localize.png"), dpi=120)


if __name__ == "__main__":
    main()
