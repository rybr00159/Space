"""TGLC (Han & Brandt) light curves of the six Gaia stars near the Mira (Sectors 12, 39): which source is dimming?

TGLC subtracts a STATIC model of every neighbour (fluxes fixed from Gaia) and then fits the target's ePSF plus a plane
background to the 5x5 decontaminated cut (psf_flux), or sums 3x3 pixels / portion (aperture_flux). A flux drop dF on
star S therefore leaks into the TGLC light curve of every star X as  dF_X = dF * c_X(S), where c_X(S) is the
coefficient obtained by fitting [PRF_X, x, y, 1] to PRF_S on X's 5x5 cut (source code checked: tglc 0.8.1
effective_psf.fit_lc). We compute c_X(S) with the TESS PRF (PRF package, same as loc5.py) at the Gaia positions
(epoch-corrected) in the TESScut WCS, place each 5x5 cut using TGLC's own STAR_X/STAR_Y header values, and fit the
observed dF vector (e-/s) for each source hypothesis with a GLS fit whose covariance comes from depths measured at
400 common random phases (so the correlation between overlapping stars is included).

Also: a free-position grid search (where must the source be to reproduce TGLC's dF pattern?).
Outputs: out/tglc_leak.json, out/tglc_leak_grid_s{12,39}.npz
"""
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import PRF
import common as C

PRFDIR = "<local>/tic713/localise/prf_files/s0004/"
GAIA = "/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/gaia_150arcsec_all.csv"
STARS = {"M3V": 5925373272850784768, "pairA": 5925373272850782208, "pairB": 5925373272827253632,
         "Mira": 5925373268528684288, "ref1": 5925373204131295872, "ref2": 5925373204131294336}
EXPT = {12: 1800 / 86400, 39: 600 / 86400}
SECTORS = (39, 12)   # S12 TGLC light curves are dominated by leaked Mira/background systematics (see fig/tglc_raw_psf_flux.png)


def positions(sector):
    f = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{sector:04d}-*_astrocut.fits")[0]
    with fits.open(f) as h:
        w = WCS(h[2].header); h1 = h[1].header; cam, ccd = h[0].header["CAMERA"], h[0].header["CCD"]
        col0, row0 = int(h1["1CRV4P"]), int(h1["2CRV4P"]); tmid = np.nanmedian(h[1].data["TIME"])
    g = pd.read_csv(GAIA); ep = 2000.0 + (tmid + 2457000.0 - 2451545.0) / 365.25; dt = ep - 2016.0
    g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec)); g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
    pos = {}
    for k, sid in STARS.items():
        r = g[g.source_id == sid].iloc[0]; x, y = w.all_world2pix([[r.ra_ep, r.dec_ep]], 0)[0]; pos[k] = (float(x), float(y))
    # arcsec per pixel Jacobian at the Mira
    xm, ym = pos["Mira"]; r0 = w.all_pix2world([[xm, ym]], 0)[0]
    def off(p): r = w.all_pix2world([p], 0)[0]; return np.array([(r[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
    J = np.column_stack([off([xm + 1, ym]), off([xm, ym + 1])])
    prf = PRF.TESS_PRF(cam, ccd, sector, col0 + xm, row0 + ym, localdatadir=PRFDIR)
    return pos, prf, J, w


def tglc_star_xy(sector):
    out = {}
    for k, sid in STARS.items():
        f = glob.glob(f"lc/tglc/hlsp_tglc_tess_ffi_gaiaid-{sid}-s{sector:04d}-*.fits")[0]
        h = fits.getheader(f, 0); out[k] = (h["STAR_X"], h["STAR_Y"])
    return out


def leak_coeffs(pos, prf, sxy, src_xy):
    """c_X for a unit flux drop at src_xy, for every TGLC star X (psf and aperture versions)."""
    cpsf, cap = {}, {}
    yy, xx = np.mgrid[0:5, 0:5] - 2.0
    for k, (x, y) in pos.items():
        ox, oy = x - sxy[k][0], y - sxy[k][1]                      # TESScut coords of the cut's pixel (0,0) centre
        # evaluate on a 25x25 stamp with the 5x5 cut at [10:15, 10:15] (avoids edge effects of PRF.locate)
        R = prf.locate(src_xy[0] - ox + 10, src_xy[1] - oy + 10, (25, 25))[10:15, 10:15]   # residual from the source
        big = prf.locate(x - ox + 10, y - oy + 10, (25, 25)); PX = big[10:15, 10:15]; PX11 = big[7:18, 7:18].sum()
        A = np.column_stack([PX.ravel() / PX11, xx.ravel(), yy.ravel(), np.ones(25)])
        cpsf[k] = float(np.linalg.lstsq(A, R.ravel(), rcond=None)[0][0])
        cx, cy = int(round(sxy[k][0])), int(round(sxy[k][1]))
        portion = PX[cy - 1:cy + 2, cx - 1:cx + 2].sum() / PX11
        cap[k] = float(R[cy - 1:cy + 2, cx - 1:cx + 2].sum() / portion)
    return cpsf, cap


def observed(sector, col):
    """dF (e-/s) for every star at phase 0 and at 400 common random phases (for the covariance)."""
    rng = np.random.default_rng(11 + sector); ph = []
    while len(ph) < 400:
        p = rng.uniform(0, 1)
        if min(abs(p), abs(p - 0.5), abs(1 - p)) > 2.5 * C.T14 / C.P: ph.append(p)
    obs, rnd, lvl = {}, {}, {}
    for k, sid in STARS.items():
        f = glob.glob(f"lc/tglc/hlsp_tglc_tess_ffi_gaiaid-{sid}-s{sector:04d}-*.fits")[0]
        d = fits.getdata(f, 1); ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0) & np.isfinite(d[col])
        t, fl = d["time"][ok], np.array(d[col][ok], float)
        # work in e-/s: subtract a masked biweight trend (additive), so the fitted amplitude is dF directly
        mask = C.in_window(t) | C.in_window(t, phase0=0.5)
        med = np.nanmedian(fl); r = fl - C.detrend_additive(t, fl, mask, window=1.0)
        s = C.robust_std(r[~mask]); good = ~((~mask) & (np.abs(r) > 5 * s)); t, r = t[good], r[good]
        tm = C.template(t, EXPT[sector]); obs[k] = -np.sum(tm * r) / np.sum(tm * tm)
        rnd[k] = np.array([-np.sum((q := C.template(t, EXPT[sector], phase0=p)) * r) / np.sum(q * q) for p in ph])
        lvl[k] = med
    return obs, rnd, lvl


def gls(obs_vec, cov, model_vec):
    ci = np.linalg.inv(cov); a = model_vec @ ci @ model_vec; dF = model_vec @ ci @ obs_vec / a
    r = obs_vec - dF * model_vec; return dF, 1 / np.sqrt(a), float(r @ ci @ r)


def main():
    res = {}
    names = list(STARS)
    for col in ("psf_flux", "aperture_flux"):
        res[col] = {}
        tot = {}
        for sector in SECTORS:
            pos, prf, J, w = positions(sector); sxy = tglc_star_xy(sector)
            obs, rnd, lvl = observed(sector, col)
            ov = np.array([obs[k] for k in names]); R = np.array([rnd[k] for k in names]); cov = np.cov(R)
            # random-phase means (bias) and correlations
            corr = np.corrcoef(R)
            hyp = {}
            srcs = {k: pos[k] for k in ("M3V", "pairA", "pairB", "Mira")}
            for hname, sxy_src in srcs.items():
                cpsf, cap = leak_coeffs(pos, prf, sxy, sxy_src)
                c = cpsf if col == "psf_flux" else cap
                mv = np.array([c[k] for k in names]); dF, edF, chi2 = gls(ov, cov, mv)
                hyp[hname] = dict(coeffs={k: round(c[k], 4) for k in names}, dF=dF, dF_err=edF, chi2=chi2, dof=len(names) - 1)
                tot.setdefault(hname, []).append((ov, cov, mv))
            # null hypothesis: no signal anywhere
            chi0 = float(ov @ np.linalg.inv(cov) @ ov)
            # free-position grid (offsets in arcsec from the Mira, converted with the Jacobian)
            xm, ym = pos["Mira"]; Ji = np.linalg.inv(J); grid = np.arange(-45, 45.1, 1.5); chi = np.full((len(grid), len(grid)), np.nan)
            for i, dn in enumerate(grid):
                for j, de in enumerate(grid):
                    dx, dy = Ji @ np.array([de, dn]); cpsf, cap = leak_coeffs(pos, prf, sxy, (xm + dx, ym + dy))
                    c = cpsf if col == "psf_flux" else cap; mv = np.array([c[k] for k in names])
                    chi[i, j] = gls(ov, cov, mv)[2]
            i, j = np.unravel_index(np.nanargmin(chi), chi.shape)
            np.savez_compressed(f"out/tglc_leak_grid_{col}_s{sector}.npz", grid=grid, chi2=chi)
            star_off = {k: (Ji @ np.zeros(2)).tolist() for k in []}
            offs = {k: (J @ np.array([pos[k][0] - xm, pos[k][1] - ym])).round(1).tolist() for k in names}
            res[col][f"S{sector}"] = dict(observed_dF={k: round(obs[k], 2) for k in names},
                                          emp_err={k: round(float(np.sqrt(cov[n, n])), 2) for n, k in enumerate(names)},
                                          random_mean={k: round(float(rnd[k].mean()), 2) for k in names},
                                          corr=np.round(corr, 2).tolist(), star_offsets_arcsec_EN=offs,
                                          tglc_median_flux={k: round(lvl[k], 1) for k in names},
                                          hypotheses=hyp, chi2_null=chi0,
                                          grid_best_EN_arcsec=[float(grid[j]), float(grid[i])], grid_chi2_min=float(np.nanmin(chi)))
            print(col, sector, "obs dF:", {k: round(obs[k], 1) for k in names})
            print("   emp err:", {k: round(float(np.sqrt(cov[n, n])), 1) for n, k in enumerate(names)})
            for h, v in hyp.items():
                print(f"   {h:6s} dF={v['dF']:6.1f}+-{v['dF_err']:4.1f} chi2={v['chi2']:6.2f}/{v['dof']}  coeffs {v['coeffs']}")
            print(f"   null chi2 {chi0:.2f}/6; free-position best (E,N)=({grid[j]:.1f},{grid[i]:.1f})\" chi2 {np.nanmin(chi):.2f}; offsets {offs}")
        # joint over sectors with a common dF
        joint = {}
        for h, lst in tot.items():
            ci = [np.linalg.inv(c) for _, c, _ in lst]
            a = sum(m @ q @ m for (_, _, m), q in zip(lst, ci)); b = sum(m @ q @ o for (o, _, m), q in zip(lst, ci)); dF = b / a
            chi2 = sum(float((o - dF * m) @ q @ (o - dF * m)) for (o, _, m), q in zip(lst, ci))
            joint[h] = dict(dF=dF, dF_err=1 / np.sqrt(a), chi2=chi2, dof=2 * len(names) - 1)
            print(f"   JOINT {col} {h:6s} dF={dF:6.1f}+-{1/np.sqrt(a):4.1f} chi2={chi2:6.2f}/{2*len(names)-1}")
        res[col]["joint_common_dF"] = joint
    json.dump(res, open("out/tglc_leak.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
