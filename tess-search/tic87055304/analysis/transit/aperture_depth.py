"""Homogeneous aperture photometry from the TESScut FFI cutouts: transit depth versus aperture, all sectors.

Purpose: (1) check whether the shallower QLP Sector 104 depth comes from its different aperture/pipeline
(QLP S104 = TGLC-based '3x3 square'; S13-S93 = classic QLP 2.5-px circle); (2) compare the depth-vs-aperture
trend with the dilution expected if the signal is on the target or on the 9" neighbour (TIC 87055315), using
the TESS PRF (Bell & Higgins 2022 'PRF' package) and Gaia DR3 fluxes.

Photometry (identical in every sector): SPOC-calibrated FFI cutouts (TESScut, 21x21 px); QUALITY == 0;
S104 scattered-light event 4178.80-4179.20 and +-1 h around momentum dumps removed; per-frame background = median of
'sky' pixels (>= 5 px from the target and in the faintest 35 % of the median image), subtracted from every pixel.
Apertures: 1x1, 2x2 (brightest 2x2 containing the target), 3x3 (QLP-S104-like), circles r = 1.5, 2.0, 2.5 (QLP S13-S93-like),
3.0 px around the target position, 5x5.
Detrending identical to prep.py (transit-masked biweight, window 1 d, segments split at gaps and dumps).
Depth: template scaling (free-density MCMC shape, exposure-integrated) with a per-sector offset, data within +-0.45 d
of transits; error = white (MAD) x beta (time-averaging, 10-120 min).
Predicted dilution: for each aperture, F_target,ap / sum_i F_i,ap from PRF-weighted Gaia stars (T from Gaia
G, BP-RP via Stassun et al. 2019), and the same for the 9" neighbour as the host.
Writes aperture_depth.json, aperture_depth.png.
"""
import os, json, glob, warnings
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wotan import flatten, transit_mask
import tmodel as tm
from transit_fits import Template
from vshape import beta_time_avg

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = "/home/user/work/tic87/data"
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th); T0, P = th[0], th[1]
DEPTH0 = tmpl.depth * 1e6
RA_T, DEC_T = 272.60647748, -43.99546646   # Gaia DR3 (epoch 2016)
raw = np.load(f"{DATA}/qlp_raw_columns.npy", allow_pickle=True).item()
gaia = pd.read_csv(f"{DATA}/gaia_dr3_90arcsec.csv")


def tess_mag(r):
    g = r.phot_g_mean_mag
    c = r.phot_bp_mean_mag - r.phot_rp_mean_mag
    if np.isfinite(c):
        return g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return g - 0.430


gaia["Tmag"] = [tess_mag(r) for r in gaia.itertuples()]
gaia = gaia[np.isfinite(gaia.Tmag)].copy()
it = int(np.argmin(gaia.sep.values))
TARGET_ID = gaia.iloc[it].source_id
nb = gaia[(gaia.sep > 7) & (gaia.sep < 11) & (gaia.phot_g_mean_mag < 15)]
NB_ID = nb.iloc[0].source_id if len(nb) else None


def apertures(tx, ty, shape, med):
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    ix, iy = int(np.round(tx)), int(np.round(ty))
    ap = {}
    m = np.zeros(shape, bool); m[iy, ix] = True; ap["1x1"] = m
    # brightest 2x2 containing target pixel
    best = None
    for dx in (-1, 0):
        for dy in (-1, 0):
            m = np.zeros(shape, bool); m[iy + dy:iy + dy + 2, ix + dx:ix + dx + 2] = True
            v = med[m].sum()
            if best is None or v > best[0]:
                best = (v, m)
    ap["2x2"] = best[1]
    m = np.zeros(shape, bool); m[iy - 1:iy + 2, ix - 1:ix + 2] = True; ap["3x3"] = m
    for r in (1.5, 2.0, 2.5, 3.0):
        ap[f"r{r}"] = np.hypot(xx - tx, yy - ty) <= r
    m = np.zeros(shape, bool); m[iy - 2:iy + 3, ix - 2:ix + 3] = True; ap["5x5"] = m
    return ap


def detrend(t, f, dumps):
    lab = np.zeros(len(t), int); k = 0
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > 0.5 or np.any((dumps > t[i - 1]) & (dumps <= t[i])):
            k += 1
        lab[i] = k
    mask = transit_mask(time=t, period=P, duration=0.30, T0=T0)
    trend = np.full_like(f, np.nan)
    for k in np.unique(lab):
        m = lab == k
        ff = f[m].copy(); ff[mask[m]] = np.nan
        if np.isfinite(ff).sum() < 5:
            continue
        tr = flatten(t[m], ff, method="biweight", window_length=1.0, return_trend=True, edge_cutoff=0)[1]
        g = np.isfinite(tr)
        if g.sum() >= 2:
            trend[m] = np.interp(t[m], t[m][g], tr[g])
    return f / trend, mask


def depth_fit(t, f, E, mask):
    ep = np.round((t - T0) / P); tau = t - T0 - ep * P
    w = (np.abs(tau) < 0.45) & np.isfinite(f)
    s = 1 - tmpl(tau[w], E)
    X = np.vstack([np.ones(w.sum()), -s]).T
    c, *_ = np.linalg.lstsq(X, f[w], rcond=None)
    r = f[w] - X @ c
    oot = ~mask & np.isfinite(f)
    sig = 1.4826 * np.median(np.abs(f[oot] - np.median(f[oot])))
    beta, _ = beta_time_avg(r, t[w], E)
    cov = np.linalg.inv(X.T @ X) * sig ** 2
    A = c[1] / c[0]; eA = np.sqrt(cov[1, 1]) / c[0]
    return A * DEPTH0, eA * DEPTH0 * beta, sig, beta


def prf_model(sector, cam, ccd, col0, row0, w, shape):
    """Return dict source_id -> PRF image (sum 1 over infinite plane, truncated to the stamp)."""
    import PRF
    prf = PRF.TESS_PRF(cam, ccd, sector, col0 + shape[1] / 2, row0 + shape[0] / 2,
                       localdatadir=os.path.join(HERE, "prf_cache") if os.path.exists(os.path.join(HERE, "prf_cache", f"cam{cam}_ccd{ccd}")) else None)
    out = {}
    sel = gaia[gaia.sep < 120]
    for r in sel.itertuples():
        x, y = w.world_to_pixel_values(r.ra, r.dec)
        if -3 < x < shape[1] + 2 and -3 < y < shape[0] + 2:
            # PRF.locate expects positions relative to the stamp, integers at pixel centres
            img = prf.locate(float(x), float(y), shape)
            out[r.source_id] = (img, 10 ** (-0.4 * (r.Tmag - 10)))
    return out


if __name__ == "__main__":
    res = dict(template_depth_ppm=DEPTH0, target_source_id=int(TARGET_ID), neighbour_source_id=int(NB_ID) if NB_ID else None, sectors={})
    files = sorted(glob.glob(f"{DATA}/tesscut/*.fits"))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for fpath in files:
        s = int(os.path.basename(fpath).split("-")[1][1:])
        h = fits.open(fpath)
        d = h[1].data; w = WCS(h[2].header)
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
        bkg = np.median(flux[:, sky], axis=1)
        fs = flux - bkg[:, None, None]
        meds = np.median(fs, axis=0)
        aps = apertures(float(tx), float(ty), med.shape, meds)
        try:
            prfs = prf_model(s, cam, ccd, col0, row0, w, med.shape)
        except Exception as e:
            print("PRF failed", s, e); prfs = None
        E = tm.EXPT[s]
        sres = dict(target_px=[float(tx), float(ty)], n_cad=int(len(t)), bkg_med_per_px=float(np.median(bkg)), apertures={})
        for name, m in aps.items():
            fl = fs[:, m].sum(axis=1)
            fn = fl / np.median(fl)
            fd, mask = detrend(t, fn, dumps)
            dep, edep, sig, beta = depth_fit(t, fd, E, mask)
            row = dict(npix=int(m.sum()), flux_med=float(np.median(fl)), bkg_frac=float(np.median(bkg) * m.sum() / np.median(fl)),
                       depth_ppm=float(dep), err_ppm=float(edep), sigma_ppm=float(sig * 1e6), beta=float(beta))
            if prfs is not None and TARGET_ID in prfs:
                tot = sum(img[m].sum() * F for img, F in prfs.values())
                ft = prfs[TARGET_ID][0][m].sum() * prfs[TARGET_ID][1]
                row["frac_target"] = float(ft / tot)
                if NB_ID in prfs:
                    row["frac_neighbour"] = float(prfs[NB_ID][0][m].sum() * prfs[NB_ID][1] / tot)
                row["frac_target_ap_of_psf"] = float(prfs[TARGET_ID][0][m].sum())
            sres["apertures"][name] = row
        res["sectors"][s] = sres
        print(s, {k: (v["npix"], round(v["depth_ppm"]), round(v["err_ppm"]), round(v.get("frac_target", np.nan), 3), round(v.get("frac_neighbour", np.nan), 4)) for k, v in sres["apertures"].items()}, flush=True)
    # hypothesis test: depth_ap = delta * frac_host,ap; fit delta per sector for host = target / neighbour; chi2
    for s, sres in res["sectors"].items():
        ap = sres["apertures"]
        names = [n for n in ap if "frac_target" in ap[n]]
        if not names:
            continue
        dv = np.array([ap[n]["depth_ppm"] for n in names]); ev = np.array([ap[n]["err_ppm"] for n in names])
        hyp = {}
        for host in ("frac_target", "frac_neighbour"):
            if not all(host in ap[n] for n in names):
                continue
            fr = np.array([ap[n][host] for n in names])
            # apertures share most pixels -> errors strongly correlated; use the chi2 only as a shape indicator,
            # and report the implied intrinsic depth on the host
            dlt = np.sum(dv * fr / ev ** 2) / np.sum(fr ** 2 / ev ** 2)
            chi2 = float(np.sum(((dv - dlt * fr) / ev) ** 2))
            hyp[host] = dict(intrinsic_depth_ppm=float(dlt), chi2=chi2, n=len(names))
        sres["hypotheses"] = hyp
        # slope: depth(5x5) / depth(1x1) and r1.5 vs r3.0
        sres["ratio_5x5_over_1x1"] = float(ap["5x5"]["depth_ppm"] / ap["1x1"]["depth_ppm"])
        if "frac_target" in ap["5x5"]:
            sres["pred_ratio_target"] = float(ap["5x5"]["frac_target"] / ap["1x1"]["frac_target"])
            if "frac_neighbour" in ap["5x5"]:
                sres["pred_ratio_neighbour"] = float(ap["5x5"]["frac_neighbour"] / ap["1x1"]["frac_neighbour"])
    json.dump(res, open(os.path.join(HERE, "aperture_depth.json"), "w"), indent=1)
    # figure
    cols = {13: "C0", 39: "C1", 66: "C2", 93: "C3", 104: "C4"}
    order = ["1x1", "2x2", "r1.5", "3x3", "r2.0", "r2.5", "r3.0", "5x5"]
    for s, sres in res["sectors"].items():
        ap = sres["apertures"]
        x = np.arange(len(order)) + (list(cols).index(int(s)) - 2) * 0.08
        axes[0].errorbar(x, [ap[n]["depth_ppm"] for n in order], [ap[n]["err_ppm"] for n in order], fmt="o-", color=cols[int(s)], ms=3, lw=0.8, label=f"S{s}")
        if "frac_target" in ap["1x1"]:
            dl = sres["hypotheses"]["frac_target"]["intrinsic_depth_ppm"]
            axes[1].plot(np.arange(len(order)), [ap[n]["depth_ppm"] / dl for n in order], "o", color=cols[int(s)], ms=4)
            axes[1].plot(np.arange(len(order)), [ap[n]["frac_target"] for n in order], "-", color=cols[int(s)], lw=1)
            if "frac_neighbour" in ap["1x1"] and "frac_neighbour" in sres["hypotheses"]:
                dn = sres["hypotheses"]["frac_neighbour"]["intrinsic_depth_ppm"]
                axes[1].plot(np.arange(len(order)), [ap[n]["frac_neighbour"] * dn / dl for n in order], ":", color=cols[int(s)], lw=1)
    axes[0].set_xticks(range(len(order))); axes[0].set_xticklabels(order); axes[0].set_ylabel("depth (ppm, own TESScut photometry)")
    axes[0].legend(fontsize=7); axes[0].set_title("depth vs aperture (identical photometry in all sectors)", fontsize=9)
    axes[1].set_xticks(range(len(order))); axes[1].set_xticklabels(order)
    axes[1].set_title("points: depth / fitted intrinsic depth; solid: predicted target flux fraction;\n dotted: prediction if on 9\" neighbour (scaled)", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "aperture_depth.png"), dpi=130)
