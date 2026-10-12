"""Depth in the Sector 13 QLP and eleanor-lite light curves of the 9" neighbour TIC 87055315 versus the target's own
S13 light curves, with the depth ratios expected (TESS PRF + Gaia fluxes, 2.5-px circular apertures centred on each
star) if the signal is on the target or on the neighbour. Same detrending and template scaling as elsewhere.
Writes nb_lc_check.json."""
import os, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from wotan import flatten, transit_mask
import tmodel as tm
from transit_fits import Template
from aperture_depth import gaia, TARGET_ID, NB_ID

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = "/home/user/work/tic87/data"
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th); T0, P = th[0], th[1]


def depth(t, f, E=1800.0):
    ok = np.isfinite(t) & np.isfinite(f) & (f > 0)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    mask = transit_mask(time=t, period=P, duration=0.30, T0=T0)
    ff = f.copy(); ff[mask] = np.nan
    tr = flatten(t, ff, method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
    g = np.isfinite(tr); tr = np.interp(t, t[g], tr[g]); fd = f / tr
    tau = t - T0 - np.round((t - T0) / P) * P
    w = np.abs(tau) < 0.45
    s = 1 - tmpl(tau[w], E)
    X = np.vstack([np.ones(w.sum()), -s]).T
    c, *_ = np.linalg.lstsq(X, fd[w], rcond=None)
    r = fd[w] - X @ c
    sig = 1.4826 * np.median(np.abs(fd[~mask] - np.median(fd[~mask])))
    cov = np.linalg.inv(X.T @ X) * sig ** 2
    return float(c[1] / c[0] * tmpl.depth * 1e6), float(np.sqrt(cov[1, 1]) / c[0] * tmpl.depth * 1e6)


out = {}
for tic in ("87055304", "87055315"):
    q = glob.glob(f"{DATA}/mast_dl/**/hlsp_qlp*s0013-00000000{tic}*.fits", recursive=True)[0]
    d = fits.open(q)[1].data; ok = d["QUALITY"] == 0
    out[f"qlp_sap_{tic}"] = depth(d["TIME"][ok], d["SAP_FLUX"][ok])
    e = glob.glob(f"{DATA}/mast_dl/**/hlsp_gsfc-eleanor-lite*s0013-00000000{tic}*.fits", recursive=True)[0]
    d = fits.open(e)[1].data; ok = d["QUALITY"] == 0
    for col in ("RAW_FLUX", "CORR_FLUX"):
        out[f"eleanor_{col}_{tic}"] = depth(d["TIME"][ok], d[col][ok])
# expected ratios with PRF (S13 cutout geometry)
import PRF
f = sorted(glob.glob(f"{DATA}/tesscut/*s0013*.fits"))[0]
h = fits.open(f); w = WCS(h[2].header)
col0, row0 = h[1].header["1CRV4P"] - h[1].header["1CRPX4"], h[1].header["2CRV4P"] - h[1].header["2CRPX4"]
shape = h[1].data["FLUX"].shape[1:]
tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]; nbr = gaia[gaia.source_id == NB_ID].iloc[0]
xt, yt = map(float, w.world_to_pixel_values(tgt.ra, tgt.dec)); xn, yn = map(float, w.world_to_pixel_values(nbr.ra, nbr.dec))
prf = PRF.TESS_PRF(h[0].header["CAMERA"], h[0].header["CCD"], 13, col0 + xt, row0 + yt)
yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
imgs = {}
for r in gaia.itertuples():
    x, y = w.world_to_pixel_values(r.ra, r.dec)
    if -3 < x < shape[1] + 2 and -3 < y < shape[0] + 2:
        imgs[r.source_id] = prf.locate(float(x), float(y), shape) * 10 ** (-0.4 * (r.Tmag - 10))
tot = sum(imgs.values())
for lab, (cx, cy) in (("ap_on_target", (xt, yt)), ("ap_on_neighbour", (xn, yn))):
    m = np.hypot(xx - cx, yy - cy) <= 2.5
    out[lab] = dict(frac_target=float(imgs[TARGET_ID][m].sum() / tot[m].sum()), frac_neighbour=float(imgs[NB_ID][m].sum() / tot[m].sum()))
out["pred_ratio_nbLC_over_tLC_if_target"] = out["ap_on_neighbour"]["frac_target"] / out["ap_on_target"]["frac_target"]
out["pred_ratio_nbLC_over_tLC_if_neighbour"] = out["ap_on_neighbour"]["frac_neighbour"] / out["ap_on_target"]["frac_neighbour"]
for k in ("qlp_sap", "eleanor_RAW_FLUX", "eleanor_CORR_FLUX"):
    a, ea = out[f"{k}_87055315"]; b, eb = out[f"{k}_87055304"]
    out[f"measured_ratio_{k}"] = [a / b, abs(a / b) * np.hypot(ea / a, eb / b)]
json.dump(out, open(os.path.join(HERE, "nb_lc_check.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
