"""Difference-image localisation with TESScut FFI cutouts.

For each sector: in-transit frames vs. out-of-transit frames just before and after each transit.
The difference image (out - in) shows where the flux was lost; its PRF-free centroid is compared with
the target position and with nearby Gaia stars. Uncertainty comes from bootstrapping over transits.
"""
import os, io, sys, zipfile, warnings, time
import numpy as np, requests
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings("ignore")
TC = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"
CACHE = os.path.join(os.path.dirname(__file__), "..", "tesscut")


def cutout(ra, dec, sector, size=15):
    os.makedirs(CACHE, exist_ok=True)
    path = f"{CACHE}/tc_{ra:.5f}_{dec:.5f}_s{sector}_{size}.fits"
    if not os.path.exists(path):
        for a in range(4):
            try:
                r = requests.get(TC, params=dict(ra=ra, dec=dec, y=size, x=size, sector=sector), timeout=600)
                r.raise_for_status()
                z = zipfile.ZipFile(io.BytesIO(r.content)); name = z.namelist()[0]
                open(path, "wb").write(z.read(name)); break
            except Exception:
                time.sleep(5 * 2 ** a)
        else:
            return None
    return fits.open(path)


def diff_image(hdul, P, t0, dur):
    d = hdul[1].data; t = d["TIME"]; flux = d["FLUX"]; q = d["QUALITY"]
    ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
    t, flux = t[ok], flux[ok]
    ep = np.round((t - t0) / P); dt = t - t0 - ep * P
    intr = np.abs(dt) < 0.4 * dur
    oot = (np.abs(dt) > 0.75 * dur) & (np.abs(dt) < 1.75 * dur)
    diffs, outs = [], []
    for e in np.unique(ep[intr]):
        mi = intr & (ep == e); mo = oot & (ep == e)
        if mi.sum() < 2 or mo.sum() < 4: continue
        # remove a local linear trend per pixel using the out-of-transit frames on both sides
        A = np.vstack([np.ones(mo.sum()), dt[mo]]).T
        coef, *_ = np.linalg.lstsq(A, flux[mo].reshape(mo.sum(), -1), rcond=None)
        pred = (np.vstack([np.ones(mi.sum()), dt[mi]]).T @ coef).reshape(mi.sum(), *flux.shape[1:])
        diffs.append(np.mean(pred - flux[mi], axis=0)); outs.append(np.mean(flux[mo], axis=0))
    if not diffs: return None
    return np.array(diffs), np.array(outs), WCS(hdul[2].header)


def centroid(img, x0, y0, r=3):
    ny, nx = img.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    m = (np.hypot(xx - x0, yy - y0) <= r)
    w = np.clip(img, 0, None) * m
    s = w.sum()
    if s <= 0: return np.nan, np.nan
    return (w * xx).sum() / s, (w * yy).sum() / s


def localise(ra, dec, sector, P, t0, dur, nboot=300, size=15):
    h = cutout(ra, dec, sector, size)
    if h is None: return None
    r = diff_image(h, P, t0, dur)
    if r is None: return None
    diffs, outs, w = r
    tx, ty = w.world_to_pixel_values(ra, dec)
    D = diffs.mean(0); O = outs.mean(0)
    cx, cy = centroid(D, tx, ty); ox, oy = centroid(O - np.median(O), tx, ty)
    rng = np.random.default_rng(0); bx, by = [], []
    for _ in range(nboot):
        k = rng.integers(0, len(diffs), len(diffs))
        x, y = centroid(diffs[k].mean(0), tx, ty); bx.append(x); by.append(y)
    ex, ey = np.nanstd(bx), np.nanstd(by)
    # offset of the flux-loss centroid from the catalogue position of the target; 0.15 px systematic floor
    dx, dy = cx - tx, cy - ty
    sx, sy = np.hypot(ex, 0.15), np.hypot(ey, 0.15)
    sig = np.hypot(dx / sx, dy / sy)
    return dict(sector=sector, n_transits=len(diffs), target_px=(float(tx), float(ty)), diff_centroid=(float(cx), float(cy)),
                direct_centroid=(float(ox), float(oy)), offset_px=(float(dx), float(dy)), offset_arcsec=float(np.hypot(dx, dy) * 21.0),
                err_px=(float(ex), float(ey)), offset_sigma=float(sig), diff_img=D, direct_img=O, wcs=w)
