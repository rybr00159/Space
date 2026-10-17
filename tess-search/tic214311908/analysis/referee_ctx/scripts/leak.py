"""Referee re-implementation of the TGLC 'leak' operator (independent code; written from tglc 0.8.1 fit_lc):

  psf_flux_k(t)      = first coefficient of a least-squares fit of [ePSF_k / sum(ePSF_k), x, y, 1] to the
                       neighbour-subtracted 5x5 cut around round(position_k)
  aperture_flux_k(t) = sum of the 3x3 pixels around round(position_k) / portion_k

Both are linear in the pixel data, so a flux change dI(pix) produces dF_k = w_k . dI.
Positions: Gaia DR3 propagated to the sector epoch, TESScut WCS, optionally plus the loc agent's normal-star WCS
calibration offset. Cut grids are aligned to the detector pixel grid: origin_k = round(x_raw_k - STAR_X_k) where
STAR_X_k is TGLC's header value (TGLC and TESScut share the SPOC FFI WCS to << 0.5 px).
"""
import glob, json
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import PRF

PRFDIR = "<local>/tic713/localise/prf_files/s0004/"
GAIA = "/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/gaia_150arcsec_all.csv"
TGLC = "tic214311908/analysis/ctx/lc/tglc"
WCSCAL = "tic214311908/analysis/loc/out/wcs_cal.json"
STARS = {"M3V": 5925373272850784768, "pairA": 5925373272850782208, "pairB": 5925373272827253632,
         "Mira": 5925373268528684288, "ref1": 5925373204131295872, "ref2": 5925373204131294336}
EXTRA = {"G17.5_3.1": 5925373272827260544, "G17.7_5.5": 5925373272827262208}
NAMES = list(STARS)


class Field:
    def __init__(self, sector, use_cal=True):
        f = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{sector:04d}-*_astrocut.fits")[0]
        self.file = f
        with fits.open(f) as h:
            self.w = WCS(h[2].header); h1 = h[1].header
            self.cam, self.ccd = h[0].header["CAMERA"], h[0].header["CCD"]
            self.col0, self.row0 = int(h1["1CRV4P"]), int(h1["2CRV4P"])
            self.tmid = float(np.nanmedian(h[1].data["TIME"]))
        self.sector = sector
        g = pd.read_csv(GAIA)
        ep = 2000.0 + (self.tmid + 2457000.0 - 2451545.0) / 365.25; dt = ep - 2016.0
        g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
        g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
        xy = self.w.all_world2pix(np.column_stack([g.ra_ep, g.dec_ep]), 0)
        g["x_raw"], g["y_raw"] = xy[:, 0], xy[:, 1]
        off = np.zeros(2)
        if use_cal:
            off = np.array(json.load(open(WCSCAL))[str(sector)]["offset_px"])
        self.off = off
        g["x"], g["y"] = g.x_raw + off[0], g.y_raw + off[1]
        self.g = g
        m = g[g.source_id == STARS["Mira"]].iloc[0]; self.xm, self.ym = m.x, m.y
        # Jacobian arcsec(E,N) per pixel
        r0 = self.w.all_pix2world([[m.x_raw, m.y_raw]], 0)[0]
        def en(p):
            r = self.w.all_pix2world([p], 0)[0]
            return np.array([(r[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
        self.J = np.column_stack([en([m.x_raw + 1, m.y_raw]), en([m.x_raw, m.y_raw + 1])]); self.Ji = np.linalg.inv(self.J)
        self.prf = PRF.TESS_PRF(self.cam, self.ccd, sector, self.col0 + self.xm, self.row0 + self.ym, localdatadir=PRFDIR)
        # TGLC cut origins
        self.cut = {}
        for k, sid in STARS.items():
            ff = glob.glob(f"{TGLC}/hlsp_tglc_tess_ffi_gaiaid-{sid}-s{sector:04d}-*.fits")
            if not ff:
                continue
            H = fits.getheader(ff[0], 0)
            r = g[g.source_id == sid].iloc[0]
            ox = r.x_raw - H["STAR_X"]; oy = r.y_raw - H["STAR_Y"]
            self.cut[k] = dict(ox_float=ox, oy_float=oy, ox=int(np.round(ox)), oy=int(np.round(oy)),
                               star_x=H["STAR_X"], star_y=H["STAR_Y"])

    def xy(self, sid):
        r = self.g[self.g.source_id == sid].iloc[0]; return float(r.x), float(r.y)

    def en_to_xy(self, de, dn):
        d = self.Ji @ np.array([de, dn]); return self.xm + d[0], self.ym + d[1]

    def xy_to_en(self, x, y):
        return self.J @ np.array([x - self.xm, y - self.ym])

    def stamp(self, x, y, shape=(21, 21)):
        """Unit-flux PRF image in TESScut pixel coordinates (integers = pixel centres)."""
        pad = 8
        big = self.prf.locate(x + pad, y + pad, (shape[0] + 2 * pad, shape[1] + 2 * pad))
        return big[pad:pad + shape[0], pad:pad + shape[1]]

    def weights(self, k, kind="psf", star_shift=(0.0, 0.0), mask=None):
        """Linear weights w_k (21x21 image) such that dF_k = sum(w_k * dI). ePSF_k at the calibrated position."""
        c = self.cut[k]; ox, oy = c["ox"], c["oy"]
        x, y = self.xy(STARS[k]); x += star_shift[0]; y += star_shift[1]
        W = np.zeros((21, 21))
        ys, xs = slice(oy, oy + 5), slice(ox, ox + 5)
        if kind == "psf":
            full = self.stamp(x, y, (21, 21))
            # TGLC normalises by the sum of the 11x11 ePSF support around round(x)
            cx, cy = int(np.round(x)), int(np.round(y))
            big = self.prf.locate(x - cx + 5, y - cy + 5, (11, 11)); s11 = big.sum()
            PX = full[ys, xs] / s11
            yy, xx = np.mgrid[0:5, 0:5] - 2.0
            A = np.column_stack([PX.ravel(), xx.ravel(), yy.ravel(), np.ones(25)])
            keep = np.ones(25, bool) if mask is None else ~mask
            pinv = np.linalg.pinv(A[keep])          # rows: params, cols: pixels kept
            wv = np.zeros(25); wv[keep] = pinv[0]
            W[ys, xs] = wv.reshape(5, 5)
        else:
            cx, cy = c["ox"] + 2, c["oy"] + 2      # round(position) in TESScut pixels
            full = self.stamp(x, y, (21, 21))
            big = self.prf.locate(x - int(np.round(x)) + 5, y - int(np.round(y)) + 5, (11, 11))
            portion = big[4:7, 4:7].sum() / big.sum()
            W[cy - 1:cy + 2, cx - 1:cx + 2] = 1.0 / portion
        return W

    def coeff_vector(self, src_xy, kind="psf", names=NAMES, Ws=None):
        img = self.stamp(*src_xy)
        Ws = Ws or {k: self.weights(k, kind) for k in names}
        return np.array([np.sum(Ws[k] * img) for k in names])


def tglc_positions(F):
    """Where TGLC itself centres each star's ePSF, in calibrated TESScut pixel coordinates.
    TGLC writes STAR_X = x_TGLC - (round(x_TGLC) - 2); x_TGLC uses Gaia DR3 positions (the HLSP apparently without
    proper-motion propagation: the high-PM M3V and pair A deviate from the PM-propagated positions by the PM shift).
    A common WCS offset between TGLC and the calibrated frame is absorbed by TGLC's ePSF, so we remove the median
    offset of the stars with negligible proper motion."""
    pos = {}
    for k, c in F.cut.items():
        pos[k] = (c["ox"] + c["star_x"], c["oy"] + c["star_y"])
    lowpm = [k for k in ("Mira", "ref1", "ref2", "pairB") if k in pos]
    d = np.median([np.array(F.xy(STARS[k])) - np.array(pos[k]) for k in lowpm], axis=0)
    return {k: (v[0] + d[0], v[1] + d[1]) for k, v in pos.items()}, d


def weights_at(F, k, xy, kind="psf"):
    """Like Field.weights but with star k's ePSF placed at an explicit position xy (cut grid from TGLC)."""
    sh = (xy[0] - F.xy(STARS[k])[0], xy[1] - F.xy(STARS[k])[1])
    return F.weights(k, kind, star_shift=sh)
