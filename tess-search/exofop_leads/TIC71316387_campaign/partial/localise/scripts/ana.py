"""Stacking and fitting of per-window difference images (real and fake), shared by s02 (real) and s03 (injections)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loclib import *


class Ctx:
    def __init__(self, sec, psfkind=None, tag=None):
        tag = os.environ.get('LOC_TAG', '') if tag is None else tag
        self.sec = sec
        z = np.load(f'{L}/results/win_s{sec}{tag}.npz', allow_pickle=True)
        self.z = {k: z[k] for k in z.files}
        self.jac = self.z['jac']; self.jinv = np.linalg.inv(self.jac)
        self.S = Sector(sec) if False else None
        # light-weight sector geometry (no cube): rebuild Gaia + WCS
        with fits.open(cutout_file(sec)) as h:
            hdr0, hdr1, hdr2 = h[0].header, h[1].header, h[2].header
            ny, nx = h[1].data['FLUX'].shape[1:]
        self.ny, self.nx = ny, nx
        self.xt0, self.yt0 = float(self.z['xt_wcs']), float(self.z['yt_wcs'])
        g = load_gaia()
        ra, dec = propagate(g, float(self.z['epoch_year']))
        w = WCS(hdr2)
        x, y = w.all_world2pix(ra, dec, 0)
        it = int(np.where(g.source_id == TARGET_ID)[0][0])
        self.gaia = g; self.gx = x; self.gy = y; self.it = it
        g['de'] = (ra - ra[it]) * np.cos(np.deg2rad(dec[it])) * 3600; g['dn'] = (dec - dec[it]) * 3600
        g['sep_t'] = np.hypot(g.de, g.dn)
        psfkind = os.environ.get('LOC_PSF', 'prf') if psfkind is None else psfkind
        self.psfkind = psfkind

        class _S: pass
        s_ = _S(); s_.sector = sec; s_.camera = int(hdr0['CAMERA']); s_.ccd = int(hdr0['CCD'])
        s_.col0 = int(hdr1['1CRV4P']); s_.row0 = int(hdr1['2CRV4P']); s_.ny = ny; s_.nx = nx
        s_.xt_wcs = self.xt0; s_.yt_wcs = self.yt0
        self.Sg = s_
        if psfkind == 'prf':
            self.psf = TessPRF(s_)
        else:
            sx = float(os.environ.get('LOC_GSIG', 0.0))
            self.psf = None  # set after width fit
        self.mask = circle_mask((ny, nx), self.xt0, self.yt0, RFIT)
        # other stars (tied fluxes): T < 18.5 within RFIT + 4 px of the target
        d_pix = np.hypot(x - self.xt0, y - self.yt0)
        sel = np.where((g.tmag.values < 18.5) & (d_pix < RFIT + 4.0) & (np.arange(len(g)) != it))[0]
        self.others = [(float(x[i] - self.xt0), float(y[i] - self.yt0),
                        float(10 ** (-0.4 * (g.tmag.values[i] - g.tmag.values[it]))), int(g.source_id.values[i]))
                       for i in sel]
        self.kind = self.z['kind']

    def set_gauss(self, sx, sy):
        self.psf = GaussPSF(self.Sg, sx, sy)

    def others_img(self, x, y):
        m = np.zeros((self.ny, self.nx))
        for dx, dy, r, _ in self.others:
            m += r * self.psf(x + dx, y + dy)
        return m

    # ------------------------------------------------------------------ stacking
    def stack(self, idx, ampnorm):
        """idx: window indices; ampnorm: per-window relative amplitude (1 for deep, r for shallow).
        Returns D_stack (in units of the deep peak amplitude), sigma_stack, B_stack."""
        D = self.z['D'][idx]; sg = self.z['sig'][idx]; B = self.z['B'][idx]; sb = self.z['sbar'][idx]
        a = np.asarray(ampnorm, float) * sb
        m = self.mask
        wk = a ** 2 / np.array([np.median(s[m] ** 2) for s in sg])
        Dn = D / a[:, None, None]
        W = wk.sum()
        Ds = (wk[:, None, None] * Dn).sum(0) / W
        ss = np.sqrt((wk[:, None, None] ** 2 * (sg / a[:, None, None]) ** 2).sum(0)) / W
        wb = 1.0 / np.array([np.median(s[m] ** 2) for s in sg])
        Bs = (wb[:, None, None] * B).sum(0) / wb.sum()
        return Ds, ss, Bs

    # ------------------------------------------------------------------ fits
    def fit_direct(self, B, sig):
        """Target position in the direct (baseline) image: PSF(x,y)*F_T + tied others + const."""
        r = fit_free(B, sig, self.psf, self.xt0, self.yt0, self.mask, fixed_img=self.others_img,
                     box=1.0, ngrid=3, grid=0.3)
        return r

    def fit_diff(self, D, sig, x0, y0):
        return fit_free(D, sig, self.psf, x0, y0, self.mask, box=2.0, ngrid=5, grid=0.6)

    def fit_fixed(self, D, sig, x, y):
        m = self.mask
        w = 1 / sig[m]
        A = np.column_stack([self.psf(x, y)[m], np.ones(m.sum())]) * w[:, None]
        c, *_ = np.linalg.lstsq(A, D[m] * w, rcond=None)
        return c, float(((A @ c - D[m] * w) ** 2).sum())

    def analyse(self, D, sig, B, direct=None):
        rd = self.fit_direct(B, sig) if direct is None else direct
        xT, yT = rd['x'], rd['y']
        r = self.fit_diff(D, sig, xT, yT)
        off = self.jac @ np.array([r['x'] - xT, r['y'] - yT])
        cs = self.jac @ r['cov'] @ self.jac.T
        cT, chiT = self.fit_fixed(D, sig, xT, yT)
        return dict(xT=xT, yT=yT, FT=rd['amp'], bT=float(rd['coef'][-1]), xD=r['x'], yD=r['y'], amp=r['amp'],
                    dE=float(off[0]), dN=float(off[1]), cov_sky=cs, chi2_free=r['chi2'], chi2_target=chiT,
                    dchi2_target=chiT - r['chi2'], amp_target=float(cT[0]), npix=r['npix'], at_bound=r['at_bound'],
                    depth_peak=float(r['amp'] / rd['amp']))

    def empirical_target_image(self, B, sig, rd):
        """(B - const - tied others)/F_T : empirical target PSF (unit flux)."""
        oth = self.others_img(rd['x'], rd['y']) * rd['amp']
        E = (B - rd['coef'][-1] - oth) / rd['amp']
        return E
