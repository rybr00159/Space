"""Per-sector model: image products for real and fake windows, noise map, hypothesis PRF matrices and a sky grid,
closed-form weighted fits of  d = a * PRF(position) + b  (a = total flux lost by the source, e-/s; b = constant)."""
import numpy as np, json, os
from l0lib import *
from l1img import Data, window_sets

HYP_RADIUS = 80.0      # arcsec from the Mira
HYP_TMAX = 17.6        # a 15 e/s dip needs F >= 15 e/s (T <= 17.5) even for a 100% eclipse
R_FIT = 5.5            # px, fit region radius around the Mira's apparent position
FAKE_SEP = 0.7         # d, minimum separation of fake windows within one draw


class SectorModel:
    def __init__(self, sec, method='DIFF', grid_colour=2.76, grid_half=100.0, grid_step=1.0, S=None, Dd=None,
                 psf_override=None, hyp_extra=None, mira_mode='measured', w_oot=None, r_fit=None):
        self.sec = sec; self.method = method
        self.S = S if S is not None else Sector(sec)
        self.Dd = Dd if Dd is not None else Data(self.S)
        S, Dd = self.S, self.Dd
        self.psf = psf_override if psf_override is not None else S.prf()
        tr, fk = window_sets(Dd, step=0.05)
        if w_oot is None:
            prod = (lambda tc: Dd.diff_window(tc)) if method == 'DIFF' else (lambda tc: Dd.amp_window(tc))
        else:
            prod = (lambda tc: Dd.diff_window(tc, w_oot=w_oot)) if method == 'DIFF' else (lambda tc: Dd.amp_window(tc))
        self.real = [(tc, p) for tc, p in ((tc, prod(tc)) for tc in tr) if p is not None]
        fake = [(tc, p) for tc, p in ((tc, prod(tc)) for tc in fk) if p is not None]
        self.fake_tc = np.array([tc for tc, _ in fake]); self.fake = [p for _, p in fake]
        # Mira apparent position (measured, c02_mira.py)
        mp = json.load(open(f'{BASE}/out/mira_pos.json'))[str(sec)]['base']['dEN']
        self.mira_dEN = np.array(mp)
        self.mira_xy = S.sky_to_xy(self.mira_dEN)[0]
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]
        self.mask = (np.hypot(xx - self.mira_xy[0], yy - self.mira_xy[1]) <= (R_FIT if r_fit is None else r_fit)).ravel()
        self.npix = int(self.mask.sum())
        # per-window vectors on the region
        if method == 'DIFF':
            self.rv = [(p['D'][self.mask], p['V'][self.mask]) for _, p in self.real]
            self.fv = [(p['D'][self.mask], p['V'][self.mask]) for p in self.fake]
        else:
            sa = Dd.sig_amp[self.mask]
            self.rv = [(-p['num'][self.mask], p['den'], sa) for _, p in self.real]
            self.fv = [(-p['num'][self.mask], p['den'], sa) for p in self.fake]
        # empirical inflation of the formal errors, from single fake windows
        z = np.array([self._single(v)[0] / self._single(v)[1] for v in self.fv])
        rr = z.std(0)
        self.infl_pix = np.sqrt(0.5 * rr ** 2 + 0.5 * np.median(rr) ** 2)
        self.infl_med = float(np.median(rr))
        d, s = self.combine_vecs(self.rv)
        self.d_real, self.sig = d, s * self.infl_pix
        self.w = 1 / self.sig ** 2
        # hypotheses: catalogued stars + measured Mira position
        g = stars()
        xy = S.star_xy(g)
        sep = np.hypot(g.dE.values, g.dN.values)   # Gaia 2016 offsets from the Mira (arcsec)
        sel = (sep < HYP_RADIUS) & (g['T'].values < HYP_TMAX)
        self.hyp = g[sel].copy().reset_index(drop=True)
        hxy = xy[sel]
        # replace the Mira's predicted position by the measured apparent position
        im = np.where(self.hyp.source_id.values == IDS['Mira'])[0][0]
        if mira_mode == 'measured':
            hxy[im] = self.mira_xy                      # measured apparent position (c02_mira.py)
        elif mira_mode == 'gaia':
            hxy[im] = S.star_xy(g[g.source_id == IDS['Mira']], colour=False)[0]
        # 'colourcurve': keep Gaia + extrapolated colour-curve shift (BP-RP 4.68)
        self.hyp['x'] = hxy[:, 0]; self.hyp['y'] = hxy[:, 1]
        self.hyp['name'] = [next((k for k, v in IDS.items() if v == s), str(s)) for s in self.hyp.source_id]
        if hyp_extra is not None:       # list of (name, x, y)
            for nm, x, y in hyp_extra:
                self.hyp.loc[len(self.hyp)] = {**{c: np.nan for c in self.hyp.columns}, 'name': nm, 'x': x, 'y': y, 'source_id': -1}
        self.H = np.array([self.psf(x, y).ravel()[self.mask] for x, y in zip(self.hyp.x, self.hyp.y)])
        # sky grid (dE, dN) relative to the M3V's position AT THIS EPOCH (co-moving with the M3V), for a source of
        # colour grid_colour.  self.m3v_off = M3V offset from the Mira's Gaia position at this epoch.
        r3 = g[g.source_id == IDS['M3V']].iloc[0]; rm = g[g.source_id == IDS['Mira']].iloc[0]
        dt = S.epoch - 2016.0
        self.m3v_off = np.array([(r3.ra - rm.ra) * np.cos(np.radians(rm.dec)) * 3600 + (r3.pmra - rm.pmra) * dt / 1000,
                                 (r3.dec - rm.dec) * 3600 + (r3.pmdec - rm.pmdec) * dt / 1000])
        ax = np.arange(-grid_half, grid_half + 1e-9, grid_step)
        GE, GN = np.meshgrid(ax, ax)
        self.grid = np.column_stack([GE.ravel(), GN.ravel()])
        self.grid_axis = ax
        gxy = S.sky_to_xy(self.grid + self.m3v_off[None, :], C=np.full(len(self.grid), grid_colour) if grid_colour is not None else None)
        self.grid_xy = gxy
        self.G = np.array([self.psf(x, y).ravel()[self.mask] for x, y in gxy])
        self._prep()

    def _single(self, v):
        if self.method == 'DIFF':
            return v[0], np.sqrt(v[1])
        return v[0] / v[1], v[2] / np.sqrt(v[1])

    def combine_vecs(self, vs):
        if self.method == 'DIFF':
            w = np.array([1 / v[1] for v in vs]); D = np.array([v[0] for v in vs])
            return (D * w).sum(0) / w.sum(0), 1 / np.sqrt(w.sum(0))
        num = sum(v[0] for v in vs); den = sum(v[1] for v in vs)
        return num / den, vs[0][2] / np.sqrt(den)

    def _prep(self):
        w = self.w
        def pre(M):
            return dict(Sw=w.sum(), Sp=(M * w).sum(1), Spp=(M * M * w).sum(1), Mw=M * w)
        self.pH = pre(self.H); self.pG = pre(self.G)

    def fit(self, d, which='H'):
        """chi2 and amplitude for each column set (hypotheses 'H' or grid 'G') given data vector d."""
        p = self.pH if which == 'H' else self.pG
        w = self.w
        Sd = (w * d).sum(); Sdd = (w * d * d).sum()
        Spd = p['Mw'] @ d
        det = p['Spp'] * p['Sw'] - p['Sp'] ** 2
        a = (Spd * p['Sw'] - p['Sp'] * Sd) / det
        b = (p['Spp'] * Sd - p['Sp'] * Spd) / det
        chi2 = Sdd - a * Spd - b * Sd
        sa = np.sqrt(p['Sw'] / det)
        return chi2, a, b, sa

    def null_chi2(self, d):
        w = self.w; b = (w * d).sum() / w.sum()
        return float((w * (d - b) ** 2).sum())

    def fake_draw(self, rng, n=None):
        n = len(self.rv) if n is None else n
        for _ in range(200):
            idx = rng.permutation(len(self.fv))
            sel = []
            for i in idx:
                if all(abs(self.fake_tc[i] - self.fake_tc[j]) >= FAKE_SEP for j in sel):
                    sel.append(i)
                if len(sel) == n: break
            if len(sel) == n: return sel
        raise RuntimeError('cannot draw fakes')

    def fake_combo(self, sel):
        return self.combine_vecs([self.fv[i] for i in sel])[0]

    def grid_refine(self, chi2g):
        """Sub-grid minimum by a quadratic fit to the 3x3 neighbourhood of the grid minimum."""
        n = len(self.grid_axis); i = int(np.argmin(chi2g)); iy, ix = divmod(i, n)
        ix = min(max(ix, 1), n - 2); iy = min(max(iy, 1), n - 2)
        Z = chi2g.reshape(n, n)[iy - 1:iy + 2, ix - 1:ix + 2]
        h = self.grid_axis[1] - self.grid_axis[0]
        gx = (Z[1, 2] - Z[1, 0]) / 2; gy = (Z[2, 1] - Z[0, 1]) / 2
        hxx = Z[1, 2] - 2 * Z[1, 1] + Z[1, 0]; hyy = Z[2, 1] - 2 * Z[1, 1] + Z[0, 1]
        hxy = (Z[2, 2] - Z[2, 0] - Z[0, 2] + Z[0, 0]) / 4
        Hm = np.array([[hxx, hxy], [hxy, hyy]])
        try:
            dlt = -np.linalg.solve(Hm, [gx, gy]) if np.all(np.linalg.eigvalsh(Hm) > 0) else np.zeros(2)
        except np.linalg.LinAlgError:
            dlt = np.zeros(2)
        dlt = np.clip(dlt, -1, 1)
        return np.array([self.grid_axis[ix] + dlt[0] * h, self.grid_axis[iy] + dlt[1] * h]), float(Z.min())
