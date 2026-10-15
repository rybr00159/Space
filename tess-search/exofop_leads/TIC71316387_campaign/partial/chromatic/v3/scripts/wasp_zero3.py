"""SuperWASP flux zero-level check (additive flux bias -> depth scale) from the tile330042 catalogue.

For WASP sources within 1.5 deg of the target, predict the 48"-aperture flux from Gaia DR3 (all sources
G<17.5 within 2 arcmin, WASP-band colour term from the warped BT-Settl spectra, calibrated PSF/aperture
of dilution3: R = 3.5 px, sigma = 25.4").  Fit  TAMMAG_median = ZP - 2.5 log10(F_pred + b)  with an
additive flux offset b (in units of the predicted flux of a G=12 star), plus a colour term, robustly.
An additive offset b (sky over/under-subtraction) rescales eclipse depths by F/(F+b).
"""
import sys, json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import dilution3 as DL3

D = L.V3 + '/data/wasp_tile/'
w = pd.read_csv(D + 'tile330042.csv')
g = pd.read_csv(D + 'gaia_1p6deg_G17p5.csv')
dc = json.load(open(L.V3 + '/results/dilution3.json'))
R, SIG = 3.5 * 13.7, dc['wasp_cal']['3.5']['sigma_arcsec']
ep = 2007.0 - 2016.0
g['ra7'] = g.ra + g.pmra.fillna(0) / 3.6e6 / np.cos(np.radians(g.dec)) * ep
g['de7'] = g.dec + g.pmdec.fillna(0) / 3.6e6 * ep
c = g.bp_rp.values.copy(); c[~np.isfinite(c)] = np.nanmedian(c)
ct = np.array([DL3.colour_term(x, 'WASP') for x in c])
g['Fw'] = 10 ** (-0.4 * (g.phot_g_mean_mag.values + ct - 12.0))       # WASP-band flux in units of a G=12 (ct=0) star
cosd = np.cos(np.radians(L.DEC))
sep0 = np.hypot((w.ra - L.RA) * cosd, w.dec - L.DEC) * 3600
w = w[(sep0 < 1.45 * 3600) & np.isfinite(w.tm_median)].copy()
GX, GY = g.ra7.values * cosd, g.de7.values
rows = []
for i, r in w.iterrows():
    dx = (GX - r.ra * cosd) * 3600; dy = (GY - r.dec) * 3600
    s = np.hypot(dx, dy)
    near = s < R + 5 * SIG
    if not near.any():
        continue
    ph = np.array([DL3.phi(x, R, SIG) for x in s[near]])
    contrib = g.Fw.values[near] * ph
    j = np.argmin(s[near])
    if s[near][j] > 6:             # no Gaia counterpart within 6"
        continue
    rows.append((r.sourceid, r.ra, r.dec, r.tm_median, r.wasp_mag, contrib.sum(), contrib[j] / contrib.sum(), c[near][j], s[near][j], r.npts))
T = pd.DataFrame(rows, columns=['id', 'ra', 'dec', 'tm', 'wmag', 'Fpred', 'purity', 'bprp', 'sep', 'npts'])
print('WASP sources with Gaia counterparts: %d' % len(T))
sel = (T.purity > 0.9) & (T.tm > 9.0) & (T.tm < 14.8) & (T.npts > 5000)
S = T[sel]
print('isolated (purity>0.9), 9<mag<14.8: %d' % len(S))


def model(p, F, col):
    zp, b, k = p
    return zp - 2.5 * np.log10(np.maximum(F + b, 1e-6)) + k * (col - 0.8)


def fit(S, with_b=True):
    def res(p):
        if not with_b:
            p = [p[0], 0.0, p[1]]
        return S.tm.values - model(p, S.Fpred.values, S.bprp.values)
    p0 = [12.0, 0.0, 0.0] if with_b else [12.0, 0.0]
    r = least_squares(res, p0, loss='soft_l1', f_scale=0.05)
    for _ in range(3):
        rr = res(r.x); s = 1.4826 * np.median(np.abs(rr))
        m = np.abs(rr) < 4 * s
        def res2(p, m=m):
            return res(p)[m]
        r = least_squares(res2, r.x, loss='soft_l1', f_scale=0.05)
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * (1.4826 * np.median(np.abs(r.fun))) ** 2
    return r.x, np.sqrt(np.diag(cov)), 1.4826 * np.median(np.abs(r.fun)), m.sum()


p, pe, sc, n = fit(S)
print('fit with additive offset: ZP %.3f  b = %.4f +- %.4f (flux units of a G=12 star)  colour term %.3f; robust scatter %.3f mag (n=%d)' % (p[0], p[1], pe[1], p[2], sc, n))
p0, pe0, sc0, n0 = fit(S, False)
print('fit without offset: ZP %.3f colour %.3f scatter %.3f' % (p0[0], p0[1], sc0))
# target
tg = T.loc[T.id == '1SWASP J220333.07-425007.9'].iloc[0]
pred = model(p, tg.Fpred, tg.bprp)
print('target: TAMMAG median %.3f, predicted %.3f (residual %+.3f mag); Fpred %.3f; purity %.3f' % (tg.tm, pred, tg.tm - pred, tg.Fpred, tg.purity))
# depth scale factor: depth_true = depth_meas * (F + b)/F, with F = the target's own aperture flux (excl. neighbours)
Ft = tg.Fpred * tg.purity
fac = (Ft + p[1] + (tg.Fpred - Ft)) / Ft    # (target + neighbours + offset) / target ; the neighbour part is the dilution D
fac_off = (tg.Fpred + p[1]) / tg.Fpred
print('depth scale from the additive offset alone: x%.3f +- %.3f   (dilution part from neighbours separately: D = %.3f)' % (
    fac_off, pe[1] / tg.Fpred, tg.Fpred / Ft - 1))
# binned residuals vs magnitude
S2 = S.copy(); S2['res'] = S2.tm - model([p0[0], 0, p0[1]], S2.Fpred, S2.bprp)
bins = np.arange(9, 15.01, 0.5)
for a, b_ in zip(bins[:-1], bins[1:]):
    m = (S2.tm >= a) & (S2.tm < b_)
    if m.sum() > 5:
        print('   mag %.1f-%.1f: n=%3d median residual (no-offset fit) %+.3f' % (a, b_, m.sum(), np.median(S2.res[m])))
json.dump(dict(b=float(p[1]), b_err=float(pe[1]), zp=float(p[0]), colour=float(p[2]), scatter=float(sc), n=int(n),
               target_resid=float(tg.tm - pred), depth_scale_offset=float(fac_off), depth_scale_offset_err=float(pe[1] / tg.Fpred),
               target_Fpred=float(tg.Fpred), target_purity=float(tg.purity)),
          open(L.V3 + '/results/wasp_zero3.json', 'w'), indent=1)
