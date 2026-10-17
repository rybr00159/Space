"""Quick-look of DSS plates (SERC-J 1975, SERC-I 1979, SERC-QV 1987, AAO-SES R 1992) around TIC 214311908 with the
star's position at the plate epoch, its 2019-2026 (TESS) positions and Gaia sources (propagated) overplotted.
Output: figs/dss_look.png"""
import glob, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/star'
g = pd.read_csv(f'{B}/data/gaia_120.csv')
M3V = 5925373272850784768
r = g[g.source_id == M3V].iloc[0]
fs = ['dss_poss2ukstu_blue.fits', 'dss_poss2ukstu_ir.fits', 'dss_phase2_gsc1.fits', 'dss_poss2ukstu_red.fits']
fig, axs = plt.subplots(1, 4, figsize=(16, 4.4))
def pos(ra, dec, pmra, pmdec, ep):
    dt = ep - 2016.0
    return ra + np.nan_to_num(pmra) * dt / 3.6e6 / np.cos(np.radians(dec)), dec + np.nan_to_num(pmdec) * dt / 3.6e6
for ax, f in zip(axs, fs):
    h = fits.open(f'{B}/data/img/dss/{f}')[0]; w = WCS(h.header)
    ep = Time(h.header['DATE-OBS'][:10]).jyear
    im = h.data.astype(float)
    x0, y0 = w.all_world2pix(*pos(r.ra, r.dec, r.pmra, r.pmdec, ep), 0)
    sub = 30 if im.shape[0] > 150 else 18
    xs, ys = int(round(float(x0))), int(round(float(y0)))
    cut = im[ys - sub:ys + sub + 1, xs - sub:xs + sub + 1]
    v1, v2 = np.percentile(cut, [5, 99.5])
    ax.imshow(cut, origin='lower', cmap='gray_r', vmin=v1, vmax=v2, extent=(xs - sub - 0.5, xs + sub + 0.5, ys - sub - 0.5, ys + sub + 0.5))
    gx, gy = w.all_world2pix(*pos(g.ra.values, g.dec.values, g.pmra.values, g.pmdec.values, ep), 0)
    k = (np.abs(gx - xs) < sub) & (np.abs(gy - ys) < sub)
    ax.scatter(gx[k], gy[k], s=60, facecolors='none', edgecolors='C0', lw=0.6)
    for k_, G in zip(np.where(k)[0], g.phot_g_mean_mag.values[k]):
        ax.text(gx[k_] + 0.8, gy[k_] + 0.8, f'{G:.1f}', color='C0', fontsize=6)
    for ept, col in [(2019.4, 'C1'), (2026.3, 'C3')]:
        xt, yt = w.all_world2pix(*pos(r.ra, r.dec, r.pmra, r.pmdec, ept), 0)
        ax.plot(xt, yt, '+', color=col, ms=12, mew=1.5)
    ax.plot(x0, y0, 'x', color='C2', ms=10)
    ax.set_title(f"{h.header.get('SURVEY')} {ep:.1f} ({h.header.get('EMULSION')}/{h.header.get('FILTER')})\nx: M3V then; +: M3V in 2019 (orange), 2026 (red)", fontsize=8)
fig.tight_layout(); fig.savefig(f'{B}/figs/dss_look.png', dpi=110)
