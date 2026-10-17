"""f6: localisation.png - which star dims?
(a) TESS S39 TESScut median image (zoom) with Gaia DR3 stars at the sector epoch (proper motion applied, referee_loc
    WCS offset), and the Mira's measured apparent image position (referee_loc/out/wcs_mira.json; lateral colour);
(b) raw-pixel dip S/N map for S39, recomputed with the model-free method of referee2_search r07 (common mode = median
    of pixels > 6 px from the Mira; per pixel and transit: local linear baseline at 0.06 < |dt| < 0.30 d minus the
    mean at |dt| < 0.022 d; errors from 80 fake 4-transit 'sectors' at random epochs >= 0.4 d from any transit),
    with the adopted ephemeris;
(c) combined six-sector localisation map of referee_loc (chi2 of a free-position source relative to its minimum, M-dwarf
    colour frame, arcsec from the M3V; inj_base.npz real_gsum) with the injection-calibrated 68/95/99.7% regions
    (locregion_base.json), and the catalogued stars;
(d) hypothesis tests of referee_loc (fit_base.json): chi2(star) - chi2(M3V) observed, and the 2.5-97.5% range
    expected if the M3V is the source and if that star is the source (1,000 injections each)."""
import sys, os, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Circle

style()
SEC = 39
RL = f'{B}/referee_loc'
STARS = [('M3V', 5925373272850784768, 'M3V (target)'),
         ('Mira', 5925373268528684288, 'Mira'),
         ('pairA', 5925373272850782208, 'pair A'),
         ('pairB', 5925373272827253632, 'pair B'),
         ('G175', 5925373272827260544, '3.1" star (G 17.5)'),
         ('G177', 5925373272827262208, '5.5" star (G 17.7)'),
         ('G170', 5925373272827256576, '12" star (G 17.0)')]
COL = {'M3V': INK, 'Mira': ORANGE, 'pairA': BLUE, 'pairB': BLUE, 'G175': VIOLET, 'G177': VIOLET, 'G170': VIOLET}
MK = {'M3V': '*', 'Mira': 'o', 'pairA': 's', 'pairB': 's', 'G175': 'D', 'G177': 'D', 'G170': 'D'}


def dip_map(t, F, tcs, rng):
    def epoch_depth(tc):
        k = np.abs(t - tc) < 0.30
        inn = np.abs(t - tc) < 0.022
        base = k & (np.abs(t - tc) > 0.06)
        if inn.sum() < 1 or base.sum() < 8: return None
        X = np.column_stack([np.ones(base.sum()), t[base] - tc])
        c, *_ = np.linalg.lstsq(X, F[base], rcond=None)
        return c[0] - F[inn].mean(0)
    D = [v for v in (epoch_depth(tc) for tc in tcs) if v is not None]
    Dm = np.mean(D, 0)
    ph = lambda tc: ((tc - T0_ADOPT) / P_ADOPT + 0.5) % 1 - 0.5
    fk = []
    while len(fk) < 80:
        vals = []
        for _ in range(len(D)):
            tc2 = rng.uniform(t.min() + 0.3, t.max() - 0.3)
            if abs(ph(tc2)) * P_ADOPT < 0.4: continue
            v = epoch_depth(tc2)
            if v is not None: vals.append(v)
        if len(vals) == len(D): fk.append(np.mean(vals, 0))
    fk = np.array(fk)
    E = 1.4826 * np.median(np.abs(fk - np.median(fk, 0)), 0)
    return Dm, E, len(D), fk


def main():
    rng = np.random.default_rng(707)
    fn = glob.glob(f'{CUBES}/tess-s{SEC:04d}-*_astrocut.fits')[0]
    with fits.open(fn) as h:
        d = h[1].data; w = WCS(h[2].header)
        t = np.array(d['TIME'], float); F = np.array(d['FLUX'], np.float64); q = np.array(d['QUALITY'], int)
    ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    t, F = t[ok], F[ok]
    ny, nx = F.shape[1:]
    med = np.median(F, 0)
    wm = json.load(open(f'{RL}/out/wcs_mira.json'))[str(SEC)]
    off = np.array(wm['offset_px'])
    g = pd.read_csv(f'{RL}/data/gaia_400.csv')
    ep = 2016.0 + (np.median(t) + 2457000 - 2457389.0) / 365.25
    dtt = ep - 2016.0
    ra = g.ra.values + np.nan_to_num(g.pmra.values) * dtt / 3.6e6 / np.cos(np.radians(g.dec.values))
    de = g.dec.values + np.nan_to_num(g.pmdec.values) * dtt / 3.6e6
    px, py = w.all_world2pix(ra, de, 0)
    g['px'], g['py'] = px + off[0], py + off[1]
    pos = {k: g[g.source_id == sid][['px', 'py']].values[0] for k, sid, _ in STARS}
    # local Jacobian (arcsec E,N per pixel) for the Mira's apparent offset
    r0, d0 = w.all_pix2world([[pos['M3V'][0] - off[0], pos['M3V'][1] - off[1]]], 0)[0]
    def sky(x, y):
        rr, dd = w.all_pix2world([[x, y]], 0)[0]
        return np.array([(rr - r0) * np.cos(np.radians(d0)) * 3600, (dd - d0) * 3600])
    x0, y0 = pos['M3V'][0] - off[0], pos['M3V'][1] - off[1]
    J = np.column_stack([(sky(x0 + 0.01, y0) - sky(x0, y0)) / 0.01, (sky(x0, y0 + 0.01) - sky(x0, y0)) / 0.01])
    mira_dEN = np.mean([v['dEN'] for v in wm['mira'].values()], 0)
    mira_app = pos['Mira'] + np.linalg.solve(J, mira_dEN)
    # time correction is irrelevant at the 0.3-d window scale (target barycentre differs by < 31 s from TESScut TIME)
    tcs = [T0_ADOPT + e * P_ADOPT for e in range(-450, 100) if t.min() + 0.3 < T0_ADOPT + e * P_ADOPT < t.max() - 0.3]
    yy, xx = np.mgrid[0:ny, 0:nx]
    cm = np.hypot(xx - pos['Mira'][0], yy - pos['Mira'][1]).ravel() > 6
    Ff = F.reshape(len(t), -1)
    Ff = Ff - np.median(Ff[:, cm], axis=1)[:, None]
    Dm, E, ntr, fk = dip_map(t, Ff, tcs, rng)
    SN = (Dm / E).reshape(ny, nx)
    box = ((np.abs(xx - pos['M3V'][0]) <= 1.5) & (np.abs(yy - pos['M3V'][1]) <= 1.5)).ravel()
    box_sum = float(Dm[box].sum()); box_err = float(np.std(fk[:, box].sum(1)))
    far = np.hypot(xx - pos['M3V'][0], yy - pos['M3V'][1]).ravel() > 2.5
    print(f'S{SEC}: {ntr} transits; 3x3-box dip {box_sum:.1f} +- {box_err:.1f} e/s; max far-pixel S/N {SN.ravel()[far].max():.2f}; '
          f'max |S/N| overall {np.abs(SN).max():.2f}')
    r07 = json.load(open(f'{B}/referee2_search/out/r07_pixdepth.json'))
    # ---------------- figure
    fig = plt.figure(figsize=(10.6, 9.8))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.34, height_ratios=[1, 0.92])
    ax = fig.add_subplot(gs[0, 0]); bx = fig.add_subplot(gs[0, 1]); cx = fig.add_subplot(gs[1, 0]); dx = fig.add_subplot(gs[1, 1])
    # (a) image
    cx0, cy0 = pos['M3V']
    hw = 2.3
    im = ax.imshow(np.clip(med, 50, None), origin='lower', cmap='Greys', norm=LogNorm(vmin=150, vmax=med.max()), interpolation='nearest')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.set_label('median flux (e$^-$/s per pixel)', fontsize=8)
    ax.set_xlim(cx0 - hw, cx0 + hw); ax.set_ylim(cy0 - hw, cy0 + hw)
    named = [sid for _, sid, _ in STARS]
    gg = g[(np.abs(g.px - cx0) < hw) & (np.abs(g.py - cy0) < hw) & (g.phot_g_mean_mag < 18.5) & (~g.source_id.isin(named))]
    ax.plot(gg.px, gg.py, 'o', color=MUTED, mfc='none', ms=4, mew=0.8, alpha=0.9)
    lab_off = {'M3V': (-1.2, -1.25), 'Mira': (1.0, 1.25), 'pairA': (-0.3, -1.75), 'pairB': (1.3, -0.75), 'G175': (-1.35, 0.55),
               'G177': (-0.9, 1.45), 'G170': (0.6, 1.75)}
    for k, sid, name in STARS:
        p = pos[k]
        ax.plot(p[0], p[1], MK[k], color=COL[k], ms=12 if k == 'M3V' else 6.5, mec='white', mew=0.8, zorder=5)
        lo = lab_off[k]
        ax.annotate(name, xy=p, xytext=(p[0] + lo[0], p[1] + lo[1]), fontsize=7.6, color=COL[k] if k != 'M3V' else INK,
                    fontweight='bold' if k == 'M3V' else None, ha='center',
                    arrowprops=dict(arrowstyle='-', color=COL[k], lw=0.7), zorder=6,
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    ax.plot(mira_app[0], mira_app[1], 'o', mfc='none', mec=ORANGE, ms=9, mew=1.4, zorder=5)
    ax.annotate("Mira's apparent image\n(lateral colour)", xy=mira_app, xytext=(mira_app[0] + 0.3, mira_app[1] - 1.5), fontsize=7.2, color=ORANGE,
                ha='center', arrowprops=dict(arrowstyle='-', color=ORANGE, lw=0.7), bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    ax.set_xlabel('TESScut column (pixel, 21")'); ax.set_ylabel('TESScut row (pixel)')
    ax.set_title(f'$\\bf{{a}}$   TESS S{SEC} image, Gaia DR3 stars', loc='left')
    ax.grid(False)
    # (b) dip S/N
    lim = 6
    im2 = bx.imshow(SN, origin='lower', cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), interpolation='nearest')
    cb2 = fig.colorbar(im2, ax=bx, fraction=0.046, pad=0.02); cb2.set_label('dip S/N per pixel (+ = dimming in transit)', fontsize=8)
    bx.plot(*pos['M3V'], '*', color=INK, ms=11, mec='white', mew=0.8)
    bx.plot(*pos['Mira'], 'o', color=ORANGE, ms=6, mec='white', mew=0.7)
    bx.plot(*pos['pairA'], 's', color=BLUE, ms=5, mec='white', mew=0.7)
    bx.add_patch(Circle(pos['M3V'], 2.5, fill=False, ec=INK2, lw=1.0, ls='--'))
    bx.text(0.02, 0.02, f'S{SEC}, {ntr} transits, raw pixels, no scene model\n3x3 box at the M3V (also holds the Mira core): {box_sum:.1f} $\\pm$ {box_err:.1f} e$^-$/s'
            f'\nall six sectors: 13.9 $\\pm$ 1.4 e$^-$/s (referee2_search r07)\nno pixel beyond 2.5 px (dashed) dims\nin any sector (r07: p 0.11-0.90)',
            transform=bx.transAxes, fontsize=7.3, color=INK, va='bottom', bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85))
    bx.set_xlabel('TESScut column (pixel)'); bx.set_ylabel('TESScut row (pixel)')
    bx.set_title('$\\bf{b}$   Where the pixels dim (whole 21x21 cutout)', loc='left')
    bx.grid(False)
    # (c) localisation map
    Z = np.load(f'{RL}/out/inj_base.npz')
    GE, GN = Z['GE'], Z['GN']
    m = Z['real_gsum'].reshape(len(GN), len(GE)); m = m - m.min()
    LR = json.load(open(f'{RL}/out/locregion_base.json'))
    thr = [LR['68_injection']['threshold'], LR['95_injection']['threshold'], LR['99.7_injection']['threshold']]
    im3 = cx.imshow(np.sqrt(np.clip(m, 0, 400)), origin='lower', extent=[GE[0] - 0.25, GE[-1] + 0.25, GN[0] - 0.25, GN[-1] + 0.25],
                    cmap='cividis_r', aspect='equal', interpolation='nearest')
    cb3 = fig.colorbar(im3, ax=cx, orientation='horizontal', fraction=0.05, pad=0.2, shrink=0.8)
    cb3.set_label('$\\sqrt{\\Delta\\chi^2}$ of a source at that position', fontsize=8)
    cs = cx.contour(GE, GN, m, levels=thr, colors=['white', 'white', 'white'], linewidths=[1.6, 1.1, 0.8], linestyles=['-', '--', ':'])
    off_ = LR['star_offsets_from_M3V_2025']
    names = {'M3V': 'M3V', 'pairA': 'pair A', 'pairB': 'pair B', 'Mira': 'Mira', 'G175': '3.1"', 'G177': '5.5"', 'G170': '12"'}
    for k, (e_, n_) in off_.items():
        cx.plot(e_, n_, MK[k], color=COL[k] if k != 'M3V' else 'white', ms=11 if k == 'M3V' else 6, mec=INK if k == 'M3V' else 'white', mew=0.8, zorder=5)
        dy = {'M3V': -2.3, 'pairA': -2.6, 'pairB': -2.6, 'Mira': 1.4, 'G175': 1.6, 'G177': 1.6, 'G170': 1.5}[k]
        dxl = {'M3V': 0, 'pairA': -1.5, 'pairB': 2.5, 'Mira': 0, 'G175': -1.4, 'G177': 0.4, 'G170': 0}[k]
        cx.text(e_ + dxl, n_ + dy, names[k], fontsize=7.4, color='white' if k in ('M3V',) else INK, ha='center', va='center',
                fontweight='bold', bbox=dict(boxstyle='round,pad=0.12', fc=('white' if k != 'M3V' else INK), ec='none', alpha=0.8))
    cx.set_xlim(22, -12); cx.set_ylim(-9, 9)
    cx.set_xlabel('offset East from the M3V (arcsec; East to the left)'); cx.set_ylabel('offset North (arcsec)')
    cx.set_title('$\\bf{c}$   Six-sector localisation (sky)', loc='left')
    from matplotlib.lines import Line2D
    cx.legend([Line2D([], [], color=INK2, lw=1.6), Line2D([], [], color=INK2, lw=1.1, ls='--'), Line2D([], [], color=INK2, lw=0.8, ls=':')],
              ['68%', '95%', '99.7%'], loc='lower right', fontsize=7.2, ncol=3, title='injection-calibrated regions', title_fontsize=7.2, handlelength=1.6, columnspacing=1.0,
              facecolor='white', framealpha=0.85, frameon=True, edgecolor='none')
    cx.grid(False)
    # (d) hypothesis tests
    FB = json.load(open(f'{RL}/out/fit_base.json'))
    main_ = FB['main']; real = FB['real']['comb_dchi2']
    hyps = ['pairA', 'pairB', 'Mira', 'G175', 'G177', 'G170']
    hname = {'pairA': 'pair A (11.7")', 'pairB': 'pair B (12.3")', 'Mira': 'Mira (18.8")', 'G175': '3.1" star', 'G177': '5.5" star', 'G170': '12" star'}
    yv = np.arange(len(hyps))[::-1]
    for i, k in enumerate(hyps):
        j = main_.index(k)
        qm = [FB['inj']['M3V']['q'][n][j] for n in (1, 5)]       # 2.5, 97.5% if the M3V is the source
        qt = [FB['inj'][k]['q'][n][j] for n in (1, 5)]           # if that star is the source
        dx.plot(qm, [yv[i] + 0.15] * 2, '-', color=INK2, lw=5, alpha=0.35, solid_capstyle='butt', label='expected if the M3V is the source' if i == 0 else None)
        dx.plot(qt, [yv[i] - 0.15] * 2, '-', color=ORANGE, lw=5, alpha=0.55, solid_capstyle='butt', label='expected if that star is the source' if i == 0 else None)
        dx.plot(real[k], yv[i], 'D', color=INK, ms=6.5, mec='white', mew=0.7, zorder=4, label='observed' if i == 0 else None)
        dx.text(real[k] + 12, yv[i] + 0.02, f'+{real[k]:.0f}', fontsize=7.6, va='center', color=INK)
    dx.axvline(0, color=MUTED, lw=0.9)
    dx.set_yticks(yv); dx.set_yticklabels([hname[k] for k in hyps], fontsize=8.2)
    dx.set_xlabel('$\\chi^2$(source on that star) $-$ $\\chi^2$(source on the M3V)\n(positive: the M3V fits better)')
    dx.set_xlim(-260, 460)
    dx.legend(loc='upper center', bbox_to_anchor=(0.45, -0.2), fontsize=7.4, handlelength=1.6, ncol=2)
    dx.set_title('$\\bf{d}$   Pixel-level hypothesis tests (2.5-97.5%)', loc='left')
    dx.grid(axis='y', visible=False)
    save(fig, 'localisation.png')
    print('thresholds', thr, 'M3V map value', round(float(m[np.argmin(np.abs(GN)), np.argmin(np.abs(GE))]), 2))
    print('Mira apparent offset (E,N) arcsec', mira_dEN.round(2), 'r07 S39 sum', r07[str(SEC)]['sum_depth_3x3_e_s'], r07[str(SEC)]['sum_depth_3x3_err'])


if __name__ == '__main__':
    main()
