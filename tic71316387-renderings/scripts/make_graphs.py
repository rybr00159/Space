"""Simple, labelled graphs of the TIC 71316387 (CTOI 71316387.01) findings.

Writes six PNGs to ../graphs/. The numbers printed on the graphs are the adopted values of the campaign report,
tess-search/exofop_leads/TIC71316387_campaign/README.md (2026-10-15), and its partial/*.json results.
The TESS light curves are downloaded from MAST into a temporary cache; nothing but the PNGs is written to the repo.

Requires numpy, scipy, matplotlib and astropy.  Usage:  python make_graphs.py
"""
import os
import tempfile
import urllib.request

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle, Circle
from matplotlib.ticker import FuncFormatter, MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(HERE, '..', 'graphs')
CACHE = os.environ.get('TIC713_CACHE', os.path.join(tempfile.gettempdir(), 'tic71316387_lc'))
CAMPAIGN = os.path.join(REPO, 'tess-search', 'exofop_leads', 'TIC71316387_campaign')

# ---- adopted values (campaign report / synthesis_result.json) --------------------------------------------
P = 1.30367228                      # eclipsing-binary period (d)
T0 = 3165.670416                    # primary (deep) eclipse, BTJD = BJD_TDB - 2457000
P_CTOI = P / 2                      # catalogued planet-candidate period
DEEP, SHALLOW = 5770, 2600          # diluted TESS depths (ppm)
M_G = 1.03                          # G-star mass (Msun)
# empirical eclipse shapes (1 - |x|^p)^q, x = dt / (T/2)  (chromatic_result.json, template shape)
SHAPE = {'deep': (0.0909, 3.38, 10.0), 'shallow': (0.0851, 4.36, 10.0)}

# ---- style (validated categorical slots 1-2; ink and chrome tokens) --------------------------------------
SURFACE, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
GRID, AXIS = '#e1e0d9', '#c3c2b7'
BLUE, ORANGE = '#2a78d6', '#eb6834'
BLUE_STEPS = ['#cde2fb', '#9ec5f4', '#5598e7']          # sequential 100 / 200 / 350
matplotlib.rcParams.update({
    'font.family': ['Inter', 'DejaVu Sans'], 'font.size': 10, 'text.color': INK,
    'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': AXIS, 'axes.linewidth': 0.8, 'axes.labelcolor': INK2, 'axes.labelsize': 10,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
    'xtick.color': AXIS, 'ytick.color': AXIS, 'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.frameon': False, 'legend.fontsize': 9,
})
PPM = FuncFormatter(lambda v, _: f'{v:,.0f}'.replace('-', '−'))


def frame(title, subtitle, foot, w=9.0, h=5.4):
    fig = plt.figure(figsize=(w, h))
    fig.text(0.055, 0.945, title, fontsize=15, weight='semibold', color=INK, va='top')
    fig.text(0.055, 0.885, subtitle, fontsize=10, color=INK2, va='top', wrap=True)
    fig.text(0.055, 0.025, foot, fontsize=7.5, color=MUTED, va='bottom', wrap=True)
    return fig


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name), dpi=200)
    plt.close(fig)
    print('wrote', name)


# ---- TESS light curves ----------------------------------------------------------------------------------
def load(sector):
    """TESS-SPOC full-frame-image light curve of TIC 71316387 (PDCSAP flux), normalised and detrended."""
    os.makedirs(CACHE, exist_ok=True)
    name = f'hlsp_tess-spoc_tess_phot_0000000071316387-s{sector:04d}_tess_v1_lc.fits'
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        uri = f'mast:HLSP/tess-spoc/s{sector:04d}/target/0000/0000/7131/6387/{name}'
        urllib.request.urlretrieve('https://mast.stsci.edu/api/v0.1/Download/file?uri=' + uri, path)
    with fits.open(path) as h:
        d = h[1].data
        t = np.asarray(d['TIME'], float)
        f = np.asarray(d['PDCSAP_FLUX'], float)
        q = np.asarray(d['QUALITY'])
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    t, f = t[ok], f[ok] / np.median(f[ok])
    return dict(sector=sector, t=t, f=detrend(t, f), texp=np.median(np.diff(t)))


def dt_to(t, which):
    """Time from the nearest deep (which='deep') or shallow eclipse, and the cycle number."""
    tc = T0 + (0.0 if which == 'deep' else 0.5) * P
    n = np.round((t - tc) / P)
    return t - (tc + n * P), n


def detrend(t, f, w=0.75):
    """Divide by a running median of the out-of-eclipse flux (removes spot modulation and drifts)."""
    m = (np.abs(dt_to(t, 'deep')[0]) < 0.07) | (np.abs(dt_to(t, 'shallow')[0]) < 0.07)
    tt, ff = t[~m], f[~m]
    grid = np.arange(t.min(), t.max() + 0.02, 0.02)
    lo, hi = np.searchsorted(tt, grid - w / 2), np.searchsorted(tt, grid + w / 2)
    med = np.array([np.median(ff[a:b]) if b - a > 10 else np.nan for a, b in zip(lo, hi)])
    good = np.isfinite(med)
    return f / np.interp(t, grid[good], med[good])


def shape(dt, which, texp):
    """Eclipse template (1 at mid-eclipse) averaged over the exposure time."""
    T, p, q = SHAPE[which]
    ns = max(1, int(round(texp * 1440 / 2)))
    off = ((np.arange(ns) + 0.5) / ns - 0.5) * texp
    x = np.abs(dt[:, None] + off[None, :]) / (T / 2)
    return np.where(x < 1, (1 - np.clip(x, 0, 1) ** p) ** q, 0.0).mean(axis=1)


# ---- graph 1: the raw alternation ------------------------------------------------------------------------
def graph_alternation(lc):
    t, f = lc['t'], lc['f']
    # a 4.6-day stretch with no gaps longer than 30 min, starting just before a deep eclipse
    span = 4.6
    best = None
    for n in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P)):
        a = T0 + n * P - 0.3
        s = (t >= a) & (t <= a + span)
        if s.sum() < 0.9 * span / lc['texp']:
            continue
        if np.max(np.diff(t[s])) < 30 / 1440:
            best = a
            break
    s = (t >= best) & (t <= best + span)
    tb = np.arange(best, best + span + 1e-9, 15 / 1440)
    idx = np.digitize(t[s], tb)
    xs = np.array([t[s][idx == i].mean() for i in np.unique(idx)]) - best
    ys = np.array([np.median(f[s][idx == i]) for i in np.unique(idx)]) * 1e6 - 1e6
    tm = np.linspace(best, best + span, 4000)
    model = -(DEEP * shape(dt_to(tm, 'deep')[0], 'deep', 0) + SHALLOW * shape(dt_to(tm, 'shallow')[0], 'shallow', 0))

    fig = frame('Every other dip is deeper',
                'Four and a half days of TESS data. The catalogued "planet" dips every 0.65 days, '
                'but the dips alternate between two depths.',
                'TESS Sector 68 (July–August 2023), TESS-SPOC light curve in 15-minute medians, slow brightness '
                'variations removed. Line and depths: adopted eclipse model of the campaign report.')
    ax = fig.add_axes([0.085, 0.17, 0.885, 0.585])
    ax.plot(tm - best, model, color=INK, lw=1.1, alpha=0.85, zorder=2, label='eclipse model')
    ax.plot(xs, ys, 'o', color=MUTED, ms=3.0, mew=0, zorder=3, label='TESS, 15-min medians')
    ax.yaxis.set_major_formatter(PPM)
    ax.set_xlim(0, span)
    ax.set_ylim(-8400, 4300)
    ax.set_yticks([-6000, -4000, -2000, 0, 2000])
    ax.set_xlabel('Days')
    ax.set_ylabel('Brightness change (ppm)')
    # a row of markers under the dips: which eclipse is which
    tcs = []
    for which, col in (('deep', BLUE), ('shallow', ORANGE)):
        n0 = np.ceil((best - (T0 + (which == 'shallow') * 0.5 * P)) / P)
        for n in np.arange(n0, n0 + 5):
            tc = T0 + (n + (which == 'shallow') * 0.5) * P - best
            if 0 < tc < span:
                ax.plot(tc, -7600, marker='^', ms=10, color=col, mec=SURFACE, mew=1.5, zorder=4)
                tcs.append((tc, which))
    tcs.sort()
    ax.plot([], [], marker='^', ms=10, color=BLUE, mec=SURFACE, mew=1.5, ls='none',
            label=f'deeper eclipse, {DEEP:,} ppm (0.58%)')
    ax.plot([], [], marker='^', ms=10, color=ORANGE, mec=SURFACE, mew=1.5, ls='none',
            label=f'shallower eclipse, {SHALLOW:,} ppm (0.26%)')
    ax.legend(loc='lower left', ncol=4, bbox_to_anchor=(-0.01, 1.0), handletextpad=0.3, columnspacing=1.4)
    # period brackets
    (a, _), (b, _), (c, _) = tcs[0], tcs[1], tcs[2]

    def bracket(x0, x1, y, text):
        ax.plot([x0, x0, x1, x1], [y - 300, y, y, y - 300], color=INK2, lw=1.0, solid_capstyle='round')
        ax.text((x0 + x1) / 2, y + 180, text, ha='center', va='bottom', fontsize=9, color=INK)
    bracket(a, b, 1950, 'catalogued period 0.65 d')
    bracket(a, c, 3250, 'true orbit 1.30 d: one deep + one shallow eclipse')
    save(fig, '01_every_other_dip_is_deeper.png')


# ---- graph 2: folded eclipses ----------------------------------------------------------------------------
def graph_fold(lcs):
    fig = frame('The two kinds of dip, stacked',
                'All 60 dips from two TESS sectors, lined up on their mid-times: the odd dips are 2.2 times deeper '
                'than the even ones. A planet would block the same amount of light every time.',
                'TESS Sectors 28 and 68 (2020, 2023), TESS-SPOC light curves, 6-minute medians with 1σ errors. '
                'Reference lines: adopted depths from the campaign eclipse model (±8% between TESS data products).',
                h=5.6)
    ax = fig.add_axes([0.085, 0.19, 0.765, 0.555])
    nev = {}
    for which, col, lab in (('deep', BLUE, 'odd dips: deeper eclipse'), ('shallow', ORANGE, 'even dips: shallower eclipse')):
        dts, ys, evs = [], [], 0
        for lc in lcs:
            dt, n = dt_to(lc['t'], which)
            s = np.abs(dt) < 0.14
            dts.append(dt[s]); ys.append(lc['f'][s] - 1); evs += np.unique(n[s]).size
        dt, y = np.concatenate(dts) * 24, np.concatenate(ys) * 1e6
        nev[which] = evs
        edges = np.arange(-3.2, 3.2001, 0.1)
        k = np.digitize(dt, edges)
        xc, yc, ye = [], [], []
        for i in range(1, edges.size):
            v = y[k == i]
            if v.size > 10:
                xc.append(0.5 * (edges[i - 1] + edges[i])); yc.append(np.median(v))
                ye.append(1.2533 * 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(v.size))
        ax.errorbar(xc, yc, yerr=ye, color=col, lw=1.6, ms=4.5, marker='o', mec=SURFACE, mew=1.0,
                    elinewidth=0.8, capsize=0, label=f'{lab} ({evs} dips)', zorder=3)
    for d, txt in ((DEEP, f'{DEEP:,} ppm\n(0.58% dimmer)'), (SHALLOW, f'{SHALLOW:,} ppm\n(0.26% dimmer)')):
        ax.axhline(-d, color=MUTED, lw=0.8, zorder=1)
        ax.text(3.3, -d, txt, ha='left', va='center', fontsize=9, color=INK, clip_on=False, linespacing=1.3)
    ax.yaxis.set_major_formatter(PPM)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-6900, 900)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel('Hours from mid-eclipse')
    ax.set_ylabel('Brightness change (ppm)')
    ax.legend(loc='lower left', ncol=2, bbox_to_anchor=(-0.01, 1.0), columnspacing=1.6)
    save(fig, '02_two_kinds_of_dip.png')
    return nev


# ---- graph 3: every sector -------------------------------------------------------------------------------
def graph_sectors():
    """Per-sector depths from the campaign's own TESScut photometry (one pipeline for all six sectors,
    neighbour contamination removed): chromatic/v3/results/tesscut_all6.json = [deep, shallow, err_deep, err_shallow, ...]."""
    import json
    years = {1: 'Jul\u2013Aug 2018', 28: 'Jul\u2013Aug 2020', 68: 'Jul\u2013Aug 2023', 95: 'Jul\u2013Aug 2025',
             105: 'Jun\u2013Jul 2026', 106: 'Jul\u2013Aug 2026'}
    js = json.load(open(os.path.join(CAMPAIGN, 'partial', 'chromatic', 'v3', 'results', 'tesscut_all6.json')))
    rows = []
    for sec in sorted(js, key=int):
        dd, ds, ed, es = (v * 1e6 for v in js[sec][:4])
        r = ds / dd
        rows.append(dict(sector=int(sec), deep=dd, deep_e=ed, shallow=ds, shallow_e=es, ratio=r,
                         ratio_e=r * np.hypot(ed / dd, es / ds)))
        print(f"  S{sec}: deep {dd:.0f}+-{ed:.0f}, shallow {ds:.0f}+-{es:.0f}, ratio {r:.3f}")
    fig = frame('The same pattern in every TESS sector, 2018–2026',
                'Depth of the two eclipses in each sector. The deeper one is always about twice the shallower one; '
                'for a planet the two dots would coincide.',
                'Campaign photometry of the TESS full-frame images (one method for all six sectors, light from '
                'neighbouring stars removed), eclipse-template fits with 1σ errors (chromatic analysis, '
                'tesscut_all6.json).',
                h=5.6)
    ax = fig.add_axes([0.21, 0.19, 0.60, 0.59])
    for i, r in enumerate(rows[::-1]):
        ax.plot([r['shallow'], r['deep']], [i, i], color=AXIS, lw=2, solid_capstyle='round', zorder=2)
        ax.errorbar(r['deep'], i, xerr=r['deep_e'], fmt='o', color=BLUE, ms=8, mec=SURFACE, mew=1.5,
                    elinewidth=1, capsize=0, zorder=3)
        ax.errorbar(r['shallow'], i, xerr=r['shallow_e'], fmt='o', color=ORANGE, ms=8, mec=SURFACE, mew=1.5,
                    elinewidth=1, capsize=0, zorder=3)
        ax.text(7350, i, f"{r['ratio']:.2f} ± {r['ratio_e']:.2f}", ha='left', va='center', fontsize=9.5,
                color=INK, clip_on=False)
    ax.text(7350, len(rows) - 0.35, 'shallow / deep', ha='left', va='bottom', fontsize=8.5, color=INK2,
            clip_on=False)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"Sector {r['sector']}\n{years[r['sector']]}" for r in rows[::-1]], fontsize=9,
                       color=INK2, linespacing=1.3)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='y', visible=False)
    ax.set_xlim(0, 7000)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.xaxis.set_major_formatter(PPM)
    ax.set_xlabel('Eclipse depth (ppm)')
    ax.spines['left'].set_visible(False)
    ax.plot([], [], 'o', color=BLUE, ms=8, mec=SURFACE, mew=1.5, label='deeper eclipse')
    ax.plot([], [], 'o', color=ORANGE, ms=8, mec=SURFACE, mew=1.5, label='shallower eclipse')
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.12), ncol=2, handletextpad=0.2, columnspacing=1.5)
    save(fig, '03_same_pattern_every_sector.png')
    return rows


# ---- graph 4: Gaia RV ------------------------------------------------------------------------------------
def graph_rv():
    G, MSUN, DAY = 6.674e-11, 1.989e30, 86400.0
    m2 = np.linspace(0.08, 1.0, 400)
    K = (2 * np.pi * G / (P * DAY)) ** (1 / 3) * m2 * MSUN / ((M_G + m2) * MSUN) ** (2 / 3) / 1e3
    k01 = np.interp(0.1, m2, K)
    kmax = np.sqrt(2) * 2.63           # all of Gaia's per-epoch rms taken as orbital motion (circular orbit)
    fig = frame('The bright G star is not one of the eclipsing stars',
                'How fast the G star would swing back and forth if it were itself in the 1.30-day eclipsing pair, '
                'against what Gaia measured.',
                'Curve: semi-amplitude K for an edge-on circular orbit with P = 1.30367 d and a 1.03 M☉ G star. '
                'Gaia DR3: 15 velocities over 946 days (2014–17), RV = −6.07 ± 0.86 km/s, scatter consistent '
                'with noise (p = 0.28). Even if all of the 2.6 km/s rms scatter were orbital motion, K would be at '
                'most 3.7 km/s.', h=5.5)
    ax = fig.add_axes([0.085, 0.2, 0.86, 0.58])
    ax.fill_between([0, 1.05], 0, kmax, color='#cde2fb', lw=0, zorder=1)
    ax.text(0.985, kmax + 2.5, 'largest swing Gaia allows: 3.7 km/s', ha='right', va='bottom', fontsize=9,
            color=INK)
    ax.plot(m2, K, color=BLUE, lw=2, solid_capstyle='round', zorder=3)
    ax.plot(0.1, k01, 'o', color=BLUE, ms=8, mec=SURFACE, mew=1.5, zorder=4)
    ax.text(0.125, 10.5, f'even a 0.1 M☉ partner: {k01:.0f} km/s', fontsize=9.5, color=INK, ha='left',
            va='center')
    ax.text(0.5, 108, 'speed needed if the G star had an\neclipsing partner of this mass',
            ha='center', va='bottom', fontsize=9, color=INK2, linespacing=1.3)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 140)
    ax.set_xlabel('Mass of a hypothetical eclipsing partner (solar masses)')
    ax.set_ylabel('Orbital speed of the G star, K (km/s)')
    save(fig, '04_g_star_is_not_the_eclipsing_star.png')


# ---- graph 5: eclipse colour -----------------------------------------------------------------------------
def graph_colour():
    fig = frame('The eclipsed star is red, not Sun-like',
                'Depth of the deeper eclipse in visible and infrared light, as a multiple of its depth in TESS '
                '(red-optical) light.',
                'Measured: SuperWASP 2006–07 and NEOWISE W1 2010–24 (campaign colour analysis). Predictions: '
                'the adopted M2V + M4.5V pair (campaign synthesis) and a generic G/F eclipsed star of 5,600–7,000 K '
                '(campaign colour models). TESS = 1 by definition.', h=5.5)
    ax = fig.add_axes([0.085, 0.2, 0.62, 0.58])
    groups = [('Visible light\nSuperWASP, 0.4–0.7 µm', 0.42, 0.15, 0.33, (0.818, 1.310)),
              ('Infrared\nWISE W1, 3.4 µm', 3.7, 1.1, 3.3, (0.696, 1.148))]
    for i, (lab, meas, err, mpred, gf) in enumerate(groups):
        ax.add_patch(Rectangle((i - 0.30, gf[0]), 0.10, gf[1] - gf[0], fc=ORANGE, alpha=0.35, lw=0, zorder=2))
        ax.plot([i - 0.25] * 2, gf, color=ORANGE, lw=2, solid_capstyle='round', zorder=2)
        ax.plot(i - 0.02, mpred, marker='D', ms=8, color=BLUE, mec=SURFACE, mew=1.5, ls='none', zorder=3)
        ax.errorbar(i + 0.22, meas, yerr=err, fmt='o', color=INK, ms=8, mec=SURFACE, mew=1.5, elinewidth=1.2,
                    capsize=0, zorder=4)
        ax.text(i + 0.30, meas, f'{meas:.2f} ± {err:.2f}' if meas < 1 else f'{meas:.1f} ± {err:.1f}',
                ha='left', va='center', fontsize=9.5, color=INK)
    ax.axhline(1, color=MUTED, lw=0.8, zorder=1)
    ax.text(1.47, 1.0, 'same depth\nas in TESS', ha='right', va='bottom', fontsize=8.5, color=INK2)
    ax.set_yscale('log')
    ax.set_ylim(0.15, 7)
    ax.set_yticks([0.2, 0.5, 1, 2, 5])
    ax.set_yticklabels(['0.2×', '0.5×', '1×', '2×', '5×'])
    ax.minorticks_off()
    ax.set_xlim(-0.55, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[0] for g in groups], fontsize=9.5, color=INK2)
    ax.tick_params(axis='x', length=0)
    ax.grid(axis='x', visible=False)
    ax.set_ylabel('Depth relative to TESS')
    # legend outside, right
    lg = fig.add_axes([0.73, 0.2, 0.25, 0.58]); lg.axis('off')
    items = [('o', INK, 'measured'), ('D', BLUE, 'predicted: red M-dwarf pair\n(the adopted model)'),
             ('|', ORANGE, 'predicted: Sun-like (G/F)\neclipsed star, range')]
    for j, (m, c, t) in enumerate(items):
        y = 0.86 - j * 0.2
        if m == '|':
            lg.add_patch(Rectangle((0.03, y - 0.05), 0.04, 0.10, fc=ORANGE, alpha=0.35, lw=0,
                                   transform=lg.transAxes, clip_on=False))
            lg.plot([0.05, 0.05], [y - 0.05, y + 0.05], color=ORANGE, lw=2, transform=lg.transAxes, clip_on=False)
        else:
            lg.plot(0.05, y, marker=m, ms=8, color=c, mec=SURFACE, mew=1.5, transform=lg.transAxes, clip_on=False)
        lg.text(0.13, y, t, transform=lg.transAxes, va='center', fontsize=9, color=INK, linespacing=1.3)
    lg.text(0.0, 0.18, 'Shallower in visible light and\ndeeper in the infrared: the star\nbeing eclipsed is '
            'much cooler\nthan the Sun-like G star.', transform=lg.transAxes, fontsize=9, color=INK2, va='top')
    save(fig, '05_eclipsed_star_is_red.png')


# ---- graph 6: localisation -------------------------------------------------------------------------------
def graph_position():
    import csv
    rows = list(csv.DictReader(open(os.path.join(CAMPAIGN, 'partial', 'localise', 'gaia', 'neighbours_120arcsec.csv'))))
    ra0, de0 = float(rows[0]['ra']), float(rows[0]['dec'])
    stars = []
    for r in rows:
        dE = (float(r['ra']) - ra0) * np.cos(np.radians(de0)) * 3600
        dN = (float(r['dec']) - de0) * 3600
        stars.append((dE, dN, float(r['phot_g_mean_mag'])))
    fig = frame('The eclipses come from within 1.24″ of the G star',
                'A TESS pixel is 21″ wide, so TESS cannot see directly which star dims. '
                'The position of the dimming light narrows it down.',
                'Left: Gaia DR3 stars within 70″ (dot size by brightness), TESS pixel grid illustrative. Right: '
                'TESS difference images from six sectors, corrected for the colour shift in the TESS optics '
                '(campaign localisation). At 297 pc, 1″ = 297 au.', w=9.6, h=5.8)
    # left: wide field
    ax = fig.add_axes([0.06, 0.15, 0.40, 0.64])
    ax.set_aspect('equal')
    for k in range(-4, 5):
        ax.axvline(k * 21 + 10.5, color=GRID, lw=0.8, zorder=1)
        ax.axhline(k * 21 + 10.5, color=GRID, lw=0.8, zorder=1)
    ax.add_patch(Rectangle((-10.5, -10.5), 21, 21, fill=False, ec=INK2, lw=1.2, zorder=2))
    for dE, dN, g in stars:
        if abs(dE) < 70 and abs(dN) < 70:
            s = max(6, 260 * (10 ** (-0.4 * (g - 11.8))) ** 0.5)
            ax.scatter(dE, dN, s=s, color=INK if g < 12 else MUTED, edgecolors=SURFACE, linewidths=1.2, zorder=3)
    ax.grid(False)
    ax.set_xlim(70, -70); ax.set_ylim(-70, 70)
    ax.set_xticks([-60, -30, 0, 30, 60]); ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_xlabel('← East   offset (arcsec)')
    ax.set_ylabel('offset (arcsec)   North →')
    nb = [s for s in stars if 10 < np.hypot(s[0], s[1]) < 20][0]
    ax.annotate('TIC 71316387 (G = 11.8)', xy=(0, 0), xytext=(-14, 30), fontsize=8.5, color=INK, ha='left',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate(f'neighbour, 18″ (G = {nb[2]:.1f})\nruled out', xy=(nb[0], nb[1]), xytext=(60, -45),
                fontsize=8.5, color=INK, ha='left', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.text(10.5, -13, 'one TESS pixel', fontsize=8.5, color=INK2, ha='left', va='top')
    # right: zoom
    az = fig.add_axes([0.56, 0.15, 0.40, 0.64])
    az.set_aspect('equal')
    cE, cN, sE, sN = -0.17, 0.53, 0.18, 0.28            # colour-corrected best fit and empirical sd (arcsec)
    for dchi, col, lab in ((11.85, BLUE_STEPS[0], '99.7%'), (6.54, BLUE_STEPS[1], '95%'), (2.28, BLUE_STEPS[2], '68%')):
        az.add_patch(Ellipse((cE, cN), 2 * np.sqrt(dchi) * sE, 2 * np.sqrt(dchi) * sN, color=col, lw=0, zorder=2))
    for dchi, lab in ((11.85, '99.7%'), (6.54, '95%'), (2.28, '68%')):
        az.text(cE, cN + np.sqrt(dchi) * sN - 0.07, lab, ha='center', va='top', fontsize=7.5, color=INK, zorder=6)
    for r in (0.2, 0.6):
        az.add_patch(Circle((0, 0), r, fill=False, ec=INK2, lw=0.9, zorder=4))
    az.plot(0, 0, marker='*', ms=17, color=INK, mec=SURFACE, mew=1.2, zorder=5)
    az.plot(cE, cN, marker='x', ms=7, color=INK, mew=1.6, zorder=5)
    az.set_xlim(2, -2); az.set_ylim(-2, 2)
    az.set_xticks([-2, -1, 0, 1, 2]); az.set_yticks([-2, -1, 0, 1, 2])
    az.set_xlabel('← East   offset (arcsec)')
    az.set_ylabel('offset (arcsec)   North →')
    az.text(0, -0.25, 'G star', ha='center', va='top', fontsize=8.5, color=INK, zorder=6)
    az.text(cE - 0.14, cN, 'TESS best fit', ha='left', va='center', fontsize=8.5, color=INK, zorder=6)
    az.text(1.92, 1.92, 'Eclipse source (TESS, 6 sectors):\n95% within 1.24″ = 370 au', fontsize=8.5,
            color=INK, ha='left', va='top', linespacing=1.35)
    az.annotate('rings: most likely separation\n0.2–0.6″ (all evidence combined)', xy=(0.42, -0.42),
                xytext=(1.92, -1.35), fontsize=8.5, color=INK2, ha='left', va='center', linespacing=1.35,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    save(fig, '06_where_the_eclipses_come_from.png')


def main():
    spoc = [load(s) for s in (28, 68)]
    graph_alternation(spoc[1])
    nev = graph_fold(spoc)
    print('  fold events:', nev)
    graph_sectors()
    graph_rv()
    graph_colour()
    graph_position()


if __name__ == '__main__':
    main()
