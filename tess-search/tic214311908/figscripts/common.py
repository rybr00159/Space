"""Shared settings for the TIC 214311908 write-up figures (style, paths, adopted values, small helpers).

Nothing here changes any analysis folder: all inputs are read-only, all outputs go to ../figures/.
Colour palette: categorical slots blue / orange / aqua (validated colour-blind safe as a set of three, all pairs),
plus near-black ink for models and a violet fourth slot used only with a distinct marker/line style."""
import os
for _k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_k, '2')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the write-up folder (writeup/ or tic214311908/)
# analysis root: the folder above writeup/ in the working layout, or tic214311908/analysis/ in the repository
B = os.path.dirname(W) if os.path.isdir(os.path.join(os.path.dirname(W), 'clean')) else os.path.join(W, 'analysis')
FIG = os.path.join(W, 'figures')
CUBES = '/home/user/work/tesscut5/214311906'

# ---------------------------------------------------------------- palette and ink
BLUE, ORANGE, AQUA, VIOLET = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
AQUA_DARK = '#13865c'          # darker step of aqua for text / thin lines (contrast)
INK, INK2, MUTED, GRID, BASE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
C_INV, C_REF, C_AP = BLUE, ORANGE, AQUA

# ---------------------------------------------------------------- adopted values (referee3_clean / writeup)
T0_ADOPT, P_ADOPT = 3840.9150, 5.691506          # BTJD_TDB, d (referee3_clean README item 4)
T0_ERR, P_ERR, RHO_T0P = 0.0013, 0.000011, 0.35   # writeup/ephemeris_adopted.json
T0_CLEAN, P_CLEAN = 3840.91488, 5.6915020         # clean/out/ephemeris_final.json (combined)
F_M3V = 979.0                                     # e-/s, PRF units
# adopted transit model: referee3_clean q08 posterior medians (INV/REF averaged: b 0.52, u1 0.148, u2 0.542),
# a/R* 25.28 (clean posterior median), k chosen so that the mid-transit depth is the adopted 15.8 ppt
G = 6.6743e-8
K_ADOPT, B_ADOPT, AR_ADOPT = 0.1188, 0.522, 25.28
RHO_ADOPT = 3 * np.pi * AR_ADOPT ** 3 / (G * (P_ADOPT * 86400.0) ** 2)
U1_ADOPT, U2_ADOPT = 0.148, 0.542


def aR_of_rho(rho, per=P_ADOPT):
    return (rho * G * (per * 86400.0) ** 2 / (3 * np.pi)) ** (1 / 3.)


def style():
    plt.rcParams.update({
        'figure.dpi': 150, 'savefig.dpi': 150, 'font.size': 9.5, 'axes.titlesize': 10.5, 'axes.labelsize': 9.5,
        'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'legend.fontsize': 8.2, 'font.family': 'DejaVu Sans',
        'axes.edgecolor': BASE, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'axes.titlecolor': INK, 'text.color': INK, 'axes.linewidth': 0.8, 'axes.grid': True, 'grid.color': GRID,
        'grid.linewidth': 0.6, 'axes.axisbelow': True, 'legend.frameon': False, 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'xtick.major.size': 3, 'ytick.major.size': 3,
        'lines.linewidth': 1.6, 'errorbar.capsize': 0, 'axes.spines.top': False, 'axes.spines.right': False})


def save(fig, name):
    p = f'{FIG}/{name}'
    fig.savefig(p, dpi=150, bbox_inches='tight', pad_inches=0.08)
    # keep files small: re-quantise to an 8-bit palette if larger than 380 kB
    try:
        from PIL import Image
        if os.path.getsize(p) > 380_000:
            im = Image.open(p).convert('RGB').quantize(colors=256, method=Image.Quantize.MEDIANCUT)
            im.save(p, optimize=True)
    except Exception as e:      # pragma: no cover
        print('quantise failed', e)
    print(name, os.path.getsize(p) // 1024, 'kB')


def panel_label(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight='bold', ha='right', va='bottom', color=INK)


def bin_curve(x, y, w, edges, min_n=3):
    """inverse-variance weighted bins; returns centres, mean, error (max of formal and scatter-based)"""
    xc, ym, ye = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x < b)
        if m.sum() < min_n: continue
        ww = w[m]; yy = y[m]
        mu = np.sum(ww * yy) / np.sum(ww)
        e_formal = 1 / np.sqrt(np.sum(ww))
        e_scatter = np.std(yy, ddof=1) / np.sqrt(m.sum())
        xc.append(np.mean(x[m])); ym.append(mu); ye.append(max(e_formal, e_scatter))
    return np.array(xc), np.array(ym), np.array(ye)


# ---------------------------------------------------------------- transit model with exposure integration
def transit_model(dt, k=K_ADOPT, b=B_ADOPT, rho=RHO_ADOPT, u1=U1_ADOPT, u2=U2_ADOPT, exp_s=None, per=P_ADOPT):
    """flux of a planet transit at times dt (d) from mid-transit; midpoint-rule exposure integration"""
    import batman
    dt = np.atleast_1d(np.asarray(dt, float))
    aR = aR_of_rho(rho, per)
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a = 0., per, k, aR
    p.inc = np.degrees(np.arccos(min(b / aR, 1.0))); p.ecc, p.w, p.limb_dark, p.u = 0., 90., 'quadratic', [u1, u2]
    if exp_s is None or exp_s < 60:
        return batman.TransitModel(p, dt).light_curve(p)
    ns = 5 if exp_s < 300 else (11 if exp_s < 900 else 31)
    off = ((np.arange(ns) + 0.5) / ns - 0.5) * exp_s / 86400.
    tt = (dt[:, None] + off[None, :]).ravel()
    return batman.TransitModel(p, tt).light_curve(p).reshape(-1, ns).mean(1)


def model_by_cadence(dt, cad, **kw):
    out = np.ones(len(dt))
    for c in np.unique(cad):
        m = cad == c
        out[m] = transit_model(dt[m], exp_s=float(c), **kw)
    return out


def fit_window_baseline(dt, y, w, s, order=1, extra=None):
    """weighted fit of y = poly(dt) + a*s (+ extra regressors); returns baseline (poly + extras) and a"""
    cols = [s] + [dt ** o for o in range(order + 1)]
    if extra is not None:
        cols += [extra[:, j] for j in range(extra.shape[1])]
    M = np.column_stack(cols)
    sw = np.sqrt(w)
    c, *_ = np.linalg.lstsq(M * sw[:, None], y * sw, rcond=None)
    base = M[:, 1:] @ c[1:]
    return base, c[0]
