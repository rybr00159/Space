"""TIC 90086266: independent per-event eclipse fits, ephemeris / cycle count,
odd-even depths, secondary eclipse, and per-event instrumental diagnostics.
Usage: python -I t1_events.py
"""
import sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from reflib import *

TIC = 90086266
P0, T00 = 2.339469, 4180.4740
RES = {}


def my_aperture_lc(sector, r=1.5, prod_tag='CUT'):
    cut = load_cut(TIC, sector)
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    tgt = g.iloc[0]
    x, y = cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])
    ny, nx = cut['flux'].shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.hypot(xx - x, yy - y)
    ap = rr <= r
    sky = (rr > 6) & (rr < 9)
    fl = cut['flux']
    bk = np.nanmedian(fl[:, sky], axis=1)
    f = np.nansum(fl[:, ap], axis=1) - bk * ap.sum()
    # centroids (flux-weighted, in aperture of radius 2.5)
    ap2 = rr <= 2.5
    w = fl[:, ap2] - bk[:, None]
    cx = np.sum(w * xx[ap2], axis=1) / np.sum(w, axis=1)
    cy = np.sum(w * yy[ap2], axis=1) / np.sum(w, axis=1)
    med = np.median(f)
    extra = dict(SAP_BKG=bk + np.nanmedian(cut['bkg'][:, ap], axis=1), MOM_CENTR1=cx, MOM_CENTR2=cy,
                 POS_CORR1=cut['pc1'], POS_CORR2=cut['pc2'], col='CUT_r%.1f' % r, median=med)
    e = np.full_like(f, robust_std(np.diff(f / med)) / np.sqrt(2))
    return cut['time'], f / med, e, extra


def texp_of(t):
    return np.median(np.diff(t))


def model_events(p, t, ev_idx, ep_list, texp, shape):
    """p = [tc_k, depth_k, a0_k, a1_k, a2_k for each event] with fixed shape (T14,T12)."""
    T14, T12 = shape
    m = np.empty_like(t)
    for k, ep in enumerate(ep_list):
        s = ev_idx == k
        tc, dep, a0, a1, a2 = p[5 * k:5 * k + 5]
        dt = t[s] - tc
        base = a0 + a1 * dt + a2 * dt ** 2
        m[s] = base * supersample(trapezoid, t[s], texp, t0=tc, depth=dep, T14=T14, T12=T12)
    return m


def fit_sector(t, f, e, P, T0, win=0.30, shape=None, poly=2):
    """Global shape fit then per-event fits. Returns dict."""
    texp = texp_of(t)
    ep = np.round((t - T0) / P)
    ph = t - (T0 + ep * P)
    sel = np.abs(ph) < win
    eps = []
    for E in np.unique(ep[sel]):
        s = sel & (ep == E)
        # require coverage on both sides and in eclipse
        if (np.sum(np.abs(ph[s]) < 0.04) >= max(2, int(0.02 / texp))) and np.sum(ph[s] < -0.1) >= 3 and np.sum(ph[s] > 0.1) >= 3:
            eps.append(E)
    eps = np.array(eps)
    # assemble arrays
    tt, ff, ee, idx = [], [], [], []
    for k, E in enumerate(eps):
        s = sel & (ep == E)
        tt.append(t[s]); ff.append(f[s]); ee.append(e[s]); idx.append(np.full(s.sum(), k))
    tt, ff, ee, idx = map(np.concatenate, (tt, ff, ee, idx))
    # stage 1: global shape fit (common depth, T14, T12, T0, P) with per-event quadratic baselines
    nE = len(eps)

    def resid_global(q):
        T0g, Pg, dep, T14, T12 = q[:5]
        m = np.empty_like(tt)
        for k, E in enumerate(eps):
            s = idx == k
            tc = T0g + E * Pg
            a0, a1, a2 = q[5 + 3 * k:8 + 3 * k]
            dt = tt[s] - tc
            m[s] = (a0 + a1 * dt + a2 * dt ** 2) * supersample(trapezoid, tt[s], texp, t0=tc, depth=dep, T14=T14, T12=T12)
        return (ff - m) / ee
    q0 = [T0, P, 0.03, 0.13, 0.05] + [1, 0, 0] * nE
    lo = [T0 - 0.02, P - 0.002, 0.0, 0.05, 0.005] + [0.9, -1, -10] * nE
    hi = [T0 + 0.02, P + 0.002, 0.2, 0.30, 0.15] + [1.1, 1, 10] * nE
    rg = least_squares(resid_global, q0, bounds=(lo, hi), x_scale='jac')
    T0g, Pg, depg, T14g, T12g = rg.x[:5]
    chi2g = np.sum(rg.fun ** 2)
    if shape is None:
        shape = (T14g, T12g)
    # stage 2: per-event fits with fixed shape
    out = []
    for k, E in enumerate(eps):
        s = idx == k
        ts, fs, es = tt[s], ff[s], ee[s]
        tc0 = T0g + E * Pg

        def r1(q):
            tc, dep, a0, a1, a2 = q
            dt = ts - tc
            base = a0 + a1 * dt + (a2 * dt ** 2 if poly == 2 else 0)
            return (fs - base * supersample(trapezoid, ts, texp, t0=tc, depth=dep, T14=shape[0], T12=shape[1])) / es
        r = least_squares(r1, [tc0, depg, 1, 0, 0], x_scale='jac')
        J = r.jac
        red = np.sum(r.fun ** 2) / max(1, len(ts) - 5)
        try:
            cov = np.linalg.inv(J.T @ J) * max(red, 1.0)
            err = np.sqrt(np.diag(cov))
        except np.linalg.LinAlgError:
            err = np.full(5, np.nan)
        out.append(dict(epoch=int(E), tc=r.x[0], tc_err=err[0], depth=r.x[1], depth_err=err[1], redchi2=red, n=len(ts)))
    return dict(T0=T0g, P=Pg, depth=depg, T14=T14g, T12=T12g, chi2=chi2g, n=len(tt), events=out, shape=shape)


def linfit(x, y, s):
    A = np.vstack([np.ones_like(x), x]).T
    W = 1 / s ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    b = C @ (A.T @ (W * y))
    chi2 = np.sum(((y - A @ b) / s) ** 2)
    return b, np.sqrt(np.diag(C)), chi2


datasets = []
for (sec, prod, col) in [(13, 'TESS-SPOC', 'SAP_FLUX'), (13, 'TESS-SPOC', 'PDCSAP_FLUX'), (13, 'QLP', 'SAP_FLUX'),
                         (104, 'QLP', 'SAP_FLUX'), (104, 'QLP', 'DET_FLUX')]:
    r = load_lc(TIC, sec, prod, col)
    datasets.append((sec, prod + ':' + col, r))
for sec in (13, 104):
    for rad in (1.0, 2.0):
        datasets.append((sec, 'CUT:r%.1f' % rad, my_aperture_lc(sec, r=rad)))

for sec, name, r in datasets:
    t, f, e, ex = r
    fit = fit_sector(t, f, e, P0, T00 + np.round((np.median(t) - T00) / P0) * P0)
    ev = fit['events']
    E = np.array([v['epoch'] for v in ev], float)
    tc = np.array([v['tc'] for v in ev])
    tce = np.array([v['tc_err'] for v in ev])
    dep = np.array([v['depth'] for v in ev])
    depe = np.array([v['depth_err'] for v in ev])
    b, be, chi2 = linfit(E, tc, tce)
    # odd/even with epoch numbering relative to T00 (S104 first event = 0)
    Eabs = np.round((tc - T00) / P0).astype(int)
    odd = Eabs % 2 == 1
    w = 1 / depe ** 2
    d_odd = np.sum(dep[odd] * w[odd]) / np.sum(w[odd]); s_odd = 1 / np.sqrt(np.sum(w[odd]))
    d_even = np.sum(dep[~odd] * w[~odd]) / np.sum(w[~odd]); s_even = 1 / np.sqrt(np.sum(w[~odd]))
    # empirical errors from scatter within each group
    se_odd = np.std(dep[odd], ddof=1) / np.sqrt(odd.sum()) if odd.sum() > 1 else np.nan
    se_even = np.std(dep[~odd], ddof=1) / np.sqrt((~odd).sum()) if (~odd).sum() > 1 else np.nan
    so = max(s_odd, se_odd); sev = max(s_even, se_even)
    diff = d_even - d_odd
    sig = diff / np.hypot(so, sev)
    RES[f'S{sec}_{name}'] = dict(T0=b[0] + 0, P=b[1], P_err=be[1], T0_epoch0=b[0], T0_err=be[0], chi2_lin=chi2, nev=len(ev),
                                 global_fit=dict(T0=fit['T0'], P=fit['P'], depth=fit['depth'], T14_h=fit['T14'] * 24, T12_h=fit['T12'] * 24),
                                 depth_even=d_even, depth_even_err=sev, depth_odd=d_odd, depth_odd_err=so,
                                 even_minus_odd_frac=diff / np.mean(dep), oe_sigma=sig,
                                 events=[dict(E=int(a), tc=c, tc_err=ce, depth=dd, depth_err=de) for a, c, ce, dd, de in zip(Eabs, tc, tce, dep, depe)])
    print(f'S{sec} {name:22s} nev={len(ev):2d} P={b[1]:.6f}+-{be[1]:.6f} T0(E0)={b[0]:.5f} chi2={chi2:.1f}/{len(ev)-2} '
          f'depth={fit["depth"]:.4f} T14={fit["T14"]*24:.2f}h T12={fit["T12"]*24:.2f}h  even={d_even:.4f}+-{sev:.4f} odd={d_odd:.4f}+-{so:.4f} '
          f'(E-O)/mean={diff/np.mean(dep)*100:.1f}% sig={sig:.1f}')

# cycle count between sectors using per-sector ephemerides from the best products
a = RES['S13_TESS-SPOC:PDCSAP_FLUX']; b_ = RES['S104_QLP:SAP_FLUX']
# reference times near the middle of each sector
gap = b_['T0_epoch0'] - a['T0_epoch0']  # both referenced to epochs counted from T00
# epochs in a are absolute (relative to T00) only through E numbering; recompute T at E=0 for each
cyc = []
for N in range(1070, 1087):
    Pn = (b_['events'][0]['tc'] - a['events'][0]['tc']) / (b_['events'][0]['E'] - a['events'][0]['E'] + (N - 1078))
    cyc.append((N, Pn, (Pn - a['P']) / a['P_err'], (Pn - b_['P']) / b_['P_err']))
RES['cycle_count'] = [dict(N=c[0], P=c[1], dev_S13_sigma=c[2], dev_S104_sigma=c[3]) for c in cyc]
print('\ncycle count test (N between first S13 and first S104 events, relative to 1078 nominal):')
for c in cyc:
    print('  N=%d P=%.6f  dev vs S13-only P: %6.1f sig ; vs S104-only P: %6.1f sig' % c)

# joint linear ephemeris through all events of best products
allE = np.array([v['E'] for v in a['events']] + [v['E'] for v in b_['events']], float)
allT = np.array([v['tc'] for v in a['events']] + [v['tc'] for v in b_['events']])
allS = np.array([v['tc_err'] for v in a['events']] + [v['tc_err'] for v in b_['events']])
bb, bbe, chi2j = linfit(allE, allT, allS)
RES['joint_ephem'] = dict(T0=bb[0], T0_err=bbe[0], P=bb[1], P_err=bbe[1], chi2=chi2j, dof=len(allE) - 2)
print('joint ephemeris: T0=%.5f+-%.5f P=%.7f+-%.7f chi2=%.1f/%d' % (bb[0], bbe[0], bb[1], bbe[1], chi2j, len(allE) - 2))
save_json(RES, f'{BASE}/results/t1_events.json')
