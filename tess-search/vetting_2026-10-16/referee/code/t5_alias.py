"""TIC 61141412: (1) odd/even and phase-0.5 depths of the 6.5285-d dips
(after removing the measured ETV), (2) ellipsoidal/beaming photometric test of
the Gaia SB1 period (5.4524 d) against the dip period (6.5285 d), (3) search
for eclipses at the SB1 period, (4) per-event O-C listing with parity."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *

TIC = 61141412
OUT = {}
ETV = json.load(open(f'{BASE}/results/t5_etv.json'))
ep = ETV['best_win0.25_poly1_detNone']
print('per-event O-C (best products, linear ephemeris):')
for v in ep['events']:
    print(f"  S{v['sector']:3d} {v['prod']:9s} n={v['n']:4d} ({'odd' if v['n']%2 else 'even'}) tc={v['tc']:.5f} +- {v['err_min']:.1f} min  O-C={v['oc_lin_min']:+6.1f} min  depth={v['depth_ppm']:.0f} ppm S/N={v['snr']:.1f}")
c2 = ep['deg2']['coef']


def tpred(n):
    return c2[0] + c2[1] * n + c2[2] * n ** 2


QLP_SHIFT = {104: 387.2 / 86400.}
dsets = [(13, 'SPOC', 'PDCSAP_FLUX'), (27, 'TESS-SPOC', 'PDCSAP_FLUX'), (67, 'TESS-SPOC', 'PDCSAP_FLUX'), (94, 'QLP', 'SAP_FLUX'), (104, 'QLP', 'SAP_FLUX')]
T14, T12 = 1.0 / 24, 0.35 / 24


def depth_at(t, f, centres, win=0.25):
    blocks = []
    nW = 0
    for c in centres:
        s = np.abs(t - c) < win
        if s.sum() > 10 and np.sum(np.abs(t[s] - c) < T14 / 2) >= 2 and np.sum(t[s] < c - T14) >= 3 and np.sum(t[s] > c + T14) >= 3:
            blocks.append((c, s)); nW += 1
    if nW == 0:
        return np.nan, np.nan, 0
    X = np.zeros((0, 1 + 2 * nW)); Y = np.zeros(0)
    for k, (c, s) in enumerate(blocks):
        dt = t[s] - c
        texp = np.median(np.diff(t))
        sh = 1 - supersample(trapezoid, t[s], texp, t0=c, depth=1.0, T14=T14, T12=T12)
        Xk = np.zeros((s.sum(), 1 + 2 * nW)); Xk[:, 0] = -sh; Xk[:, 1 + 2 * k] = 1; Xk[:, 2 + 2 * k] = dt
        X = np.vstack([X, Xk]); Y = np.concatenate([Y, f[s]])
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ b
    cov = np.linalg.inv(X.T @ X) * np.var(r, ddof=X.shape[1])
    return b[0], np.sqrt(cov[0, 0]), nW


# (1) odd / even / phase 0.5 with ETV-corrected centres (quadratic ephemeris)
od, ev_, ph5 = [], [], []
for sec, prod, col in dsets:
    t, f, e, ex = load_lc(TIC, sec, prod, col)
    if prod == 'QLP' and sec in QLP_SHIFT:
        t = t + QLP_SHIFT[sec]
    n0 = int(np.floor((t.min() - c2[0]) / c2[1])) - 1
    n1 = int(np.ceil((t.max() - c2[0]) / c2[1])) + 1
    ns = np.arange(n0, n1 + 1)
    for n in ns:
        d, de, nw = depth_at(t, f, [tpred(n)])
        if nw:
            (od if n % 2 else ev_).append((d, de))
        d5, de5, nw5 = depth_at(t, f, [tpred(n + 0.5)])
        if nw5:
            ph5.append((d5, de5))


def wmean(a):
    a = np.array(a)
    w = 1 / a[:, 1] ** 2
    return np.sum(a[:, 0] * w) / np.sum(w), 1 / np.sqrt(np.sum(w)), len(a)


mo, so, no = wmean(od); me, se, ne = wmean(ev_); m5, s5, n5 = wmean(ph5)
OUT['odd'] = [mo, so, no]; OUT['even'] = [me, se, ne]; OUT['phase05'] = [m5, s5, n5]
print(f'\nodd depth {mo*1e6:.0f}+-{so*1e6:.0f} ppm (n={no}); even {me*1e6:.0f}+-{se*1e6:.0f} (n={ne}); diff {(me-mo)/np.hypot(so,se):+.1f} sigma')
print(f'phase-0.5 (half-period) depth {m5*1e6:.0f}+-{s5*1e6:.0f} ppm (n={n5})')

# (2) ellipsoidal / beaming test on SAP-type fluxes (preserve multi-day signals)
sap_sets = [(13, 'SPOC', 'SAP_FLUX'), (13, 'TESS-SPOC', 'SAP_FLUX'), (27, 'TESS-SPOC', 'SAP_FLUX'), (67, 'TESS-SPOC', 'SAP_FLUX'),
            (13, 'QLP', 'SAP_FLUX'), (27, 'QLP', 'SAP_FLUX'), (67, 'QLP', 'SAP_FLUX'), (94, 'QLP', 'SAP_FLUX'), (104, 'QLP', 'SAP_FLUX')]
for sec in (13, 27, 67, 94, 104):
    sap_sets.append((sec, 'CUT', '2.5'))
Ptest = {'SB1_P': 5.452382, 'SB1_P/2': 2.726191, 'dip_P': 6.52854, 'dip_P/2': 3.26427}
OUT['sinus'] = {}
for sec, prod, col in sap_sets:
    if prod == 'CUT':
        t, f, e, ex = aperture_lc(TIC, sec, r=float(col))
    else:
        r = load_lc(TIC, sec, prod, col)
        if r is None:
            continue
        t, f, e, ex = r
    # mask the dips
    nn = np.round((t - c2[0]) / c2[1])
    m = np.abs(t - tpred(nn)) > 0.06
    t, f = t[m], f[m]
    # bin to 30 min
    tb = np.arange(t.min(), t.max(), 1 / 48.)
    idx = np.digitize(t, tb)
    tt = np.array([t[idx == i].mean() for i in np.unique(idx)])
    ff = np.array([f[idx == i].mean() for i in np.unique(idx)])
    # segments split at gaps > 1 d; per-segment quadratic baseline
    seg = np.concatenate([[0], np.cumsum(np.diff(tt) > 1.0)])
    cols = []
    for s in np.unique(seg):
        w = seg == s
        x = (tt - tt[w].mean())
        for k in range(3):
            cols.append(np.where(w, x ** k, 0.0))
    res = {}
    for name, Pp in Ptest.items():
        X = np.vstack(cols + [np.sin(2 * np.pi * tt / Pp), np.cos(2 * np.pi * tt / Pp)]).T
        b, *_ = np.linalg.lstsq(X, ff, rcond=None)
        rr = ff - X @ b
        cov = np.linalg.inv(X.T @ X) * np.var(rr, ddof=X.shape[1])
        amp = np.hypot(b[-2], b[-1])
        aerr = np.sqrt((cov[-2, -2] + cov[-1, -1]) / 2)
        res[name] = (amp * 1e6, aerr * 1e6)
    # also LS periodogram peak in 1.5-8 d
    from astropy.timeseries import LombScargle
    Xb = np.vstack(cols).T
    bb, *_ = np.linalg.lstsq(Xb, ff, rcond=None)
    r0 = ff - Xb @ bb
    fr = np.linspace(1 / 8., 1 / 1.5, 3000)
    pw = LombScargle(tt, r0).power(fr)
    pk = 1 / fr[np.argmax(pw)]
    OUT['sinus'][f'S{sec}_{prod}_{col}'] = dict(res, LS_peak_P=pk, span=tt.max() - tt.min())
    print(f'S{sec:3d} {prod:9s} {col:9s} semi-amp ppm: ' + '  '.join(f'{k}={v[0]:4.0f}+-{v[1]:3.0f}' for k, v in res.items()) + f'   LS peak {pk:.2f} d')

# expected ellipsoidal semi-amplitudes
G = 6.674e-11; Ms = 1.989e30; Rs = 6.957e8
M1, R1 = 1.38, 2.6
for P in (5.452382, 6.52854):
    K = 18.13e3
    f_m = P * 86400 * K ** 3 / (2 * np.pi * G) / Ms
    M2 = 0.2
    for _ in range(50):
        M2 = (f_m * (M1 + M2) ** 2) ** (1 / 3.)
    a = (G * (M1 + M2) * Ms * (P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3.) / Rs
    A_ell = 1.26 * (M2 / M1) * (R1 / a) ** 3
    A_beam = 1.0 * 4 * K / 3e8  # beaming factor ~1 in TESS band for ~6000 K
    OUT[f'expected_P{P}'] = dict(M2_min=M2, a_Rsun=a, A_ell_ppm=A_ell * 1e6, A_beam_ppm=A_beam * 1e6, transit_prob=(R1 + 0.22) / a)
    print(f'P={P}: M2(sin i=1)={M2:.3f} Msun a={a:.1f} Rsun -> ellipsoidal semi-amp ~{A_ell*1e6:.0f} ppm at P/2, beaming ~{A_beam*1e6:.0f} ppm at P; eclipse prob {(R1+0.22)/a:.2f}')

# Gaia-sampling alias frequency
f1, f2 = 1 / 5.452382, 1 / 6.52854
print(f'frequency difference 1/5.452 - 1/6.5285 = {f1-f2:.5f} c/d -> {1/(f1-f2):.1f} d; Gaia precession 63.12 d harmonics: 1/63.12={1/63.12:.5f}, 2/63.12={2/63.12:.5f}')
OUT['alias_df'] = f1 - f2
save_json(OUT, f'{BASE}/results/t5_alias.json')
