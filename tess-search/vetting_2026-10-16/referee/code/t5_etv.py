"""TIC 61141412: independent per-event eclipse times with several products and
detrending choices; linear vs quadratic ephemeris; QLP S104 barycentric fix.
"""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from wotan import flatten
from reflib import *

TIC = 61141412
P0, T0 = 6.52855, 3143.000
OUT = {}
rng = np.random.default_rng(3)

# QLP S103/S104 timestamps lack the barycentric correction (verified against SPOC 2-min
# momentum-dump flags); add the TESScut TIMECORR for this target (387.2 s in S104).
QLP_SHIFT = {104: 387.2 / 86400.}


def get_data(sec, prod, col, detrend=None):
    if prod == 'CUT':
        t, f, e, ex = aperture_lc(TIC, sec, r=float(col))
    else:
        r = load_lc(TIC, sec, prod, col)
        if r is None:
            return None
        t, f, e, ex = r
        if prod == 'QLP' and sec in QLP_SHIFT:
            t = t + QLP_SHIFT[sec]
    if detrend is not None:
        fl = flatten(t, f, method='biweight', window_length=detrend, break_tolerance=0.3)
        m = np.isfinite(fl)
        t, f, e = t[m], fl[m], e[m]
    return t, f, e


def fit_event(t, f, e, tc0, shape, win, poly, texp):
    T14, T12 = shape
    s = np.abs(t - tc0) < win
    if s.sum() < 8:
        return None
    ts, fs, es = t[s], f[s], e[s]
    if np.sum(np.abs(ts - tc0) < T14 / 2 + 0.05) < 3 or np.sum(ts < tc0 - T14 / 2) < 3 or np.sum(ts > tc0 + T14 / 2) < 3:
        return None

    def design(tc):
        cols = [np.ones_like(ts)]
        dt = ts - tc0
        for k in range(1, poly + 1):
            cols.append(dt ** k)
        sh = 1 - supersample(trapezoid, ts, texp, t0=tc, depth=1.0, T14=T14, T12=T12)
        cols.append(-sh)
        return np.vstack(cols).T

    # grid scan in tc, linear LSQ for baseline + depth at each tc
    grid = tc0 + np.arange(-0.08, 0.0801, 0.001)
    chi = []
    for tc in grid:
        X = design(tc)
        b, *_ = np.linalg.lstsq(X / es[:, None], fs / es, rcond=None)
        chi.append(np.sum(((fs - X @ b) / es) ** 2))
    chi = np.array(chi)
    j = int(np.argmin(chi))
    # parabola refinement around minimum using finer grid
    fine = grid[j] + np.arange(-0.002, 0.00201, 0.0001)
    chif = []
    for tc in fine:
        X = design(tc)
        b, *_ = np.linalg.lstsq(X / es[:, None], fs / es, rcond=None)
        chif.append(np.sum(((fs - X @ b) / es) ** 2))
    chif = np.array(chif)
    k = int(np.argmin(chif))
    tcb = fine[k]
    X = design(tcb)
    b, *_ = np.linalg.lstsq(X / es[:, None], fs / es, rcond=None)
    res = fs - X @ b
    red = np.sum((res / es) ** 2) / (len(fs) - X.shape[1] - 1)
    # timing error: delta-chi2 = 1 * max(red,1) from the scan curve (grid + fine)
    allg = np.concatenate([grid, fine]); allc = np.concatenate([chi, chif])
    o = np.argsort(allg); allg, allc = allg[o], allc[o]
    thr = allc.min() + max(red, 1.0)
    inside = allg[allc <= thr]
    err = (inside.max() - inside.min()) / 2 if len(inside) > 1 else 0.0005
    # prayer-bead (cyclic residual shifts) for red noise
    pb = []
    model = X @ b
    for sh in np.linspace(0, len(res), 25, endpoint=False).astype(int)[1:]:
        fsim = model + np.roll(res, sh)
        cc = []
        for tc in grid[max(0, j - 30):j + 31]:
            Xs = design(tc)
            bb, *_ = np.linalg.lstsq(Xs / es[:, None], fsim / es, rcond=None)
            cc.append(np.sum(((fsim - Xs @ bb) / es) ** 2))
        pb.append(grid[max(0, j - 30):j + 31][int(np.argmin(cc))])
    pb_err = np.std(pb)
    # depth error
    cov = np.linalg.inv((X / es[:, None]).T @ (X / es[:, None])) * max(red, 1)
    return dict(tc=tcb, err=max(err, pb_err), err_formal=err, err_pb=pb_err, depth=b[-1], depth_err=np.sqrt(cov[-1, -1]),
                redchi2=red, n=int(s.sum()), at_edge=bool(j in (0, len(grid) - 1)))


def run(dsets, shape, win, poly, detrend):
    evs = []
    for sec, prod, col in dsets:
        d = get_data(sec, prod, col, detrend=detrend)
        if d is None:
            continue
        t, f, e = d
        texp = np.median(np.diff(t))
        n0, n1 = int(np.floor((t.min() - T0) / P0)), int(np.ceil((t.max() - T0) / P0))
        for n in range(n0, n1 + 1):
            tc0 = T0 + n * P0
            if tc0 < t.min() - 0.05 or tc0 > t.max() + 0.05:
                continue
            r = fit_event(t, f, e, tc0, shape, win, poly, texp)
            if r is None:
                continue
            r.update(sector=sec, prod=prod, col=col, n=n)
            evs.append(r)
    return evs


def ephem(evs, snr_min=3.0, max_err=0.01):
    good = [v for v in evs if v['depth'] / v['depth_err'] > snr_min and v['err'] < max_err and not v['at_edge']]
    n = np.array([v['n'] for v in good], float)
    tc = np.array([v['tc'] for v in good]); s = np.array([v['err'] for v in good])
    out = {}
    for deg in (1, 2):
        A = np.vander(n, deg + 1, increasing=True)
        W = 1 / s
        c, *_ = np.linalg.lstsq(A * W[:, None], tc * W, rcond=None)
        res = tc - A @ c
        chi2 = float(np.sum((res / s) ** 2))
        out[f'deg{deg}'] = dict(coef=c.tolist(), chi2=chi2, dof=len(tc) - deg - 1, res_min=(res * 1440).tolist())
    # O-C wrt linear fit
    c1 = out['deg1']['coef']
    out['events'] = [dict(sector=v['sector'], prod=v['prod'], n=int(v['n']), tc=v['tc'], err_min=v['err'] * 1440, depth_ppm=v['depth'] * 1e6,
                          snr=v['depth'] / v['depth_err'], oc_lin_min=(v['tc'] - (c1[0] + c1[1] * v['n'])) * 1440) for v in good]
    return out


# ---------------- products: "best" per sector (as in our verdict) and fully independent TESScut photometry
best = [(13, 'SPOC', 'PDCSAP_FLUX'), (27, 'TESS-SPOC', 'PDCSAP_FLUX'), (67, 'TESS-SPOC', 'PDCSAP_FLUX'), (94, 'QLP', 'SAP_FLUX'), (104, 'QLP', 'SAP_FLUX')]
cut = [(13, 'CUT', '2.0'), (27, 'CUT', '2.0'), (67, 'CUT', '2.0'), (94, 'CUT', '2.0'), (104, 'CUT', '2.0')]
qlp = [(13, 'QLP', 'SAP_FLUX'), (27, 'QLP', 'SAP_FLUX'), (67, 'QLP', 'SAP_FLUX'), (94, 'QLP', 'SAP_FLUX'), (104, 'QLP', 'SAP_FLUX')]
shape = (1.0 / 24, 0.35 / 24)
configs = []
for name, ds in [('best', best), ('tesscut_r2', cut), ('qlp_all', qlp)]:
    for win, poly, det in [(0.25, 1, None), (0.25, 2, None), (0.4, 2, None), (0.15, 1, None), (0.3, 0, 0.5), (0.3, 0, 1.0)]:
        configs.append((name, ds, win, poly, det))
summary = []
for name, ds, win, poly, det in configs:
    evs = run(ds, shape, win, poly, det)
    ep = ephem(evs)
    key = f'{name}_win{win}_poly{poly}_det{det}'
    OUT[key] = ep
    bysec = {}
    for v in ep['events']:
        bysec.setdefault(v['sector'], []).append(v['oc_lin_min'])
    line = (f"{key:34s} nev={len(ep['events']):2d} lin chi2={ep['deg1']['chi2']:7.1f}/{ep['deg1']['dof']:2d}  quad chi2={ep['deg2']['chi2']:6.1f}/{ep['deg2']['dof']:2d}  "
            f"P={ep['deg1']['coef'][1]:.6f}  O-C(min) " + ' '.join(f"S{s}:{np.mean(v):+5.1f}" for s, v in sorted(bysec.items())))
    print(line)
    summary.append(line)
save_json(OUT, f'{BASE}/results/t5_etv.json')
