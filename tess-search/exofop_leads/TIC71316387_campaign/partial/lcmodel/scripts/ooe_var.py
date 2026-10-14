"""Out-of-eclipse variability at P_EB and P_EB/2.
Model per sector: per-segment Legendre polynomial trend (degree ~ duration/2.5 d) + harmonics of P:
   1 + a1 cos(2 pi phi) + b1 sin(2 pi phi) + a2 cos(4 pi phi) + b2 sin(4 pi phi),  phi=0 at primary eclipse.
Eclipses masked (+-0.06 d = +-1.45 h around both conjunctions; T14/2 = 0.032 d).
Joint fit over sectors (harmonic amplitudes common), plus per-sector fits.
Empirical null: identical fit at 400 fake periods in 0.8-3.0 d (avoiding P, P/2, 2P, 3P/2 +-3%)."""
import numpy as np, sys, json, glob
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import lclib as L
from astropy.io import fits
from numpy.polynomial import legendre as LG
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P, T0 = 1.3036712, 3165.6697
RD = L.RD
def load_spoc_col(sec, col):
    f = glob.glob(f'{RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
    d = fits.open(f)[1].data
    t = np.array(d['TIME'], float); fl = np.array(d[col], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    return dict(t=t[m], f=fl[m] / np.median(fl[m]), sector=sec, prod=col, texp=L.EXPT[sec])
def prep(lc, mask_hw=0.06):
    t, f = lc['t'], lc['f']
    ph = ((t - T0) / P) % 1
    dd = np.minimum(np.minimum(ph, 1 - ph), np.abs(ph - 0.5)) * P
    m = dd > mask_hw
    t, f = t[m], f[m]
    # 5-sigma clip vs 0.5-d running median
    import pandas as pd
    rm = pd.Series(f, index=pd.to_datetime(t, unit='D')).rolling('0.5D', center=True, min_periods=5).median().values
    r = f - rm; s = 1.4826 * np.nanmedian(np.abs(r))
    k = ~(np.abs(r) > 5 * s)
    t, f = t[k], f[k]
    # segments: split at gaps > 0.5 d
    br = np.r_[0, np.where(np.diff(t) > 0.5)[0] + 1, t.size]
    seg = np.zeros(t.size, int)
    for i in range(len(br) - 1): seg[br[i]:br[i + 1]] = i
    return t, f, seg
def design_trend(t, seg, sec_id, dper=2.5):
    cols = []
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < 20: continue
        ts = t[m]; x = 2 * (ts - ts.min()) / (ts.max() - ts.min()) - 1
        deg = int(max(2, np.ceil((ts.max() - ts.min()) / dper)))
        for j in range(deg + 1):
            c = np.zeros(t.size); c[m] = LG.legval(x, np.eye(deg + 1)[j]); cols.append(c)
    return np.array(cols).T
def harm(t, per, t0=T0):
    ph = 2 * np.pi * (t - t0) / per
    return np.vstack([np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)]).T
_cache = {}
def fit(datasets, per, return_resid=False):
    """trend projected out once per dataset combination (cached), then 4-parameter LSQ"""
    key = id(datasets[0][0]), len(datasets)
    if key not in _cache and not return_resid:
        A_tr = []; Y = []; W = []
        for (t, f, seg, sig) in datasets:
            A_tr.append(design_trend(t, seg, 0)); Y.append(f - 1); W.append(np.full(t.size, 1 / sig))
        ntot = sum(a.shape[1] for a in A_tr); N = sum(len(y) for y in Y)
        A = np.zeros((N, ntot)); i0 = 0; j0 = 0
        for a in A_tr:
            n = a.shape[0]; A[i0:i0 + n, j0:j0 + a.shape[1]] = a; i0 += n; j0 += a.shape[1]
        w = np.concatenate(W); y = np.concatenate(Y)
        Q, R = np.linalg.qr(A * w[:, None])
        yw = y * w; yp = yw - Q @ (Q.T @ yw)
        tt = np.concatenate([d[0] for d in datasets])
        _cache[key] = (Q, w, yp, tt)
    if not return_resid:
        Q, w, yp, tt = _cache[key]
        H = harm(tt, per) * w[:, None]; Hp = H - Q @ (Q.T @ H)
        c, *_ = np.linalg.lstsq(Hp, yp, rcond=None)
        cov = np.linalg.inv(Hp.T @ Hp)
        return c, np.sqrt(np.diag(cov))
    return fit_full(datasets, per, return_resid)
def fit_full(datasets, per, return_resid=False):
    """datasets: list of (t, f, seg, sigma). Joint LSQ; returns harmonic coefs and formal errors"""
    A_tr = []; H = []; Y = []; W = []
    ntr = [0]
    for (t, f, seg, sig) in datasets:
        A_tr.append(design_trend(t, seg, 0)); H.append(harm(t, per)); Y.append(f - 1); W.append(np.full(t.size, 1 / sig))
    ntot = sum(a.shape[1] for a in A_tr); N = sum(len(y) for y in Y)
    A = np.zeros((N, ntot + 4)); i0 = 0; j0 = 0
    for a, h in zip(A_tr, H):
        n = a.shape[0]; A[i0:i0 + n, j0:j0 + a.shape[1]] = a; A[i0:i0 + n, ntot:] = h; i0 += n; j0 += a.shape[1]
    y = np.concatenate(Y); w = np.concatenate(W)
    Aw = A * w[:, None]; yw = y * w
    c, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    coef = c[ntot:]; err = np.sqrt(np.diag(cov)[ntot:])
    if return_resid:
        return coef, err, y - A @ c, A[:, ntot:] @ coef
    return coef, err
def amps(c):
    return np.hypot(c[0], c[1]), np.hypot(c[2], c[3])
if __name__ == '__main__':
    sets = {}
    for sec in (1, 28, 68):
        sets[(sec, 'PDCSAP')] = load_spoc_col(sec, 'PDCSAP_FLUX')
        sets[(sec, 'SAP')] = load_spoc_col(sec, 'SAP_FLUX')
    for sec in (1, 28, 68, 95, 105, 106):
        lc = L.load_tc(sec); lc['prod'] = 'TC'; sets[(sec, 'TC')] = lc
    prepped = {}
    for k, lc in sets.items():
        t, f, seg = prep(lc)
        # white noise per point from point-to-point scatter
        sig = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
        prepped[k] = (t, f, seg, sig)
    rng = np.random.default_rng(42)
    fake = []
    while len(fake) < 400:
        pf = rng.uniform(0.8, 3.0)
        if min(abs(pf / x - 1) for x in (P, P / 2, 2 * P, 1.5 * P)) > 0.03: fake.append(pf)
    fake = np.array(fake)
    out = {}
    combos = {'default(PDC S1/28/68 + TC 95/105/106)': [(1, 'PDCSAP'), (28, 'PDCSAP'), (68, 'PDCSAP'), (95, 'TC'), (105, 'TC'), (106, 'TC')],
              'TC all six': [(s, 'TC') for s in (1, 28, 68, 95, 105, 106)],
              'SAP S1/28/68 + TC 95/105/106': [(1, 'SAP'), (28, 'SAP'), (68, 'SAP'), (95, 'TC'), (105, 'TC'), (106, 'TC')]}
    for s in (1, 28, 68): combos[f'S{s} PDCSAP'] = [(s, 'PDCSAP')]; combos[f'S{s} SAP'] = [(s, 'SAP')]
    for s in (1, 28, 68, 95, 105, 106): combos[f'S{s} TC'] = [(s, 'TC')]
    for name, keys in combos.items():
        ds = [prepped[k] for k in keys]
        c, e = fit(ds, P)
        A1, A2 = amps(c)
        nullc = np.array([fit(ds, pf)[0] for pf in fake])
        n1 = np.hypot(nullc[:, 0], nullc[:, 1]); n2 = np.hypot(nullc[:, 2], nullc[:, 3])
        # empirical per-coefficient scatter (robust) from null
        esc = 1.4826 * np.median(np.abs(nullc - np.median(nullc, 0)), 0)
        res = dict(coef_ppm=(c * 1e6).round(1).tolist(), formal_err_ppm=(e * 1e6).round(1).tolist(), emp_err_ppm=(esc * 1e6).round(1).tolist(),
                   A1_ppm=A1 * 1e6, A2_ppm=A2 * 1e6, p_A1=float(np.mean(n1 >= A1)), p_A2=float(np.mean(n2 >= A2)),
                   ell_cos2_ppm=-c[2] * 1e6, ell_z_emp=float(-c[2] / esc[2]),
                   A2_null95_ppm=float(np.percentile(n2, 95) * 1e6), A1_null95_ppm=float(np.percentile(n1, 95) * 1e6),
                   # one-sided 95% upper limit on ellipsoidal amplitude (-a2): measured + 1.645 * empirical sigma
                   ell_ul95_ppm=float((max(-c[2], 0) + 1.645 * esc[2]) * 1e6),
                   ell_lo_ppm=float((-c[2] - 1.645 * esc[2]) * 1e6))
        out[name] = res
        print(f"{name:42s} a1={c[0]*1e6:7.1f} b1={c[1]*1e6:7.1f} a2={c[2]*1e6:7.1f} b2={c[3]*1e6:7.1f} ppm | formal err {e[0]*1e6:5.1f} emp err {esc*1e6} | A1={A1*1e6:6.1f} (p={res['p_A1']:.3f}) A2={A2*1e6:6.1f} (p={res['p_A2']:.3f}) ell(-a2)={-c[2]*1e6:6.1f} UL95={res['ell_ul95_ppm']:.0f}", flush=True)
    json.dump(out, open(L.HERE + '/results/ooe_var.json', 'w'), indent=1)
    # plot binned folded OOE residual (trend removed) for default combo
    keys = combos['default(PDC S1/28/68 + TC 95/105/106)']
    ds = [prepped[k] for k in keys]
    c, e, r, hm = fit(ds, P, return_resid=True)
    tt = np.concatenate([d[0] for d in ds]); ph = ((tt - T0) / P) % 1
    y = r + hm
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    b = np.linspace(0, 1, 41); idx = np.digitize(ph, b)
    bx = [(b[j - 1] + b[j]) / 2 for j in range(1, len(b)) if (idx == j).sum() > 20]
    by = [np.mean(y[idx == j]) * 1e6 for j in range(1, len(b)) if (idx == j).sum() > 20]
    be = [np.std(y[idx == j]) / np.sqrt((idx == j).sum()) * 1e6 for j in range(1, len(b)) if (idx == j).sum() > 20]
    ax.errorbar(bx, by, be, fmt='o', ms=3)
    xx = np.linspace(0, 1, 300); ax.plot(xx, (harm(T0 + xx * P, P) @ c) * 1e6, 'r-')
    ax.axvspan(0, 0.06 / P, color='0.85'); ax.axvspan(1 - 0.06 / P, 1, color='0.85'); ax.axvspan(0.5 - 0.06 / P, 0.5 + 0.06 / P, color='0.85')
    ax.set_xlabel('phase (0 = primary eclipse)'); ax.set_ylabel('detrended flux [ppm]'); ax.set_title('TIC 71316387 out-of-eclipse, folded at P_EB (default set)')
    plt.tight_layout(); plt.savefig(L.HERE + '/plots/ooe_fold.png', dpi=110)
