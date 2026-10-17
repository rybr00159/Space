"""Referee transit-shape check from a PRF-weighted light curve at the M3V's apparent position.

Per real transit window (|t - tc| <= 0.45 d): each pixel detrended by a quadratic fitted to out-of-transit cadences
only (|t - tc| > 0.06 d, wider than half of any plausible T14); per cadence, the residual image within 4 px of the M3V
is fitted with f * PRF(M3V) + plane (weights 1/sigma_pix^2 from the OOT scatter) -> f(t) = flux lost by a source at
the M3V (e-/s, total-flux units).  Then a joint fit of all windows with an exposure-integrated trapezoid
(depth A e-/s, T14, flat fraction F = T23/T14, common time offset dT), white noise per sector from the OOT scatter,
inflated by a per-sector beta factor (time-averaged rms on 30-60 min bins vs white expectation).
Profiles: chi2 vs F (other parameters optimised), chi2 vs T14; V-shape (F=0) and box (F=1) comparisons.
Uncertainties also from a jackknife over transits.  Output: out/shape.json, figs/shape.png
"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

WM = json.load(open(f'{ROOT}/out/wcs_mira.json'))
g = gaia(); r3 = g[g.source_id == ID['M3V']]
LC = []
info = {}
for sec in SECTORS:
    S = Sec(sec); S.off = np.array(WM[str(sec)]['offset_px'])
    p3 = S.xy(r3)[0]; psf = S.prf()
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    m = (np.hypot(xx - p3[0], yy - p3[1]) <= 4.0).ravel()
    F = psf(*p3).ravel()[m]
    Pl = np.column_stack([F, np.ones(m.sum()), (xx.ravel() - p3[0])[m], (yy.ravel() - p3[1])[m]])
    gd = S.good; t = S.t[gd]; Y = S.flux[gd].reshape(gd.sum(), -1)[:, m]
    n = np.arange(np.ceil((t.min() - T0) / P) - 1, np.floor((t.max() - T0) / P) + 2)
    tr = T0 + n * P; tr = tr[(tr > t.min() - 0.2) & (tr < t.max() + 0.2)]
    beta_list = []
    for k, tc in zip(n[(T0 + n * P > t.min() - 0.2) & (T0 + n * P < t.max() + 0.2)], tr):
        w = np.abs(t - tc) <= 0.45
        if w.sum() < 20: continue
        dt = t[w] - tc
        oot = np.abs(dt) > 0.06
        inn = np.abs(dt) < 0.03
        if inn.sum() < max(1, 0.5 * 0.06 / S.dt) or (dt[oot] < 0).sum() < 0.3 * 0.39 / S.dt or (dt[oot] > 0).sum() < 0.3 * 0.39 / S.dt:
            continue
        X = np.column_stack([np.ones(w.sum()), dt, dt ** 2])
        B, *_ = np.linalg.lstsq(X[oot], Y[w][oot], rcond=None)
        R = Y[w] - X @ B
        sp = R[oot].std(0) + 1e-3
        Wt = 1 / sp ** 2
        A = Pl * np.sqrt(Wt)[:, None]
        coef = np.linalg.lstsq(A, (R * np.sqrt(Wt)[None, :]).T, rcond=None)[0]
        f = -coef[0]                                 # flux LOST (positive in transit)
        e = np.full_like(f, np.std(f[oot]))
        # beta: rms of OOT residual binned to ~1 h vs white expectation
        nb = max(1, int(round((1 / 24) / S.dt)))
        fo = f[oot][: (oot.sum() // nb) * nb]
        if nb > 1 and len(fo) >= 3 * nb:
            sb = fo.reshape(-1, nb).mean(1).std(); beta_list.append(sb / (np.std(f[oot]) / np.sqrt(nb)))
        LC.append(dict(sec=sec, epoch=int(k), tc=tc, t=t[w], f=f, e=e, dt=S.dt))
    info[sec] = dict(beta=float(np.median(beta_list)) if beta_list else 1.0, n_beta=len(beta_list))
    print(f'S{sec}: windows {sum(1 for L in LC if L["sec"]==sec)}; beta {info[sec]["beta"]:.2f}', flush=True)

for L in LC:
    L['e'] = L['e'] * max(1.0, info[L['sec']]['beta'])


def model(L, A, t14, fl, dT):
    return A * trapezoid(L['t'], L['tc'] + dT, L['dt'], t14=t14, flat=fl)


def chi2(p, LCs, fixF=None, fixT=None):
    A, t14h, fl, dTm = p
    if fixF is not None: fl = fixF
    if fixT is not None: t14h = fixT
    if not (0 <= fl <= 1 and 0.5 < t14h < 3.5): return 1e12
    return sum((((L['f'] - model(L, A, t14h / 24, fl, dTm / 1440)) / L['e']) ** 2).sum() for L in LCs)


def fit(LCs, fixF=None, fixT=None, fast=False):
    best = None
    for fl0 in ((0.5, 0.85) if fast else (0.2, 0.5, 0.8)):
        for t0h in ((1.4, 1.7) if fast else (1.4, 1.6, 1.8)):
            r = minimize(chi2, [14.0, t0h, fl0, 0.0], args=(LCs, fixF, fixT), method='Nelder-Mead',
                         options=dict(xatol=1e-4, fatol=1e-4, maxiter=4000))
            if best is None or r.fun < best.fun: best = r
    x = best.x.copy()
    if fixF is not None: x[2] = fixF
    if fixT is not None: x[1] = fixT
    return x, best.fun


xb, cb = fit(LC)
ndata = sum(len(L['f']) for L in LC)
out = dict(n_windows=len(LC), ndata=ndata, best=dict(A=xb[0], T14_h=xb[1], flat=xb[2], dT_min=xb[3], chi2=cb), info=info)
print('best', xb.round(3), cb, ndata)
prof_F = []
for F in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    x, c = fit(LC, fixF=F); prof_F.append([F, c - cb, x[1], x[0]])
    print(f'flat {F:.1f}: dchi2 {c-cb:.2f} T14 {x[1]:.2f} A {x[0]:.1f}', flush=True)
prof_T = []
for T in (1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9):
    x, c = fit(LC, fixT=T); prof_T.append([T, c - cb, x[2], x[0]])
    print(f'T14 {T:.1f}: dchi2 {c-cb:.2f} flat {x[2]:.2f} A {x[0]:.1f}', flush=True)
out['profile_flat'] = prof_F; out['profile_T14'] = prof_T
# jackknife over windows: refit, and profile delta-chi2 of a V-like shape (flat 0.16 = max for a dark body of
# k = 0.72, the G17.5-star case) and flat 0.30, each relative to the subsample's own best fit (fine grid in flat)
jk = []; jprof = []
FG = np.array([0.0, 0.15, 0.30, 0.45, 0.6, 0.75, 0.9, 1.0])
def prof_min(LCs):
    cs = [fit(LCs, fixF=F, fast=True)[1] for F in FG]
    return np.array(cs)
for i in range(len(LC)):
    sub = [L for j, L in enumerate(LC) if j != i]
    cs = prof_min(sub); ib = int(np.argmin(cs))
    jk.append(FG[ib])
    jprof.append(dict(dropped=f"S{LC[i]['sec']} E{LC[i]['epoch']}", best_flat=float(FG[ib]),
                      dchi2_flat0=float(cs[0] - cs[ib]), dchi2_flat0p15=float(cs[1] - cs[ib]), dchi2_flat0p30=float(cs[2] - cs[ib])))
    print(jprof[-1], flush=True)
jk = np.array(jk); nj = len(jk)
out['jackknife_flat_best'] = jk.tolist()
out['jackknife_flat_sd'] = float(np.sqrt((nj - 1) / nj * ((jk - jk.mean()) ** 2).sum()))
out['jackknife_profiles'] = jprof
csall = prof_min(LC); out['profile_flat_fine'] = [[float(F), float(c - csall.min())] for F, c in zip(FG, csall)]
print('fine profile', np.round(csall - csall.min(), 2).tolist())
np.savez_compressed(f'{ROOT}/out/shape_lc.npz', **{f'{L["sec"]}_{L["epoch"]}_{k}': L[k] for L in LC for k in ('t', 'f', 'e')})
out['per_window'] = [dict(sec=L['sec'], epoch=L['epoch'], tc=L['tc']) for L in LC]
jdump(out, f'{ROOT}/out/shape.json')
# figure
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ph = np.concatenate([(L['t'] - L['tc']) * 24 for L in LC]); ff = np.concatenate([L['f'] for L in LC])
o = np.argsort(ph); bins = np.arange(-10, 10.01, 0.2)
idx = np.digitize(ph, bins)
bm = np.array([ff[idx == i].mean() if (idx == i).sum() > 3 else np.nan for i in range(1, len(bins))])
ax[0].plot(ph, ff, '.', ms=1, color='0.7'); ax[0].plot(0.5 * (bins[1:] + bins[:-1]), bm, 'ko', ms=3)
tt = np.linspace(-10, 10, 2000) / 24
ax[0].plot(tt * 24, xb[0] * trapezoid(tt, xb[3] / 1440, 1e-6, t14=xb[1] / 24, flat=xb[2]), 'r-')
ax[0].set_xlim(-4, 4); ax[0].set_ylim(-20, 35); ax[0].set_xlabel('hours from mid-transit'); ax[0].set_ylabel('flux lost at M3V (e/s)')
pf = np.array(prof_F); ax[1].plot(pf[:, 0], pf[:, 1], 'o-'); ax[1].set_xlabel('T23/T14'); ax[1].set_ylabel('delta chi2')
fig.tight_layout(); fig.savefig(f'{ROOT}/figs/shape.png', dpi=90)
