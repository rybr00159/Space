"""Out-of-eclipse variability at P_EB and P_EB/2.
Per sector and product: eclipses masked (|dt| < 1.0 T14), outliers clipped; joint weighted least squares of
  flux = spline(t; knots every KS d, separately per spacecraft orbit) * [1 + a1 cos(2 pi phi) + b1 sin(2 pi phi)
                                                                       + a2 cos(4 pi phi) + b2 sin(4 pi phi)]
(linearised: flux ~ S(t) + harmonics), phi = (t - T0)/P_EB.  Conventions: ellipsoidal = -a2 (minima at the
eclipses), reflection = -a1 (maximum at phase 0.5), beaming = b1.
Empirical null: the same fit at 400 fake periods in [1.0, 1.75] d (avoiding +-3% of P_EB, 2/3, 3/2 P_EB) gives the
noise distribution of each amplitude; injection-recovery of ellipsoidal/reflection signals into the raw flux
(before detrending) calibrates suppression by the spline and gives 95% upper limits.
Products: SPOC PDCSAP and SAP (S1, S28, S68), own TESScut (all), QLP SAP (S95, S105)."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import data as D
from scipy.interpolate import BSpline
V = D.V
P = D.P0; T0 = D.T00
KS = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
rng = np.random.default_rng(42)

def segments(t, gap=1.0):
    br = np.where(np.diff(t) > gap)[0]
    return np.split(np.arange(t.size), br + 1)

def spline_basis(t, ks):
    cols = []
    for seg in segments(t):
        ts = t[seg]; span = ts[-1] - ts[0]
        nint = max(1, int(np.round(span / ks)))
        inner = np.linspace(ts[0], ts[-1], nint + 1)
        kn = np.r_[[inner[0]] * 3, inner, [inner[-1]] * 3]
        nb = len(kn) - 4
        for j in range(nb):
            c = np.zeros(nb); c[j] = 1
            col = np.zeros(t.size)
            col[seg] = BSpline(kn, c, 3, extrapolate=False)(ts)
            cols.append(np.nan_to_num(col))
    return np.array(cols).T

def harm(t, Pf):
    ph = 2 * np.pi * (t - T0) / Pf
    return np.vstack([np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)]).T

def fit(t, f, Pf, B=None):
    if B is None: B = spline_basis(t, KS)
    H = harm(t, Pf)
    A = np.hstack([B, H])
    keep = np.ones(t.size, bool)
    for _ in range(3):
        c, *_ = np.linalg.lstsq(A[keep], f[keep], rcond=None)
        r = f - A @ c; s = 1.4826 * np.median(np.abs(r[keep]))
        keep = np.abs(r) < 5 * s
    trend = B @ c[:B.shape[1]]
    amp = c[B.shape[1]:] / np.median(trend)
    return amp, r, keep

def prep(lc):
    t, f = lc['t'], lc['f']
    d1, d2 = D.phase_dist(t)
    m = (d1 > 1.0 * D.T14) & (d2 > 1.0 * D.T14) & np.isfinite(f)
    return t[m], f[m]

def analyse(lc, label, ninj=60, nnull=400):
    t, f = prep(lc)
    B = spline_basis(t, KS)
    amp, r, keep = fit(t, f, P, B)
    # binned-noise for info
    out = dict(label=label, n=int(t.size), rms_ppm=float(np.std(r[keep]) * 1e6),
               a1=float(amp[0]), b1=float(amp[1]), a2=float(amp[2]), b2=float(amp[3]),
               A1=float(np.hypot(amp[0], amp[1])), A2=float(np.hypot(amp[2], amp[3])))
    # empirical null at fake periods
    fakes = []
    while len(fakes) < nnull:
        Pf = rng.uniform(1.0, 1.75)
        if min(abs(Pf / P - x) for x in (1.0, 2 / 3, 1.5, 4 / 3, 0.75)) < 0.03: continue
        fakes.append(Pf)
    nul = np.array([fit(t, f, Pf, B)[0] for Pf in fakes])
    for j, nm in enumerate(['a1', 'b1', 'a2', 'b2']):
        out[f'null_sd_{nm}'] = float(1.4826 * np.median(np.abs(nul[:, j] - np.median(nul[:, j]))))
        out[f'null_med_{nm}'] = float(np.median(nul[:, j]))
    out['null_A2_95'] = float(np.percentile(np.hypot(nul[:, 2], nul[:, 3]), 95))
    out['null_A1_95'] = float(np.percentile(np.hypot(nul[:, 0], nul[:, 1]), 95))
    # injection-recovery of ellipsoidal (-A cos 4 pi phi) and reflection (-A cos 2 pi phi) into raw flux
    ph = 2 * np.pi * (t - T0) / P
    inj = []
    for A in (50e-6, 100e-6, 200e-6):
        for _ in range(ninj // 3):
            ff = f * (1 - A * np.cos(2 * ph) - A * np.cos(ph))
            a, *_ = fit(t, ff, P, B)
            inj.append((A, -(a[2] - amp[2]) / A, -(a[0] - amp[0]) / A))
    inj = np.array(inj)
    out['inj_recovery_ell'] = float(np.median(inj[:, 1])); out['inj_recovery_ref'] = float(np.median(inj[:, 2]))
    return out

if __name__ == '__main__':
    res = []
    for sec in (1, 28, 68):
        res.append(analyse(D.clip(D.load_spoc(sec)), f'SPOC-PDC{sec}'))
        res.append(analyse(D.clip(D.load_spoc(sec, col='SAP_FLUX')), f'SPOC-SAP{sec}'))
    for sec in (1, 28, 68, 95, 105, 106):
        res.append(analyse(D.clip(D.load_tc(sec)), f'TC{sec}'))
    for sec in (95, 105):
        res.append(analyse(D.clip(D.load_qlp(sec)), f'QLP{sec}'))
    for o in res:
        print(f"{o['label']:12s} rms {o['rms_ppm']:6.0f}  ell=-a2 {-o['a2']*1e6:+6.1f} +- {o['null_sd_a2']*1e6:5.1f} (b2 {o['b2']*1e6:+6.1f})  "
              f"refl=-a1 {-o['a1']*1e6:+6.1f} +- {o['null_sd_a1']*1e6:5.1f}  beam b1 {o['b1']*1e6:+6.1f} +- {o['null_sd_b1']*1e6:5.1f}  "
              f"A1 {o['A1']*1e6:5.1f} (null95 {o['null_A1_95']*1e6:5.1f})  A2 {o['A2']*1e6:5.1f} (null95 {o['null_A2_95']*1e6:5.1f})  inj-rec ell {o['inj_recovery_ell']:.2f} ref {o['inj_recovery_ref']:.2f}", flush=True)
    # combined (inverse-variance with null sd) over the default product set
    comb = {}
    for setname, labels in (('default', ['SPOC-PDC1', 'SPOC-PDC28', 'SPOC-PDC68', 'TC95', 'TC105', 'TC106']),
                            ('spoc_sap+tc', ['SPOC-SAP1', 'SPOC-SAP28', 'SPOC-SAP68', 'TC95', 'TC105', 'TC106']),
                            ('tc_all', ['TC1', 'TC28', 'TC68', 'TC95', 'TC105', 'TC106'])):
        sel = [o for o in res if o['label'] in labels]
        c = {}
        for nm in ('a1', 'b1', 'a2', 'b2'):
            w = np.array([1 / o[f'null_sd_{nm}'] ** 2 for o in sel]); v = np.array([o[nm] / o[f'inj_recovery_{"ell" if nm[1]=="2" else "ref"}'] if nm[0]=='a' else o[nm] for o in sel])
            c[nm] = float((w * v).sum() / w.sum()); c[nm + '_err'] = float(1 / np.sqrt(w.sum()))
            c[nm + '_chi2'] = float(((v - c[nm]) ** 2 * w).sum()); c[nm + '_n'] = len(sel)
        comb[setname] = c
        print(setname, ' '.join(f"{k}={v*1e6:+.1f}" if not k.endswith(('chi2', '_n')) else f"{k}={v:.1f}" for k, v in c.items()))
    json.dump(dict(knot_spacing_d=KS, per_product=res, combined=comb), open(f'{V}/results/ooe_ks{KS:g}.json', 'w'), indent=1)
