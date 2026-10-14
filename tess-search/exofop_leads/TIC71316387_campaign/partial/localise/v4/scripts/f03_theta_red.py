"""v4: field-angle dependence of the TESS lateral-colour term and the red-end colour curve (whole CCD, 2.6 deg).

Per FFI: residual = poly3(x, y) per component + m*(G-13) + sum_b (C-C0) * k_b * u_rad * [theta in bin b] + k_t (C-C0) u_tan.
-> k_rad(theta) per theta bin; fit k(theta) = k_t * (theta/theta_target)^alpha jointly over FFIs (alpha common).
Red-end curve: radial residual (after poly + mag terms) of every star rescaled to the target's field angle with the
fitted alpha, binned in colour, for C up to 4 (stars from the whole CCD), compared with the linear model.
Chi2r cut relaxed to < 100 (bright stars have PRF-mismatch chi2), robust 4-sigma clipping.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, json
from scipy.optimize import minimize_scalar
from v4lib import V4, to_json

C0 = 0.786
TB = [(0, 5), (5, 8), (8, 10), (10, 12), (12, 14), (14, 18)]


def load(f):
    z = np.load(f); a = z['stars']
    ok = (a['chi2r'] < 100) & (a['ruwe'] < 1.4) & (a['flux'] > 0) & (a['ex'] < 0.25) & (a['ey'] < 0.25)
    a = a[ok]; J = z['J']
    d = J @ np.vstack([a['xf'] - a['xw'], a['yf'] - a['yw']])
    e = np.sqrt((J ** 2) @ np.vstack([a['ex'] ** 2, a['ey'] ** 2]))
    pa = np.radians(a['pa_bs']); ur = np.vstack([np.sin(pa), np.cos(pa)]); ut = np.vstack([ur[1], -ur[0]])
    return dict(a=a, d=d, e=e, ur=ur, ut=ut, sec=int(z['sec']), xt=float(z['xt']), yt=float(z['yt']), tht=float(z['theta_t']), file=os.path.basename(f))


def basecols(F):
    a = F['a']; N = len(a)
    x = (a['xw'] - 1024) / 1000.0; y = (a['yw'] - 1024) / 1000.0
    poly = [np.ones(N), x, y, x * x, x * y, y * y, x ** 3, x * x * y, x * y * y, y ** 3]
    Z = np.zeros(N)
    cols = [np.r_[p, Z] for p in poly] + [np.r_[Z, p] for p in poly] + [np.r_[a['G'] - 13, Z], np.r_[Z, a['G'] - 13]]
    return cols


def robust(A, b, e, use):
    n = len(b) // 2; ok = use.copy()
    for _ in range(15):
        w = 1 / e[ok]
        p, *_ = np.linalg.lstsq(A[ok] * w[:, None], b[ok] * w, rcond=None)
        r = (b - A @ p) / e; s = 1.4826 * np.median(np.abs(r[ok]))
        g = (np.abs(r) < 4 * s); g = g[:n] & g[n:]
        new = use & np.r_[g, g]
        if np.array_equal(new, ok):
            break
        ok = new
    return p, ok


frames = [load(f) for f in sorted(glob.glob(f'{V4}/results/f01_tess*.npz'))]
kb = {}
pts = []
for F in frames:
    a = F['a']; N = len(a)
    dc = a['C'] - C0
    cols = basecols(F)
    nb = 0
    for lo, hi in TB:
        m = (a['theta'] >= lo) & (a['theta'] < hi)
        cols.append(np.r_[dc * m * F['ur'][0], dc * m * F['ur'][1]]); nb += 1
    cols.append(np.r_[dc * F['ut'][0], dc * F['ut'][1]])
    A = np.column_stack(cols); b = np.r_[F['d'][0], F['d'][1]]; e = np.sqrt(np.r_[F['e'][0], F['e'][1]] ** 2 + 0.1 ** 2)
    use = np.r_[a['G'] < 15.5, a['G'] < 15.5]
    p, ok = robust(A, b, e, use)
    # bootstrap errors
    rng = np.random.default_rng(1); bs = []
    for _ in range(150):
        i = rng.integers(0, N, N); ii = np.r_[i, i + N]
        bs.append(robust(A[ii], b[ii], e[ii], use[ii])[0][-nb - 1:-1])
    bs = np.array(bs)
    kb[F['file']] = dict(sec=F['sec'], tht=F['tht'], k=p[-nb - 1:-1].tolist(), ek=bs.std(0).tolist(),
                         n=[int(((a['theta'] >= lo) & (a['theta'] < hi) & ok[:N]).sum()) for lo, hi in TB],
                         thmed=[float(np.median(a['theta'][(a['theta'] >= lo) & (a['theta'] < hi)])) if ((a['theta'] >= lo) & (a['theta'] < hi)).any() else np.nan for lo, hi in TB])
    # residual after base (poly + mag) terms only, using the base part of p
    nbase = len(basecols(F))
    rr = b - A[:, :nbase] @ p[:nbase]
    rad = rr[:N] * F['ur'][0] + rr[N:] * F['ur'][1]; tan = rr[:N] * F['ut'][0] + rr[N:] * F['ut'][1]
    pts.append(dict(C=a['C'], G=a['G'], th=a['theta'], tht=np.full(N, F['tht']), rad=rad, tan=tan, ok=ok[:N] | (a['G'] >= 15.5), sec=np.full(N, F['sec']),
                    sid=a['sid'], plx=a['plx'], eplx=a['eplx'], e=np.hypot(F['e'][0], F['e'][1]) / np.sqrt(2)))
    print(F['file'][4:30], 'k_rad by theta bin:', ' '.join(f"[{lo}-{hi}] {k:+.2f}+-{ek:.2f} (n={n}, th~{t:.1f})" for (lo, hi), k, ek, n, t in
                                                   zip(TB, kb[F['file']]['k'], kb[F['file']]['ek'], kb[F['file']]['n'], kb[F['file']]['thmed'])), flush=True)
# joint power-law fit k(theta) = K_f * (theta/theta_t)^alpha, K_f per FFI free (profile over alpha)
def chi_alpha(al, ret=False):
    tot = 0; Ks = {}
    for f, v in kb.items():
        th = np.array(v['thmed']); k = np.array(v['k']); ek = np.array(v['ek']); n = np.array(v['n'])
        m = np.isfinite(th) & (n >= 15) & (ek > 0)
        x = (th[m] / v['tht']) ** al
        K = (x * k[m] / ek[m] ** 2).sum() / (x * x / ek[m] ** 2).sum()
        tot += (((k[m] - K * x) / ek[m]) ** 2).sum(); Ks[f] = K
    return (tot, Ks) if ret else tot
grid = np.linspace(-1.0, 2.0, 61)
cg = np.array([chi_alpha(a) for a in grid])
al = float(grid[np.argmin(cg)])
c2min, Ks = chi_alpha(al, True)
in1 = grid[cg <= c2min + 1]
print(f'\npower law k ~ theta^alpha: alpha = {al:.2f} (dchi2<1: {in1.min():.2f}..{in1.max():.2f}); chi2 {c2min:.1f}; chi2(alpha=0) {chi_alpha(0.0):.1f}, chi2(alpha=1) {chi_alpha(1.0):.1f}')
print('  K at target field angle per FFI:', {f[4:30]: round(k, 3) for f, k in Ks.items()})
# red-end curve, rescaled to the target field angle
P_ = {k: np.concatenate([p[k] for p in pts]) for k in pts[0]}
sc = (P_['th'] / P_['tht']) ** al
radT = P_['rad'] / np.where(np.abs(P_['C'] - C0) > 0, sc, 1)
key = P_['sid'].astype(np.int64) * 1000 + P_['sec']
kmean = float(np.mean(list(Ks.values())))
bins = [(0.7, 0.85), (0.85, 1.0), (1.0, 1.2), (1.2, 1.45), (1.45, 1.75), (1.75, 2.05), (2.05, 2.35), (2.35, 2.65), (2.65, 3.0), (3.0, 3.5), (3.5, 5)]
curve = []
print(f'\n  red-end radial curve at the target field angle (arcsec; negative = away from boresight); linear k = {kmean:+.3f}"/mag:')
for lo, hi in bins:
    m = P_['ok'] & (P_['C'] >= lo) & (P_['C'] < hi)
    if m.sum() < 1:
        continue
    u, inv = np.unique(key[m], return_inverse=True)
    r = np.bincount(inv, radT[m]) / np.bincount(inv); c = np.bincount(inv, P_['C'][m]) / np.bincount(inv)
    gm = np.bincount(inv, P_['G'][m]) / np.bincount(inv)
    sids = np.unique(P_['sid'][m])
    med = float(np.median(r)); err = float(1.2533 * 1.4826 * np.median(np.abs(r - med)) / np.sqrt(len(u))) if len(u) > 2 else float(np.nan)
    curve.append(dict(lo=lo, hi=hi, n_meas=int(len(u)), n_stars=int(len(sids)), C=float(np.median(c)), rad=med, err=err, lin=kmean * (float(np.median(c)) - C0),
                      G_med=float(np.median(gm))))
    print(f"   C {lo:.2f}-{hi:.2f}: {len(u):3d} star-sectors ({len(sids):3d} stars, G~{np.median(gm):.1f}) medC {np.median(c):.2f} rad {med:+.2f} +- {err:.2f}  (linear {kmean * (np.median(c) - C0):+.2f})")
red = P_['ok'] & (P_['C'] > 2.0)
print('\n  individual red stars (C > 2), mean over FFIs of the same sector, rescaled to target theta:')
reds = []
for s in np.unique(P_['sid'][red]):
    m = red & (P_['sid'] == s)
    rr = radT[m]; e_ = P_['e'][m] / sc[m]
    reds.append(dict(sid=int(s), C=float(P_['C'][m][0]), G=float(P_['G'][m][0]), n=int(m.sum()), rad=float(np.mean(rr)), err=float(np.std(rr) / np.sqrt(m.sum()) if m.sum() > 2 else np.mean(e_)),
                     theta=float(np.mean(P_['th'][m])), plx=float(P_['plx'][m][0])))
reds.sort(key=lambda r: r['C'])
for r in reds:
    print(f"   {r['sid']} C={r['C']:.2f} G={r['G']:.1f} plx={r['plx']:.1f} n={r['n']} theta~{r['theta']:.1f}  rad={r['rad']:+.2f} +- {r['err']:.2f}  (lin {kmean * (r['C'] - C0):+.2f})")
to_json(dict(theta_bins=TB, per_ffi=kb, alpha=al, alpha_1sig=[float(in1.min()), float(in1.max())], chi2_alpha=dict(best=c2min, a0=chi_alpha(0.0), a1=chi_alpha(1.0)),
             K_target=Ks, k_mean=kmean, curve=curve, reds=reds), f'{V4}/results/f03_theta_red.json')
# per-sector (camera) check of the red-end response, relative to the linear model at each star's colour
print('\n  per-sector red-end response: median (rad / [k_mean*(C-C0)]) for 1.3<C<2.0 and C>2.0 (rescaled to target theta):')
persec_red = {}
for s in sorted(set(P_['sec'])):
    for lo, hi in ((1.3, 2.0), (2.0, 3.5)):
        m = P_['ok'] & (P_['sec'] == s) & (P_['C'] >= lo) & (P_['C'] < hi)
        if m.sum() < 2:
            continue
        ratio = radT[m] / (kmean * (P_['C'][m] - C0))
        persec_red[f'{s}_{lo}-{hi}'] = dict(n=int(m.sum()), median_ratio=float(np.median(ratio)), err=float(1.2533 * 1.4826 * np.median(np.abs(ratio - np.median(ratio))) / np.sqrt(m.sum())))
        print(f"   S{s:3d} C {lo}-{hi}: n={m.sum():3d} response ratio {np.median(ratio):.2f} +- {persec_red[f'{s}_{lo}-{hi}']['err']:.2f}")
json.dump(persec_red, open(f'{V4}/results/f03_persec_red.json', 'w'), indent=1)
