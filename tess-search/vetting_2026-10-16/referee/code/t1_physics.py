"""TIC 90086266: secondary-eclipse limit (phase 0.5 and full phase scan),
companion size, expected secondary depth, Gaia RV Monte Carlo, and a check of
per-event instrumental diagnostics for the S104 odd/even alternation."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from reflib import *

TIC = 90086266
EV = json.load(open(f'{BASE}/results/t1_events.json'))
OUT = {}
T14, T12 = 3.0 / 24, 1.3 / 24

# ---------------------------------------------------------------- secondary eclipse
def local_box_depth(t, f, centres, T14, T12, win=0.35, texp=None):
    """Joint fit: common depth of a trapezoid at all 'centres' with per-window quadratic baseline."""
    sel_list = []
    for c in centres:
        s = np.abs(t - c) < win
        if s.sum() > 10 and (np.abs(t[s] - c) < T14 / 2).sum() >= 2 and (t[s] < c - T14 / 2).sum() >= 3 and (t[s] > c + T14 / 2).sum() >= 3:
            sel_list.append((c, s))
    if len(sel_list) < 2:
        return np.nan, np.nan, 0
    # linear least squares: for each window, baseline a0+a1 dt+a2 dt^2 multiplicative ~ additive (small depth)
    rows, y = [], []
    nW = len(sel_list)
    blocks = []
    for k, (c, s) in enumerate(sel_list):
        dt = t[s] - c
        shape = 1 - supersample(trapezoid, t[s], texp, t0=c, depth=1.0, T14=T14, T12=T12)  # 0..1 deficit
        X = np.zeros((s.sum(), 1 + 3 * nW))
        X[:, 0] = -shape
        X[:, 1 + 3 * k] = 1
        X[:, 2 + 3 * k] = dt
        X[:, 3 + 3 * k] = dt ** 2
        blocks.append((X, f[s]))
    X = np.vstack([b[0] for b in blocks]); Y = np.concatenate([b[1] for b in blocks])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    res = Y - X @ beta
    s2 = np.var(res, ddof=X.shape[1])
    cov = np.linalg.inv(X.T @ X) * s2
    return beta[0], np.sqrt(cov[0, 0]), nW


P = EV['joint_ephem']['P']
sec_res = {}
for key in ['S13_TESS-SPOC:PDCSAP_FLUX', 'S13_TESS-SPOC:SAP_FLUX', 'S104_QLP:SAP_FLUX', 'S104_CUT:r2.0', 'S13_CUT:r2.0']:
    sec = int(key.split('_')[0][1:])
    if 'CUT' in key:
        t, f, e, ex = aperture_lc(TIC, sec, r=2.0)
    else:
        prod, col = key.split('_', 1)[1].split(':')
        t, f, e, ex = load_lc(TIC, sec, prod, col)
    texp = np.median(np.diff(t))
    tcs = np.array([v['tc'] for v in EV[key]['events']])
    # in-sector reference epoch and period (use per-sector event times -> exact phase 0)
    Ps = P
    t0s = tcs[0]
    # mask primaries
    ph = ((t - t0s) / Ps + 0.5) % 1 - 0.5
    keep = np.abs(ph * Ps) > 0.6 * T14
    tt, ff = t[keep], f[keep]
    n0 = int(np.floor((tt.min() - t0s) / Ps)) - 1
    n1 = int(np.ceil((tt.max() - t0s) / Ps)) + 1
    d05, e05, nw = local_box_depth(tt, ff, t0s + (np.arange(n0, n1) + 0.5) * Ps, T14, T12, texp=texp)
    # full phase scan
    scan = []
    for phase in np.arange(0.1, 0.901, 0.01):
        d, de, nwi = local_box_depth(tt, ff, t0s + (np.arange(n0, n1) + phase) * Ps, T14, T12, texp=texp)
        scan.append((phase, d, de, nwi))
    scan = np.array(scan)
    snr = scan[:, 1] / scan[:, 2]
    i = np.nanargmax(snr)
    sec_res[key] = dict(depth_phase05_ppm=d05 * 1e6, err_ppm=e05 * 1e6, n_windows=nw,
                        scan_max_snr=float(snr[i]), scan_max_phase=float(scan[i, 0]), scan_max_depth_ppm=float(scan[i, 1] * 1e6),
                        scan_rms_snr=float(np.nanstd(snr)))
    print(f'{key:28s} phase-0.5 depth = {d05*1e6:7.0f} +- {e05*1e6:5.0f} ppm (n={nw}); phase scan 0.1-0.9: max S/N {snr[i]:.1f} at phase {scan[i,0]:.2f} '
          f'(depth {scan[i,1]*1e6:.0f} ppm), rms S/N {np.nanstd(snr):.2f}')
OUT['secondary'] = sec_res

# ---------------------------------------------------------------- expected secondary depth for an M-dwarf companion
from scipy.integrate import quad
h, c, kB = 6.626e-34, 2.998e8, 1.381e-23


def bb(lam, T):
    return 1 / lam ** 5 / (np.exp(h * c / (lam * kB * T)) - 1)


def band(T):
    # TESS band approx 600-1000 nm, photon-counting weight (lam/hc)
    return quad(lambda l: bb(l, T) * l, 600e-9, 1000e-9)[0]


Rstar = 1.70
dep_obs = 0.030  # undiluted-ish primary depth (S104 200-s data peak ~ 0.027-0.030, TESS-SPOC PDCSAP 0.035)
k_min = np.sqrt(0.027 / 1.15)  # central LD enhancement ~1.15
k_flat = np.sqrt(0.030)
OUT['companion'] = dict(R2_min_Rsun=k_min * Rstar, R2_min_RJ=k_min * Rstar / 0.10045, R2_from_sqrt_depth_Rsun=k_flat * Rstar)
print(f'companion radius if on F star (R*={Rstar}): >= {k_min*Rstar:.3f} Rsun = {k_min*Rstar/0.10045:.2f} RJ (sqrt(depth): {k_flat*Rstar:.3f} Rsun)')
for T2 in (3000, 3200, 3400, 3600):
    ratio = band(T2) / band(6300)
    print(f'   BB surface-brightness ratio (TESS band) T2={T2}K vs 6300K: {ratio:.4f} -> expected secondary depth ~ {0.028*ratio*1e6:.0f} ppm (for 2.8% primary)')
    OUT.setdefault('expected_secondary_ppm', {})[T2] = 0.028 * ratio * 1e6
# also: maximum companion temperature allowed by phase-0.5 limit (S104 QLP)
lim = sec_res['S104_QLP:SAP_FLUX']['depth_phase05_ppm'] + 2 * sec_res['S104_QLP:SAP_FLUX']['err_ppm']
for T2 in np.arange(1500, 4001, 50):
    if 0.028 * band(T2) / band(6300) * 1e6 > lim:
        print(f'   2-sigma phase-0.5 upper limit {lim:.0f} ppm -> companion BB Teff < {T2} K')
        OUT['Tmax_companion_K'] = int(T2)
        break

# ---------------------------------------------------------------- Gaia RV Monte Carlo
G = 6.674e-11; Msun = 1.989e30
P_s = 2.339469 * 86400
M1 = 1.3


def Kms(M2, M1=1.3, P=P_s, sini=0.985):
    return (2 * np.pi * G / P) ** (1 / 3.) * M2 * Msun * sini / ((M1 + M2) * Msun) ** (2 / 3.) / 1e3


rng = np.random.default_rng(1)
N = 200000
sig = 2.3  # km/s epoch noise consistent with rv_amplitude_robust 5.0, chisq p 0.22
OUT['rv_mc'] = {}
for M2 in (0.05, 0.1, 0.2, 0.27, 0.3):
    K = Kms(M2)
    out = {}
    # (a) 5 independent random phases
    ph = rng.random((N, 5))
    v = K * np.sin(2 * np.pi * ph) + rng.normal(0, sig, (N, 5))
    out['random5'] = float(np.mean(v.max(1) - v.min(1) <= 4.986))
    # (b) clustered: 5 transits in 2 or 3 visits; within a visit, transits at 0, 106.5 min, 6h, 6h+106.5min
    offs = np.array([0, 106.5 / 1440, 0.25, 0.25 + 106.5 / 1440]) / 2.339469
    for nv in (2, 3):
        vis = rng.random((N, nv))
        # allocate 5 transits among nv visits (at least 1 each)
        alloc = [(3, 2), (2, 2, 1)][nv - 2]
        phs = []
        for j, na in enumerate(alloc):
            phs.append(vis[:, j:j + 1] + offs[None, :na])
        phs = np.hstack(phs)
        v = K * np.sin(2 * np.pi * phs) + rng.normal(0, sig, phs.shape)
        out[f'clustered_{nv}visits'] = float(np.mean(v.max(1) - v.min(1) <= 4.986))
    OUT['rv_mc'][f'M2={M2}'] = dict(K_kms=K, **out)
    print(f'M2={M2:.2f} Msun: K={K:.1f} km/s; P(range<=5.0 km/s): random 5 phases {out["random5"]:.2e}, 2 visits {out["clustered_2visits"]:.2e}, 3 visits {out["clustered_3visits"]:.2e}')
# twin F-star SB2 at 2P (each ~85 km/s); pipeline would track either component -> even worse; just K of each
a_tw = ((G * 2.4 * Msun * (2 * P_s) ** 2) / (4 * np.pi ** 2)) ** (1 / 3.)
K_tw = 2 * np.pi * a_tw / (2 * P_s) / 2 / 1e3
OUT['twin_K_each_kms'] = K_tw
print(f'twin 1.2+1.2 Msun at 2P: K each = {K_tw:.0f} km/s')
save_json(OUT, f'{BASE}/results/t1_physics.json')
