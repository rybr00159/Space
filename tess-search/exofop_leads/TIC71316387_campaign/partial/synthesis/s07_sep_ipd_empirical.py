# Synthesis s07: sensitivity of the separation posterior (and of the position factor in the odds)
# to how the Gaia IPD/BP-RP-blend signature is turned into a separation likelihood.
# s02 used an assumed soft ramp (0 below 0.35", 1 above 0.50").  Here the likelihood is taken
# directly from the localise speckle calibration (v03_gaia_like.json, 'all_blend' cells):
# match rate of target-like Gaia diagnostics among imaged pairs per separation bin, pooled over
# I-band contrast 3-5 mag (the pair's predicted dIc 3.7-4.1), Jeffreys estimate (k+0.5)/(n+1).
# Below 0.1" the single-star match rate (0/2650) is used.  PMa/RV/RUWE exclusion (astro
# pma_orbit_summary, 0.5-0.7 Msun) is applied as a weight (1 - excluded fraction) vs a.
import numpy as np, json
from scipy.interpolate import RegularGridInterpolator
rng = np.random.default_rng(7)
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d_pc = 1000 / 3.3622
z = np.load(B + 'localise/v3/results/s05_baseline_map.npz'); g = z['g']; M = z['M']
scale = 6.18 / 6.5417
tess_interp = RegularGridInterpolator((g, g), M, bounds_error=False, fill_value=M.max())
L_tess = lambda dE, dN: np.exp(-0.5 * scale * tess_interp(np.c_[dN, dE]))
h4 = json.load(open(B + 'localise/v3/results/h04_vhs_inner.json'))
seps = sorted({r['sep'] for r in h4[0]['inj']}); dms = sorted({r['dm'] for r in h4[0]['inj']})
R = np.zeros((len(seps), len(dms)))
for fr in h4:
    for r in fr['inj']:
        R[seps.index(r['sep']), dms.index(r['dm'])] += r['frac'] / len(h4)
seps = np.array(seps); dms = np.array(dms)
seps_x = np.r_[0.3, seps, 1.5, 3.0]
R_x = np.vstack([np.zeros(len(dms)), R, np.tile(np.where(dms <= 3.0, 1.0, R[-1]), (2, 1))])
rec = RegularGridInterpolator((seps_x, dms), R_x, bounds_error=False, fill_value=None)
def L_vhs(s, dK):
    s = np.clip(s, 0.3, 3.0); dK = np.clip(dK, 2.0, 4.0)
    return np.clip(1 - rec(np.c_[s, dK]), 0, 1)

cal = json.load(open(B + 'localise/v3/results/v03_gaia_like.json'))
def ipd_table(key, dmset):
    edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 1.6, 2.5]
    rates = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = n = 0
        for c in cal[key]['cells']:
            if abs(c['sep'][0] - lo) < 1e-6 and tuple(c['dm']) in dmset:
                k += c['k']; n += c['n']
        rates.append(((k + 0.5) / (n + 1), k, n))
    return edges, rates
single = cal['all_blend']['single']['L']
def make_Lipd(edges, rates):
    r = np.array([x[0] for x in rates])
    def f(s):
        s = np.asarray(s, float)
        out = np.full(s.shape, single)
        i = np.searchsorted(edges, s) - 1
        ok = (i >= 0) & (i < len(r))
        out[ok] = r[i[ok]]
        out[s >= edges[-1]] = r[-1]
        return out
    return f
soft = lambda s: np.clip((s - 0.35) / 0.15, 0, 1) * np.clip((1.45 - s) / 0.25, 0, 1)
variants = {}
e35, r35 = ipd_table('all_blend', {(3, 4), (4, 5)})
e25, r25 = ipd_table('all_blend', {(2, 3), (3, 4), (4, 5)})
variants['soft_ramp_s02'] = soft
variants['empirical_blend_dm3-5'] = make_Lipd(e35, r35)
variants['empirical_blend_dm2-5'] = make_Lipd(e25, r25)
variants['no_ipd'] = lambda s: np.ones_like(np.asarray(s, float))

# PMa/RV/RUWE exclusion vs a (astro pma_orbit_summary; 0.5-0.7 Msun companion)
a_tab = np.array([1, 5, 10, 20, 32, 50, 100, 160, 1e4]); ex_tab = np.array([1.0, 0.98, 0.70, 0.55, 0.25, 0.05, 0.011, 0.001, 0.0])
W_pma = lambda a: 1 - np.interp(np.log(a), np.log(a_tab), ex_tab)

N = 2_000_000
Mt = 1.03 + 0.66
logP = rng.normal(5.03, 2.28, N); P_yr = 10**logP / 365.25; a = (Mt * P_yr**2)**(1 / 3)
e = rng.uniform(0, 0.8, N); Mm = rng.uniform(0, 2 * np.pi, N); E = Mm.copy()
for _ in range(30):
    E = E - (E - e * np.sin(E) - Mm) / (1 - e * np.cos(E))
nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
r = a * (1 - e * np.cos(E)); inc = np.arccos(rng.uniform(-1, 1, N)); w = rng.uniform(0, 2 * np.pi, N); Om = rng.uniform(0, 2 * np.pi, N)
u = w + nu
X = r * (np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(inc))
Y = r * (np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(inc))
dN = X / d_pc; dE = Y / d_pc; s = np.hypot(dN, dE)
dK = rng.normal(2.62, 0.18, N)
sel = s < 3.0
base = np.zeros(N); base[sel] = L_tess(dE[sel], dN[sel]) * L_vhs(s[sel], dK[sel]) * W_pma(a[sel])
def wq(x, wt, qs=(0.025, 0.16, 0.5, 0.84, 0.975)):
    o = np.argsort(x); c = np.cumsum(wt[o]); c /= c[-1]
    return [round(float(np.interp(q, c, x[o])), 3) for q in qs]
# background: uniform in area, same position likelihood (no PMa weight: unrelated star)
xx = np.linspace(-3, 3, 601); XX, YY = np.meshgrid(xx, xx); dA = (xx[1] - xx[0])**2
SS = np.hypot(XX, YY).ravel()
LT_grid = L_tess(XX.ravel(), YY.ravel())
dKg = np.linspace(2.26, 2.98, 9); wK = np.exp(-0.5 * ((dKg - 2.62) / 0.18)**2); wK /= wK.sum()
LV_grid = sum(wk * L_vhs(SS, np.full(SS.shape, dk)) for dk, wk in zip(dKg, wK))
fsep_ref = 0.2716531478730537; circ = np.pi * 1.24**2
out = {'ipd_rates_dm3-5': [[lo, hi, round(x[0], 4), x[1], x[2]] for lo, hi, x in zip(e35[:-1], e35[1:], r35)],
       'ipd_rates_dm2-5': [[lo, hi, round(x[0], 4), x[1], x[2]] for lo, hi, x in zip(e25[:-1], e25[1:], r25)],
       'single_rate': single, 'variants': {}}
for name, f in variants.items():
    Wb = base * f(s)
    m = Wb > 0
    q = wq(s[m], Wb[m])
    pb = {k: round(float(Wb[(s >= lo) & (s < hi)].sum() / Wb.sum()), 3) for k, (lo, hi) in
          {'<0.1': (0, 0.1), '0.1-0.3': (0.1, 0.3), '0.3-0.45': (0.3, 0.45), '0.45-0.8': (0.45, 0.8), '0.8-1.0': (0.8, 1.0), '>1.0': (1.0, 9)}.items()}
    Lbg = float((LT_grid * LV_grid * f(SS)).sum() * dA)
    PF = (float(Wb.mean()) / fsep_ref) / (Lbg / circ)
    PA = np.degrees(np.arctan2(dE, dN)) % 360
    inN = (PA >= 300) | (PA <= 20)
    out['variants'][name] = dict(sep_q_2p5_16_50_84_97p5=q, proj_au=[round(x * d_pc) for x in q], mass_in_bins=pb,
                                 position_factor_vs_hard_window=round(PF, 3), PA_frac_300_20=round(float(Wb[inN].sum() / Wb.sum()), 3))
json.dump(out, open(B + 'synthesis/s07_sep_ipd_empirical.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

# ---------- outer-orbit observables per IPD variant ----------
n_ = 2 * np.pi / P_yr; h = n_ * a**2 * np.sqrt(1 - e**2); p_ = a * (1 - e**2)
rdot = h / p_ * e * np.sin(nu); rf = h / r
vx = rdot * (np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(inc)) + rf * (-np.cos(Om) * np.sin(u) - np.sin(Om) * np.cos(u) * np.cos(inc))
vy = rdot * (np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(inc)) + rf * (-np.sin(Om) * np.sin(u) + np.cos(Om) * np.cos(u) * np.cos(inc))
vz = rdot * np.sin(u) * np.sin(inc) + rf * np.cos(u) * np.sin(inc)
Zc = r * np.sin(u) * np.sin(inc)
GM = 4 * np.pi**2 * Mt; acc = GM / r**2; az = -acc * Zc / r
AU_YR2 = 1.495978707e11 / 3.15576e7**2; c_l = 2.99792458e8
obs = {'a_au': a, 'P_out_yr': P_yr, '|vrel_los|_kms': np.abs(vz) * 4.74047, 'mu_rel_masyr': np.hypot(vx, vy) / d_pc * 1000,
       'ETV_curv_min_yr2': np.abs(1.03 / Mt * az * AU_YR2) / (2 * c_l) * 3.15576e7**2 / 60,
       'RV_trend_G_ms_yr': np.abs(0.66 / Mt * az * AU_YR2) * 3.15576e7,
       'acc_sky_G_masyr2': 0.66 / Mt * acc * np.hypot(X, Y) / r / d_pc * 1000,
       'PA_rate_degyr': np.abs(np.degrees((X * vy - Y * vx) / (X**2 + Y**2)))}
for name in ('soft_ramp_s02', 'empirical_blend_dm3-5', 'empirical_blend_dm2-5'):
    Wb = base * variants[name](s); m = Wb > 0
    out['variants'][name]['outer_orbit_q'] = {k: wq(v[m], Wb[m]) for k, v in obs.items()}
    out['variants'][name]['P(|vrel_los|<3)'] = round(float(Wb[obs['|vrel_los|_kms'] < 3].sum() / Wb.sum()), 3)
# Gaia calibration: among pairs matching the IPD/RUWE signature without the blend criterion,
# how many also show the BP/RP blend flag, inside vs outside 0.5"
inside = [0, 0]; outside = [0, 0]
for cb, cn in zip(cal['all_blend']['cells'], cal['all_noblend']['cells']):
    assert cb['sep'] == cn['sep'] and cb['dm'] == cn['dm']
    tgt = inside if cb['sep'][1] <= 0.5 else (outside if cb['sep'][1] <= 1.2 else None)
    if tgt is not None:
        tgt[0] += cb['k']; tgt[1] += cn['k']
from scipy.stats import fisher_exact
out['blend_given_ipd_match'] = dict(inside_0p5=inside, sep_0p5_1p2=outside,
                                    fisher_p_one_sided=float(fisher_exact([[inside[0], inside[1] - inside[0]], [outside[0], outside[1] - outside[0]]], alternative='less')[1]))
# match rates dm3-5: <0.5" vs 0.5-1.2"
k_in = sum(x[1] for x in r35[:2]); n_in = sum(x[2] for x in r35[:2]); k_out = sum(x[1] for x in r35[2:5]); n_out = sum(x[2] for x in r35[2:5])
out['match_rate_dm3-5'] = dict(inside_0p5=[k_in, n_in], sep_0p5_1p2=[k_out, n_out],
                               fisher_p_one_sided=float(fisher_exact([[k_in, n_in - k_in], [k_out, n_out - k_out]], alternative='less')[1]))
json.dump(out, open(B + 'synthesis/s07_sep_ipd_empirical.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('blend_given_ipd_match', 'match_rate_dm3-5')}, indent=1))
for name in ('soft_ramp_s02', 'empirical_blend_dm3-5', 'empirical_blend_dm2-5'):
    print(name, json.dumps(out['variants'][name]['outer_orbit_q']), out['variants'][name]['P(|vrel_los|<3)'])
