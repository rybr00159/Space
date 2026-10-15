# Final referee f02: separation posterior sensitivity (own implementation, same ingredients as f01)
# Variants: IPD treatment (soft ramp / hard / empirical dI3-5 / empirical dI2-5 / logistic fit / none),
# VHS nominal vs degraded (+0.3 mag: real PSF worse than fitted PSF), PMa+ETV weights always applied.
import numpy as np, json
from scipy.optimize import minimize
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
rng = np.random.default_rng(99)
d_pc = 1000 / 3.3622
z = np.load(B + 'localise/v3/results/s05_baseline_map.npz'); g = z['g']; Mmap = z['M']
Lt_map = np.exp(-0.5 * (6.18 / 6.5417) * (Mmap - Mmap.min()))
GN, GE = np.meshgrid(g, g, indexing='ij'); R = np.hypot(GN, GE)
rb = np.arange(0, 3.0001, 0.05); rc = 0.5 * (rb[1:] + rb[:-1])
Lt_r = np.array([Lt_map[(R >= a) & (R < b)].mean() for a, b in zip(rb[:-1], rb[1:])])
h4 = json.load(open(B + 'localise/v3/results/h04_vhs_inner.json'))
S = sorted({r['sep'] for r in h4[0]['inj']}); D = sorted({r['dm'] for r in h4[0]['inj']})
Rec = np.zeros((len(S), len(D)))
for fr in h4:
    for r in fr['inj']:
        Rec[S.index(r['sep']), D.index(r['dm'])] += r['frac'] / len(h4)
S = np.array(S); D = np.array(D)
def miss(s, dk):
    col = np.array([np.interp(dk, D, Rec[i]) for i in range(len(S))])
    return 1 - np.interp(s, np.r_[0.3, S, 1.5], np.r_[0.0, col, 1.0 if dk <= 3.0 else col[-1]])
dKs = np.linspace(2.26, 2.98, 13); wK = np.exp(-0.5 * ((dKs - 2.62) / 0.18)**2); wK /= wK.sum()
Vm = {'nominal': sum(w * miss(rc, k) for k, w in zip(dKs, wK)),
      'degraded+0.3': sum(w * miss(rc, k + 0.3) for k, w in zip(dKs, wK))}
cal = json.load(open(B + 'localise/v3/results/v03_gaia_like.json'))
edges = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 1.6, 2.5]); mid = 0.5 * (edges[1:] + edges[:-1])
def kn(dmset):
    k = np.zeros(7); n = np.zeros(7)
    for c in cal['all_blend']['cells']:
        if tuple(c['dm']) in dmset:
            i = int(np.argmin(np.abs(edges[:-1] - c['sep'][0]))); k[i] += c['k']; n[i] += c['n']
    return k, n
k25, n25 = kn({(2, 3), (3, 4), (4, 5)})
# logistic-in-log(s) rising part fitted to bins up to 1.2" (beyond, companion resolved -> rate drops; use bins<=1.2 only)
def nll(p):
    x = np.log(mid[:5]); pr = 1 / (1 + np.exp(-(p[0] + p[1] * x))); pr = np.clip(pr, 1e-9, 1 - 1e-9)
    return -np.sum(k25[:5] * np.log(pr) + (n25[:5] - k25[:5]) * np.log(1 - pr))
fit = minimize(nll, [-3, 2], method='Nelder-Mead')
def ipd(kind, s):
    if kind == 'soft': return np.clip((s - 0.35) / 0.15, 0, 1) * np.clip((1.45 - s) / 0.25, 0, 1)
    if kind == 'hard': return ((s >= 0.48) & (s <= 1.2)).astype(float)
    if kind == 'none': return np.ones_like(s)
    if kind == 'logistic25':
        pr = 1 / (1 + np.exp(-(fit.x[0] + fit.x[1] * np.log(np.clip(s, 0.05, None)))))
        return np.where(s <= 1.2, pr, pr * 0.2)
    dmset = {(3, 4), (4, 5)} if kind == 'emp35' else {(2, 3), (3, 4), (4, 5)}
    k, n = kn(dmset); rate = (k + 0.5) / (n + 1)
    i = np.clip(np.searchsorted(edges, s, side='right') - 1, -1, 6)
    o = np.where(i >= 0, rate[np.clip(i, 0, 6)], cal['all_blend']['single']['L']); return o
N = 3_000_000; Mt = 1.69
logP = rng.normal(5.03, 2.28, N); P_yr = 10**logP / 365.25; a = (Mt * P_yr**2)**(1 / 3)
e_ = rng.uniform(0, 0.8, N); Mm = rng.uniform(0, 2 * np.pi, N); E = Mm.copy()
for _ in range(25):
    E -= (E - e_ * np.sin(E) - Mm) / (1 - e_ * np.cos(E))
r3 = a * (1 - e_ * np.cos(E)); nu = 2 * np.arctan2(np.sqrt(1 + e_) * np.sin(E / 2), np.sqrt(1 - e_) * np.cos(E / 2))
w = rng.uniform(0, 2 * np.pi, N); cosi = rng.uniform(-1, 1, N)
s = r3 * np.sqrt(np.cos(w + nu)**2 + (np.sin(w + nu) * cosi)**2) / d_pc
a_tab = np.array([1, 5, 10, 20, 32, 50, 100, 160, 1e4]); ex_tab = np.array([1.0, 0.98, 0.70, 0.55, 0.25, 0.05, 0.011, 0.001, 0.0])
a_etv = np.array([1, 2, 5, 10, 14, 20, 25, 30, 40, 50, 80, 160, 1e4]); ex_etv = np.array([.955, .98, .99, .94, .84, .66, .46, .30, .10, .05, .007, 0, 0])
W0 = (1 - np.interp(np.log(a), np.log(a_tab), ex_tab)) * (1 - np.interp(np.log(a), np.log(a_etv), ex_etv)) * np.interp(s, rc, Lt_r, right=0)
def wq(x, wt, qs=(0.025, 0.16, 0.5, 0.84, 0.975)):
    o = np.argsort(x); c = np.cumsum(wt[o]); c /= c[-1]; return [round(float(np.interp(q, c, x[o])), 3) for q in qs]
out = {'logistic_fit_dI2-5': list(map(float, fit.x)), 'k25': k25.tolist(), 'n25': n25.tolist()}
for vk, vm in Vm.items():
    Wv = W0 * np.interp(s, rc, vm, right=0)
    for kind in ('soft', 'hard', 'emp35', 'emp25', 'logistic25', 'none'):
        W = Wv * ipd(kind, s); m = W > 0
        out[f'{vk}|{kind}'] = dict(q_arcsec=wq(s[m], W[m]), q_au=[round(x * d_pc) for x in wq(s[m], W[m])],
                                   P_lt_0p3=round(float(W[s < 0.3].sum() / W.sum()), 3), P_gt_0p8=round(float(W[s > 0.8].sum() / W.sum()), 3),
                                   P_gt_1p0=round(float(W[s > 1.0].sum() / W.sum()), 3))
json.dump(out, open(B + 'final_referee/f02_sep.json', 'w'), indent=1)
for k, v in out.items():
    print(k, v)
