"""Detectability of an outer companion (mass M3) from the 2018-2026 eclipse timings via the light-travel-time effect.
For each outer semi-major axis a_out: Monte Carlo over isotropic inclination, random orbital phase, eccentricity
(circular, or thermal f(e)=2e with random omega).  Each configuration's LTTE delay is added to the best linear ephemeris
plus a bootstrap noise realisation (residuals resampled within sector, scaled by empirical per-event errors).
Detection = Delta-chi2(quadratic vs linear) > 99.7% null quantile OR max Delta-chi2(LTTE sinusoid grid vs linear) > its
99.7% null quantile (both quantiles from the same bootstrap with no signal).  Excluded fraction = detected fraction.
Usage: etv_companion.py <tag> <M3> <M_EB>"""
import numpy as np, json, sys
HERE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel'
tag = sys.argv[1]; M3 = float(sys.argv[2]); MEB = float(sys.argv[3])
D = json.load(open(f'{HERE}/results/etvfit_{tag}.json'))
E = np.array(D['E']); oc = np.array(D['oc']); es = np.array(D['es']); sec = np.array(D['sec']); ph = np.array(D['ph']); t = np.array(D['t'])
tyr = (t - 2900.0) / 365.25
n = len(oc); W = 1 / es
def proj(A):
    Aw = A * W[:, None]; return Aw @ np.linalg.pinv(Aw)     # hat matrix in weighted space
Al = np.array([np.ones(n), E / 1000, (ph == 0.5).astype(float)]).T
Aq = np.column_stack([Al, tyr ** 2])
Hl = proj(Al); Hq = proj(Aq)
Pg = np.exp(np.linspace(np.log(0.15), np.log(60), 120))
HL = [proj(np.column_stack([Al, np.sin(2 * np.pi * tyr / P), np.cos(2 * np.pi * tyr / P)])) for P in Pg]
def stats(Y):
    """Y: n x m matrix of O-C realisations (min). returns dchi2_quad, max dchi2_ltte"""
    Yw = Y * W[:, None]
    rl = Yw - Hl @ Yw; cl = np.sum(rl ** 2, 0)
    rq = Yw - Hq @ Yw; cq = np.sum(rq ** 2, 0)
    best = np.zeros(Y.shape[1])
    for H in HL:
        r = Yw - H @ Yw; best = np.maximum(best, cl - np.sum(r ** 2, 0))
    return cl - cq, best
cl_fit = np.linalg.lstsq(Al * W[:, None], oc * W, rcond=None)[0]
base = Al @ cl_fit; resid_n = (oc - base) / es
rng = np.random.default_rng(5)
def noise(m):
    Z = np.empty((n, m))
    for s in np.unique(sec):
        idx = np.where(sec == s)[0]; Z[idx] = rng.choice(resid_n[idx], (idx.size, m), replace=True)
    return Z * es[:, None]
# null thresholds
N0 = 4000
dq0, dl0 = stats(base[:, None] + noise(N0))
thr_q = float(np.percentile(dq0, 99.7)); thr_l = float(np.percentile(dl0, 99.7))
dq_obs, dl_obs = stats(oc[:, None])
print('null 99.7%% thresholds: quad %.2f  ltte-max %.2f ; observed quad %.2f ltte-max %.2f' % (thr_q, thr_l, dq_obs[0], dl_obs[0]))
G = 6.674e-11; MSUN = 1.989e30; AU = 1.496e11; C = 2.998e8; YR = 3.15576e7
def kepE(M, e):
    Ea = M + e * np.sin(M)
    for _ in range(50): Ea = Ea - (Ea - e * np.sin(Ea) - M) / (1 - e * np.cos(Ea))
    return Ea
out = dict(M3=M3, MEB=MEB, thr_quad=thr_q, thr_ltte=thr_l, obs_quad=float(dq_obs[0]), obs_ltte=float(dl_obs[0]), grid=[])
aG = np.exp(np.linspace(np.log(0.3), np.log(1000), 45))
for ecc_mode in ('circular', 'thermal'):
    for a in aG:
        m = 1500
        Mt = M3 + MEB
        Pout = np.sqrt(a ** 3 / Mt)   # yr
        cosi = rng.uniform(0, 1, m); sini = np.sqrt(1 - cosi ** 2)
        aAB = a * M3 / Mt * AU
        if ecc_mode == 'circular':
            e = np.zeros(m)
        else:
            e = np.sqrt(rng.uniform(0, 1, m)) * 0.95
        w = rng.uniform(0, 2 * np.pi, m); Tp = rng.uniform(0, 1, m) * Pout
        Mmean = 2 * np.pi * ((tyr[:, None] - Tp[None, :]) / Pout) % (2 * np.pi)
        Ea = kepE(Mmean, e[None, :])
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Ea / 2), np.sqrt(1 - e) * np.cos(Ea / 2))
        rr = (1 - e ** 2) / (1 + e * np.cos(nu))
        tau = aAB * sini / C * rr * np.sin(nu + w) / 60.0   # minutes
        Y = base[:, None] + tau + noise(m)
        dq, dl = stats(Y)
        det = (dq > thr_q) | (dl > thr_l)
        semi_amp = aAB * sini / C / 60
        out['grid'].append(dict(ecc=ecc_mode, a_au=float(a), Pout_yr=float(Pout), frac_excluded=float(det.mean()),
                                median_LTTE_semiamp_min=float(np.median(semi_amp))))
        print(f'{ecc_mode:8s} a={a:7.2f} au  Pout={Pout:8.2f} yr  LTTE semi-amp (median sin i) {np.median(semi_amp):7.2f} min  excluded fraction {det.mean():.3f}', flush=True)
json.dump(out, open(f'{HERE}/results/etv_companion_{tag}_M3{M3}_MEB{MEB}.json', 'w'), indent=1)
