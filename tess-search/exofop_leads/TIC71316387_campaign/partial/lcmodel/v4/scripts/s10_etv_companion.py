"""Limits on an outer companion (mass M3, default the G star, 1.03 Msun) from the 2018-2026 eclipse timings via the
light-travel-time effect (LTTE) of the EB's barycentric motion.
For each relative semi-major axis a (au): Monte Carlo over isotropic inclination, uniform orbital phase, and eccentricity
(circular, or 'thermal' f(e)=2e, e<0.95, random omega).  Two criteria per configuration:
 (A) INCONSISTENT with the data: chi2 of the measured O-C after subtracting that configuration's LTTE curve (linear ephemeris
     and secondary offset re-fitted) exceeds the linear-only chi2 by > 9 (i.e. that specific orbit is rejected at >3 sigma,
     using the injection-calibrated per-event errors);
 (B) DETECTABLE: with the configuration's LTTE added to the linear ephemeris plus null noise (normalised residuals with
     group means removed, bootstrapped), the quadratic or best-sinusoid Delta-chi2 exceeds its 99.7% null threshold.
Reports, per a, the fraction of orientations/phases that are excluded (A) or detectable (B).
Usage: s10_etv_companion.py <etvfit json> <M3> <M_EB>"""
import numpy as np, json, sys
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
D = json.load(open(sys.argv[1])); M3 = float(sys.argv[2]); MEB = float(sys.argv[3])
E = np.array(D['E']); oc = np.array(D['oc']); es = np.array(D['es']); sec = np.array(D['sec']); ph = np.array(D['ph']); t = np.array(D['t'])
tref = D['tref']; tyr = (t - tref) / 365.25; n = len(oc); W = 1 / es
def hat(A):
    Aw = A * W[:, None]; return Aw @ np.linalg.pinv(Aw)
Al = np.array([np.ones(n), E / 1000, (ph == 0.5).astype(float)]).T
Hl = hat(Al); Hq = hat(np.column_stack([Al, tyr ** 2]))
Pg = np.exp(np.linspace(np.log(0.5), np.log(60), 100))
HL = [hat(np.column_stack([Al, np.sin(2 * np.pi * tyr / P), np.cos(2 * np.pi * tyr / P)])) for P in Pg]
def chi_lin(Y):
    Yw = Y * W[:, None]; r = Yw - Hl @ Yw; return np.sum(r ** 2, 0)
def stats(Y):
    Yw = Y * W[:, None]; cl = np.sum((Yw - Hl @ Yw) ** 2, 0); cq = np.sum((Yw - Hq @ Yw) ** 2, 0)
    best = np.zeros(Y.shape[1])
    for H in HL: best = np.maximum(best, cl - np.sum((Yw - H @ Yw) ** 2, 0))
    return cl - cq, best
cfit = np.linalg.lstsq(Al * W[:, None], oc * W, rcond=None)[0]; base = Al @ cfit
z = (oc - base) / es
groups = [np.where((sec == s) & (ph == f))[0] for s in np.unique(sec) for f in (0.0, 0.5)]
zpool = np.concatenate([(z[g] - z[g].mean()) * np.sqrt(g.size / (g.size - 1)) for g in groups if g.size > 1])
rng = np.random.default_rng(5)
def noise(m): return rng.choice(zpool, (n, m), replace=True) * es[:, None]
dq0, dl0 = stats(base[:, None] + noise(4000))
thr_q = float(np.percentile(dq0, 99.7)); thr_l = float(np.percentile(dl0, 99.7))
chi0 = float(chi_lin(oc[:, None])[0])
print(f'null 99.7% thresholds: quadratic {thr_q:.2f}, sinusoid-grid {thr_l:.2f}; linear chi2 of data {chi0:.1f} (N={n})')
AU = 1.495978707e11; C = 2.99792458e8
def kepE(Mm, e):
    Ea = Mm + e * np.sin(Mm)
    for _ in range(60): Ea = Ea - (Ea - e * np.sin(Ea) - Mm) / (1 - e * np.cos(Ea))
    return Ea
out = dict(M3=M3, MEB=MEB, thr_quad=thr_q, thr_sin=thr_l, grid=[])
aG = np.exp(np.linspace(np.log(0.5), np.log(2000), 50))
for mode in ('circular', 'thermal'):
    for a in aG:
        m = 2000; Mt = M3 + MEB; Pout = np.sqrt(a ** 3 / Mt)
        cosi = rng.uniform(0, 1, m); sini = np.sqrt(1 - cosi ** 2)
        e = np.zeros(m) if mode == 'circular' else np.minimum(np.sqrt(rng.uniform(0, 1, m)), 0.95)
        w = rng.uniform(0, 2 * np.pi, m); Tp = rng.uniform(0, 1, m) * Pout
        Mm = (2 * np.pi * (tyr[:, None] - Tp[None, :]) / Pout) % (2 * np.pi)
        Ea = kepE(Mm, e[None, :]); nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Ea / 2), np.sqrt(1 - e) * np.cos(Ea / 2))
        rr = (1 - e ** 2) / (1 + e * np.cos(nu))
        aEB = a * M3 / Mt * AU
        tau = aEB * sini / C * rr * np.sin(nu + w) / 60.0     # minutes
        incons = chi_lin(oc[:, None] - tau) - chi0 > 9.0
        dq, dl = stats(base[:, None] + tau + noise(m))
        det = (dq > thr_q) | (dl > thr_l)
        span = tau.max(0) - tau.min(0)
        # curvature measure: max |deviation from the best straight line through the LTTE curve| over 2018-2026
        lin = np.polyfit(tyr, tau, 1); dev = np.abs(tau - (np.outer(tyr, lin[0]) + lin[1])).max(0)
        out['grid'].append(dict(ecc=mode, a_au=float(a), Pout_yr=float(Pout), excluded_frac=float(incons.mean()), detect_frac=float(det.mean()),
                                ltte_semiamp_med_min=float(np.median(aEB * sini / C / 60)), dev_from_line_med_min=float(np.median(dev))))
        print(f'{mode:8s} a={a:8.2f} au Pout={Pout:9.1f} yr  LTTE semi-amp (median) {np.median(aEB*sini/C/60):8.2f} min  '
              f'median 2018-26 departure from a straight line {np.median(dev):7.2f} min | excluded (3 sigma) {incons.mean():.3f}  detectable {det.mean():.3f}', flush=True)
json.dump(out, open(f'{FD}/results/etv_companion_M3{M3}_MEB{MEB}.json', 'w'), indent=1)
for mode in ('circular', 'thermal'):
    g = [x for x in out['grid'] if x['ecc'] == mode]
    for f in (0.5, 0.9):
        aa = [x['a_au'] for x in g if x['excluded_frac'] >= f]
        print(f'{mode}: excluded fraction >= {f:.0%} for a <= {max(aa) if aa else float("nan"):.1f} au (largest such grid a)')
