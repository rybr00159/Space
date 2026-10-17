"""s10: RV design to close the BD-on-M3V channels (the analyst's b11 was not run).
For a BD of mass M on the M3V (0.370 Msun, P 5.6915020 d, sin i = 1): K = 28.4329 m/s (M/MJ)(P/yr)^-1/3 ((M*+M)/Msun)^-2/3
(1-e^2)^-1/2. Orbits: (a) the grazing configuration the light curve needs (rho_c/rho* ~0.19: e 0.5-0.6, w 250-290 deg,
transit near apoastron); (b) the general BD e prior (Ferreira+26 short-period Beta, w uniform). Epoch strategies at
fixed phases after mid-transit (final ephemeris; phase error < 0.001 in 2027): A = {0.25, 0.75}; B = {0.2, 0.5, 0.8};
C = 2 random phases; D = 4 random phases. Detection: chi2 of the epochs about their weighted mean > chi2_(N-1)(FAP 1e-3).
Output: the per-epoch precision sigma at which >= 99% (and >= 99.9%) of BDs of mass >= 13 MJ (and of the empirical
M-dwarf BD mass distribution) are detected; and the planet K for reference.
Output: out/s10_rv_design.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
from scipy.stats import chi2 as C2
rng = np.random.default_rng(10)
def K_of(M, e):
    return 28.4329*M*(PC/365.25)**(-1/3.)*((MSTAR + M*MJUP/MSUN))**(-2/3.)/np.sqrt(1 - e**2)/1e3   # km/s
def rv_at_phases(K, e, w, ph):
    """RV at times t_tr + ph*P; transit at true anomaly f_tr = pi/2 - w"""
    ftr = np.pi/2 - w
    Etr = 2*np.arctan2(np.sqrt(1 - e)*np.sin(ftr/2), np.sqrt(1 + e)*np.cos(ftr/2))
    Mtr = Etr - e*np.sin(Etr)
    M = Mtr[:, None] + 2*np.pi*ph[None, :] if np.ndim(Mtr) else Mtr + 2*np.pi*ph
    E = np.array(M, float)
    for _ in range(50): E = E - (E - (e*np.sin(E).T).T - M)/(1 - (e*np.cos(E).T).T)
    f = 2*np.arctan2((np.sqrt(1 + e)*np.sin(E/2).T).T, (np.sqrt(1 - e)*np.cos(E/2).T).T)
    return (K*(np.cos((f.T + w)) + e*np.cos(w))).T
N = 200000
emp = np.array([33.7, 68.1, 54.6, 75.5, 61.6, 64.4, 77.0, 69.0, 40.4, 68.0, 62.1, 18.5])
pops = {}
e_g = rng.uniform(0.5, 0.6, N); w_g = np.radians(rng.uniform(250, 290, N))
a = np.clip(rng.normal(0.426, 0.092, N), 0.05, None); b = np.clip(rng.normal(2.354, 0.73, N), 0.1, None)
e_p = np.clip(rng.beta(a, b), 0, 0.9); w_p = rng.uniform(0, 2*np.pi, N)
for orb, (e, w) in {'grazing_config': (e_g, w_g), 'BD_e_prior': (e_p, w_p)}.items():
    for mname, M in {'13MJ': np.full(N, 13.0), 'empirical': rng.choice(emp, N)*rng.normal(1, 0.05, N)}.items():
        pops[f'{orb}|{mname}'] = (K_of(M, e), e, w)
strat = {'A_0.25_0.75': np.array([0.25, 0.75]), 'B_0.2_0.5_0.8': np.array([0.2, 0.5, 0.8])}
out = {'K_kms': {f'{m}MJ_e{e}': float(K_of(m, e)) for m in (13, 20, 40, 60, 80) for e in (0.0, 0.55)}}
sigs = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
res = {}
for pk, (K, e, w) in pops.items():
    for sk, ph in list(strat.items()) + [('C_2_random', None), ('D_4_random', None)]:
        if ph is None:
            n = 2 if sk.startswith('C') else 4
            phs = rng.uniform(0, 1, (N, n))
            v = np.stack([rv_at_phases(K, e, w, phs[:, j])[..., 0] if False else None for j in range(0)]) if False else None
            # vectorised per-epoch evaluation
            V = []
            for j in range(n):
                ftr = np.pi/2 - w
                Etr = 2*np.arctan2(np.sqrt(1 - e)*np.sin(ftr/2), np.sqrt(1 + e)*np.cos(ftr/2)); Mtr = Etr - e*np.sin(Etr)
                Mj = Mtr + 2*np.pi*phs[:, j]; E = Mj.copy()
                for _ in range(50): E = E - (E - e*np.sin(E) - Mj)/(1 - e*np.cos(E))
                f = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
                V.append(K*(np.cos(f + w) + e*np.cos(w)))
            V = np.array(V).T
        else:
            V = []
            for p_ in ph:
                ftr = np.pi/2 - w
                Etr = 2*np.arctan2(np.sqrt(1 - e)*np.sin(ftr/2), np.sqrt(1 + e)*np.cos(ftr/2)); Mtr = Etr - e*np.sin(Etr)
                Mj = Mtr + 2*np.pi*p_; E = Mj.copy()
                for _ in range(50): E = E - (E - e*np.sin(E) - Mj)/(1 - e*np.cos(E))
                f = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
                V.append(K*(np.cos(f + w) + e*np.cos(w)))
            V = np.array(V).T
        n = V.shape[1]
        dev2 = np.sum((V - V.mean(axis=1, keepdims=True))**2, axis=1)
        thr = C2.ppf(1 - 1e-3, n - 1)
        det = {}
        for s in sigs:
            # noise realisations: chi2 = sum((v + noise - mean)^2)/s^2 ~ noncentral chi2(n-1, dev2/s^2)
            lam = dev2/s**2
            ch = rng.noncentral_chisquare(n - 1, np.maximum(lam, 1e-12))
            det[f'{s}'] = float(np.mean(ch > thr))
        sreq99 = float(np.interp(0.99, [det[f'{s}'] for s in sigs[::-1]], sigs[::-1])) if det[f'{sigs[0]}'] >= 0.99 else None
        sreq999 = float(np.interp(0.999, [det[f'{s}'] for s in sigs[::-1]], sigs[::-1])) if det[f'{sigs[0]}'] >= 0.999 else None
        res[f'{pk}|{sk}'] = dict(det_frac_by_sigma=det, sigma_for_99pct=sreq99, sigma_for_99p9pct=sreq999,
                                 min_rms_dev=float(np.sqrt(np.min(dev2)/n)), p01_rms_dev=float(np.sqrt(np.percentile(dev2, 1)/n)))
        print(pk, sk, 'sigma(99%)', sreq99, 'sigma(99.9%)', sreq999, {k: round(v, 3) for k, v in det.items()}, flush=True)
out['design'] = res
out['planet_K_m_s'] = {f'{m}ME': float(28.4329*m/317.83*(PC/365.25)**(-1/3.)*MSTAR**(-2/3.)) for m in (10, 20, 30)}
print(out['K_kms'], out['planet_K_m_s'])
jdump(dict(method=__doc__, **out), f'{R3}/out/s10_rv_design.json')
