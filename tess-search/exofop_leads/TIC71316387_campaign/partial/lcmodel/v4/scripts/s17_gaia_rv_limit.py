"""Context (approximate): if the G star (M3 = 1.03 Msun) orbits the EB (M_EB), its own reflex RV must be consistent with the
constant Gaia DR3 RV (15 transits, rv_chisq_pvalue = 0.275).  Per-transit RV error inferred from radial_velocity_error =
0.857 km/s = sqrt(pi/(2N)) * eps + 0.11 (DR3 combined-method convention) -> eps ~ 2.3 km/s.  Transit epochs are unknown;
we draw 15 epochs uniformly in the DR3 window (2014-07-25 .. 2017-05-28).  For each separation: isotropic orientation, random
phase, circular or thermal e; a configuration is 'excluded' if the noiseless chi2 of its RV curve about the weighted mean
exceeds the chi2 level that DR3 would have flagged (p < 0.001 for 14 dof, i.e. chi2 > 36.1), given that the observed p = 0.275.
This is an approximation (epoch sampling assumed), reported only as context for the ETV limits."""
import numpy as np, json, sys
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
M3 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.03; MEB = float(sys.argv[2]) if len(sys.argv) > 2 else 0.60
eps = (0.857 - 0.11) / np.sqrt(np.pi / (2 * 15)); N = 15
t0, t1 = 2014.56, 2017.41
from scipy.stats import chi2 as C2
thr = C2.ppf(1 - 1e-3, N - 1) - (N - 1)   # excess over the expectation for pure noise
rng = np.random.default_rng(9)
out = []
for mode in ('circular', 'thermal'):
    for a in np.exp(np.linspace(np.log(0.5), np.log(300), 40)):
        Mt = M3 + MEB; P = np.sqrt(a ** 3 / Mt); m = 4000
        K = 29.78 * np.sqrt(Mt / a) * MEB / Mt          # km/s, G-star reflex semi-amplitude for sin i = 1, circular
        cosi = rng.uniform(0, 1, m); sini = np.sqrt(1 - cosi ** 2)
        e = np.zeros(m) if mode == 'circular' else np.minimum(np.sqrt(rng.uniform(0, 1, m)), 0.95)
        w = rng.uniform(0, 2 * np.pi, m); tp = rng.uniform(0, P, m)
        tt = np.sort(rng.uniform(t0, t1, (m, N)), axis=1)
        Mm = 2 * np.pi * (tt - tp[:, None]) / P
        E = Mm.copy()
        for _ in range(40): E = E - (E - e[:, None] * np.sin(E) - Mm) / (1 - e[:, None] * np.cos(E))
        nu = 2 * np.arctan2(np.sqrt(1 + e[:, None]) * np.sin(E / 2), np.sqrt(1 - e[:, None]) * np.cos(E / 2))
        rv = K * sini[:, None] / np.sqrt(1 - e[:, None] ** 2) * (np.cos(nu + w[:, None]) + e[:, None] * np.cos(w[:, None]))
        dchi = np.sum((rv - rv.mean(1, keepdims=True)) ** 2, 1) / eps ** 2
        out.append(dict(ecc=mode, a_au=float(a), Pout_yr=float(P), K_sini1_kms=float(K), excluded_frac=float(np.mean(dchi > thr))))
        print(f'{mode:8s} a={a:7.2f} au P={P:8.2f} yr K(sin i=1)={K:6.2f} km/s  excluded fraction {np.mean(dchi > thr):.3f}')
json.dump(dict(eps_kms=float(eps), thr=float(thr), grid=out), open(f'{FD}/results/gaia_rv_limit.json', 'w'), indent=1)
