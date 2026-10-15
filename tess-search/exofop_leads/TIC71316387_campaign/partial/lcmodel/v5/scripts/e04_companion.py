"""Constraints on an outer companion (the G star, M3) from the eclipse timings via the light-travel-time effect.
For a grid of outer semi-major axes a_out (au) and random orbital elements (cos i uniform, e uniform in [0, emax],
omega, mean anomaly at the reference epoch uniform), the LTTE delay of the EB at each observed eclipse epoch is
  dt = (a_EB sin i / c) [ (1 - e^2)/(1 + e cos nu) sin(nu + omega) ],  a_EB = a_out M3/(M_EB + M3)
and the orbit is EXCLUDED if subtracting it from the measured O-C (with T0, P and the secondary offset refitted)
raises chi2 by more than DCHI (default 9) above the best linear fit.  The same for the 'detectable' criterion on
noise-free synthetic data (chi2 of the LTTE signal after removing the best linear fit > DCHI).
Also reported: the curvature (quadratic coefficient) each orbit would produce over 2018-2026 and the Gaia DR3 RV
change of the G star over the DR3 baseline (34 months) for each orbit, for reference.
Usage: e04_companion.py <etvfit json> <M_EB> <M3> [DCHI] [emax]"""
import numpy as np, sys, json
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
F = json.load(open(sys.argv[1])); MEB = float(sys.argv[2]); M3 = float(sys.argv[3])
DCHI = float(sys.argv[4]) if len(sys.argv) > 4 else 9.0; EMAX = float(sys.argv[5]) if len(sys.argv) > 5 else 0.8
evs = F['events_used']
n = np.array([e['n'] for e in evs], float); ph = np.array([e['ph'] for e in evs]); oc = np.array([e['oc_min'] for e in evs]) / 1440
sig = np.array([e['sig_min'] for e in evs]) / 1440
P = F['fits']['lin']['coef'][1]
T0 = F['fits']['lin']['coef'][0]
t = T0 + n * P + ph * P           # BTJD of each event (approx.)
A = np.vstack([np.ones_like(n), n, (ph == 0.5).astype(float)]).T
Aw = A / sig[:, None]
proj = A @ np.linalg.pinv(Aw) / sig[None, :]   # maps data -> best linear fit
def chi2_after_lin(y):
    r = y - proj @ y
    return float(((r / sig) ** 2).sum())
chi_lin = chi2_after_lin(oc)
# conservative: rescale errors so that the best linear fit has reduced chi2 = 1 (if > 1)
fac = max(1.0, np.sqrt(chi_lin / (n.size - 3)))
sig = sig * fac
Aw = A / sig[:, None]; proj = A @ np.linalg.pinv(Aw) / sig[None, :]
chi_lin = chi2_after_lin(oc)
print(f'error rescale factor {fac:.3f}; chi2_lin after rescale {chi_lin:.1f} for {n.size} events')
AU_D = 499.004784 / 86400   # light time of 1 au in days
rng = np.random.default_rng(11)
NORB = 4000
a_grid = np.array([3, 5, 7, 10, 14, 20, 25, 30, 40, 50, 60, 80, 100, 130, 160, 200, 300, 500])
Mtot = MEB + M3
out = []
for a in a_grid:
    Pout = np.sqrt(a ** 3 / Mtot) * 365.25   # days
    aEB = a * M3 / Mtot
    excl = 0; det = 0; curv = []; drv = []
    for k in range(NORB):
        cosi = rng.uniform(0, 1); sini = np.sqrt(1 - cosi ** 2)
        e = rng.uniform(0, EMAX); w = rng.uniform(0, 2 * np.pi); M0 = rng.uniform(0, 2 * np.pi)
        M = M0 + 2 * np.pi * (t - t.mean()) / Pout
        E = M.copy()
        for _ in range(60): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        dt = aEB * sini * AU_D * (1 - e ** 2) / (1 + e * np.cos(nu)) * np.sin(nu + w)
        c_sig = chi2_after_lin(dt)
        if c_sig > DCHI: det += 1
        if chi2_after_lin(oc - dt) - chi_lin > DCHI: excl += 1
        # quadratic coefficient of the LTTE curve over the data (min / yr^2)
        x = (t - t.mean()) / 365.25
        curv.append(np.polyfit(x, dt * 1440, 2)[0])
        # Gaia DR3 RV change of the G star over 2014.6-2017.4 (approx.) : K3 * [cos(nu+w) + e cos w] difference
        tg = np.array([-136.5, 901.5])  # BTJD of the Gaia DR3 baseline start (2014-07-25) and end (2017-05-28)
        Mg = M0 + 2 * np.pi * (tg - t.mean()) / Pout
        Eg = Mg.copy()
        for _ in range(60): Eg = Eg - (Eg - e * np.sin(Eg) - Mg) / (1 - e * np.cos(Eg))
        nug = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Eg / 2), np.sqrt(1 - e) * np.cos(Eg / 2))
        K3 = 2 * np.pi * (a * MEB / Mtot) * 1.495978707e8 * sini / (Pout * 86400 * np.sqrt(1 - e ** 2))
        v = K3 * (np.cos(nug + w) + e * np.cos(w))
        drv.append(abs(v[1] - v[0]))
    curv = np.abs(np.array(curv)); drv = np.array(drv)
    rec = dict(a_au=float(a), P_out_yr=float(Pout / 365.25), ltte_semiamp_max_min=float(aEB * AU_D * 1440),
               frac_excluded=excl / NORB, frac_detectable=det / NORB, median_abs_curv_min_per_yr2=float(np.median(curv)),
               p90_abs_curv=float(np.percentile(curv, 90)), median_dRV_gaia_kms=float(np.median(drv)), frac_dRV_gt3=float(np.mean(drv > 3)))
    out.append(rec)
    print(f"a={a:5.0f} au  P_out={Pout/365.25:7.1f} yr  max LTTE amp {rec['ltte_semiamp_max_min']:7.1f} min  excluded {rec['frac_excluded']:.3f}  "
          f"detectable {rec['frac_detectable']:.3f}  median |curv| {rec['median_abs_curv_min_per_yr2']:.3f} min/yr^2  median Gaia dRV {rec['median_dRV_gaia_kms']:.2f} km/s", flush=True)
json.dump(dict(M_EB=MEB, M3=M3, DCHI=DCHI, EMAX=EMAX, chi2_lin=chi_lin, err_rescale=fac, n_events=int(n.size), grid=out),
          open(f'{V}/results/companion_MEB{MEB:g}_M3{M3:g}_d{DCHI:g}_e{EMAX:g}.json', 'w'), indent=1)
