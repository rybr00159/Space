"""Referee: independent spot check of the RUWE-based companion exclusion (s12c claims 0.08-0.3 Msun excluded at
~0.1-4 au).  Own simplified simulation: 41 field-of-view transits (each 9 CCD AL measurements with per-CCD error
sigma_CCD) at random times over 2014.6-2017.4 with random scan angles; the target's actual parallax factors are
ignored (the 5-parameter fit is linear: offsets, PMs and a parallax-like term built from the Earth's orbit).
Photocentre orbit: circular, random inclination/orientation/phase; alpha = a * (q/(1+q) - eps/(1+eps)) * plx,
eps = G-band flux ratio from the Mamajek EEM table (M_G vs mass), M1 = 0.37, MG1 = 10.24.
RUWE_sim = sqrt(RUWE0^2 + chi2_excess/(N-5)), RUWE0 from own comparison sample (G 13.5-15.5, RUWE < 1.25).
Reports P(RUWE_sim <= 1.10) / P(RUWE0 <= 1.10) for sigma_CCD = 0.22 mas (investigator-equivalent) and 0.35 mas (conservative).
Output: out/r10_ruwe_check.json"""
import json, numpy as np, pandas as pd
R = 'tic214311908/analysis/referee_star'
rng = np.random.default_rng(3)
m = pd.read_csv(f'{R}/data/mdw_cmd.csv')
r0 = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.25) & (m.ipd_frac_multi_peak <= 2)].ruwe.values
pR0 = np.mean(r0 <= 1.10); print('RUWE0 sample', len(r0), 'median', np.median(r0).round(3), 'P(<=1.10)', pR0.round(3), 'pctile of 1.076', np.mean(r0 <= 1.076).round(3))
# Mamajek mass -> MG
hdr = [l for l in open(f'{R}/data/EEM_dwarf.txt') if l.startswith('#SpT')][0].split()
iM, iG = hdr.index('Msun'), hdr.index('M_G'); tab = []
for line in open(f'{R}/data/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > iM and p[0][:1] in 'KM' and p[0].endswith('V'):
        try: tab.append((float(p[iM]), float(p[iG])))
        except Exception: pass
tab = np.array(sorted(tab))
def eps(M2): return 0.0 if M2 < 0.08 else 10 ** (-0.4 * (np.interp(M2, tab[:, 0], tab[:, 1]) - 10.24))
plx = 15.55; M1 = 0.37; NTR = 41; NCCD = 9
def sim(M2, a, sig, nsim=400):
    e_ = eps(M2); q = M2 / M1; beta = q / (1 + q) - e_ / (1 + e_); alpha = beta * a * plx; P = np.sqrt(a ** 3 / (M1 + M2))
    out = []
    for k in range(nsim):
        t = np.sort(rng.uniform(2014.6, 2017.4, NTR)); psi = rng.uniform(0, 2 * np.pi, NTR)
        # parallax-like factor (circular Earth orbit, ecliptic-ish approximation)
        lam = 2 * np.pi * (t - 2015.0); pf = np.sin(lam) * np.sin(psi) + 0.6 * np.cos(lam) * np.cos(psi)
        A = np.column_stack([np.sin(psi), np.cos(psi), np.sin(psi) * (t - 2016), np.cos(psi) * (t - 2016), pf])
        inc = np.arccos(rng.uniform(-1, 1)); Om = rng.uniform(0, 2 * np.pi); ph = rng.uniform(0, 2 * np.pi)
        u = 2 * np.pi * t / P + ph
        x, y = np.cos(u), np.sin(u) * np.cos(inc)
        dE = alpha * (x * np.cos(Om) - y * np.sin(Om)); dN = alpha * (x * np.sin(Om) + y * np.cos(Om))
        al = dE * np.sin(psi) + dN * np.cos(psi)
        res = al - A @ np.linalg.lstsq(A, al, rcond=None)[0]
        chi = NCCD * np.sum(res ** 2) / sig ** 2
        out.append(np.sqrt(rng.choice(r0) ** 2 + chi / (NTR * NCCD - 5)))
    out = np.array(out); return float((np.mean(out <= 1.10) + 1e-3) / pR0), float(np.median(out)), float(alpha)
res = {}
for sig in (0.22, 0.35):
    for M2 in (0.08, 0.1, 0.2, 0.3, 0.33):
        for a in (0.03, 0.1, 0.3, 1.0, 2.0, 4.0, 6.0):
            L, med, al = sim(M2, a, sig)
            res[f'sig{sig}_M{M2}_a{a}'] = dict(L=L, median_ruwe=med, alpha_mas=al)
        print(f'sigma_CCD {sig} M2 {M2}: ' + '  '.join(f"a={a}: L={res[f'sig{sig}_M{M2}_a{a}']['L']:.2f} (RUWE~{res[f'sig{sig}_M{M2}_a{a}']['median_ruwe']:.2f}, alpha {res[f'sig{sig}_M{M2}_a{a}']['alpha_mas']:.2f})" for a in (0.03, 0.1, 0.3, 1.0, 2.0, 4.0, 6.0)), flush=True)
json.dump(res, open(f'{R}/out/r10_ruwe_check.json', 'w'), indent=1)
