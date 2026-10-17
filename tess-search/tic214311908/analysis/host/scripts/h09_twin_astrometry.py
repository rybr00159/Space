"""h09: Gaia astrometric constraints specifically for NEAR-EQUAL (twin) companions, q = M2/M1 = 0.80-1.00,
with both stars on the main sequence (M1 ~ 0.27 Msun, total G light fixed to the observed M_G = 10.24).
Same physics as the phase-1 s12c simulation (independently re-coded): Gaia-like DR3 epochs (36 FoV transits in 18
visibility periods, 2014.65-2017.40), Keplerian photocentre orbit with beta = q/(1+q) - eps/(1+eps), eps = flux ratio
in G from the Mamajek MS table (M_G vs mass), 5-parameter AL fits for DR3 (all) and DR2 (t < 2016.38),
per-transit AL sigma calibrated to the DR3 parallax error, RUWE_0 drawn from field M dwarfs (h03 sample,
G 13.5-15.5, RUWE < 1.25, ipd_frac_multi_peak <= 2, 2.5 < BP-RP < 3.0; falls back to the phase-1 star/data/mdwarf_sample.csv if h03 is absent).
L_RUWE = P(RUWE_sim <= 1.10) / P(RUWE_0 <= 1.10);  L_PM from the calibrated DR3-DR2 PM difference (phase-1 star/s12c
values: dmu = (0.396-0.033, -0.236-0.008) mas/yr, sigma = (0.118*1.50, 0.132*1.60)); joint L per simulated orbit.
Output: out/h09_twin_astrometry.json, figs/h09_twin_astrometry.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
B = 'tic214311908/analysis/host'
S = 'tic214311908/analysis/star'
rng = np.random.default_rng(21)
ra, dec = np.radians(260.969674), np.radians(-51.323821)
plx = 15.5494; NCCD = 328
import os
m = pd.read_csv(f'{B}/data/field_mdwarfs.csv') if os.path.exists(f'{B}/data/field_mdwarfs.csv') else pd.read_csv(f'{S}/data/mdwarf_sample.csv')
r0 = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.25) & (m.ipd_frac_multi_peak <= 2) & (m.bp_rp.between(2.5, 3.0))].ruwe.values
pR0 = np.mean(r0 <= 1.10)
print('RUWE_0 sample', len(r0), 'P(<=1.10)', pR0, 'percentile of 1.076', np.mean(r0 <= 1.076))
tab = []
for line in open(f'{S}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] in 'KM' and p[0].endswith('V'):
        try:
            tab.append((float(p[30]), float(p[13])))
        except Exception:
            pass
tab = np.array(sorted(tab))
MGm = lambda M: np.interp(M, tab[:, 0], tab[:, 1])
# choose M1 for each q such that the combined M_G equals 10.24 on the Mamajek scale shifted to the observed twin case
def comb_MG(M1, q):
    return -2.5 * np.log10(10 ** (-0.4 * MGm(M1)) + 10 ** (-0.4 * MGm(q * M1)))
target_comb = comb_MG(0.264, 1.0)   # equal twin of 0.264 Msun (the star_params twin case) as reference
def M1_for_q(q):
    g = np.linspace(0.2, 0.45, 2001); v = comb_MG(g, q); return float(g[np.argmin(np.abs(v - target_comb))])
def earth(tyr):
    p = get_body_barycentric('earth', Time(tyr, format='jyear', scale='tcb'))
    return p.x.to('au').value, p.y.to('au').value, p.z.to('au').value
def epochs():
    vp = np.sort(rng.uniform(2014.65, 2017.40, 18))
    ang = rng.uniform(0, 2 * np.pi, 18)
    t = np.concatenate([vp, vp + 106.5 / 60 / 24 / 365.25]); psi = np.concatenate([ang, ang + np.radians(rng.normal(0, 2, 18))])
    o = np.argsort(t); return t[o], psi[o]
def design(t, psi):
    X, Y, Z = earth(t)
    pa = X * np.sin(ra) - Y * np.cos(ra); pd_ = X * np.cos(ra) * np.sin(dec) + Y * np.sin(ra) * np.sin(dec) - Z * np.cos(dec)
    s, c = np.sin(psi), np.cos(psi); dt = t - 2016.0
    return np.column_stack([s, c, s * dt, c * dt, pa * s + pd_ * c])
ss = []
for k in range(50):
    t, psi = epochs(); A = design(t, psi); ss.append(np.sqrt(np.linalg.inv(A.T @ A)[4, 4]))
sig_tr = 0.0251 / np.median(ss); print('per-transit sigma', sig_tr)
dmu_obs = np.array([0.396 - 0.033, -0.236 - 0.008]); sig_cal = np.array([0.118 * 1.50, 0.132 * 1.60])
den = np.exp(-0.5 * np.sum((dmu_obs / sig_cal) ** 2))
def kepler(Mn, e):
    E = Mn.copy()
    for _ in range(30): E = E - (E - e * np.sin(E) - Mn) / (1 - e * np.cos(E))
    return E
EP = []
for _ in range(25):
    t, psi = epochs(); A3 = design(t, psi); m2 = t < 2016.38
    EP.append((t, psi, A3, np.linalg.pinv(A3), m2, np.linalg.pinv(A3[m2])))
qg = np.array([0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.98, 0.99, 1.00])
ag = np.logspace(-2, 1.5, 29)
NS = 500
out = dict(q=qg.tolist(), a_au=ag.tolist(), M1=[], beta=[], eps=[], L_ruwe=[], L_pm=[], L_joint=[], med_ruwe=[])
for q in qg:
    M1 = M1_for_q(q); M2 = q * M1
    eps = 10 ** (-0.4 * (MGm(M2) - MGm(M1))); beta = q / (1 + q) - eps / (1 + eps)
    LR, LP, LJ, MR = [], [], [], []
    for a in ag:
        P = np.sqrt(a ** 3 / (M1 + M2)); alpha = abs(beta) * a * plx
        dchi = np.zeros(NS); dmu = np.zeros((NS, 2))
        for k in range(NS):
            t, psi, A3, A3p, m2, A2p = EP[k % 25]
            e = 0.0 if P * 365.25 < 10 else rng.uniform(0, 0.6)
            inc = np.arccos(rng.uniform(-1, 1)); om = rng.uniform(0, 2 * np.pi); Om = rng.uniform(0, 2 * np.pi); tp = rng.uniform(0, P)
            E = kepler(2 * np.pi * (t - tp) / P, e); X = np.cos(E) - e; Yv = np.sqrt(1 - e ** 2) * np.sin(E)
            Ati = np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc); Bti = np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc)
            Fti = -np.sin(om) * np.cos(Om) - np.cos(om) * np.sin(Om) * np.cos(inc); Gti = -np.sin(om) * np.sin(Om) + np.cos(om) * np.cos(Om) * np.cos(inc)
            al = (alpha * (Bti * X + Gti * Yv)) * np.sin(psi) + (alpha * (Ati * X + Fti * Yv)) * np.cos(psi)
            p3 = A3p @ al; dchi[k] = np.sum((al - A3 @ p3) ** 2) / sig_tr ** 2
            dmu[k] = p3[2:4] - (A2p @ al[m2])[2:4]
        ruwe = np.sqrt(rng.choice(r0, NS) ** 2 + dchi / (NCCD - 5))
        pm = np.exp(-0.5 * np.sum(((dmu_obs - dmu) / sig_cal) ** 2, axis=1))
        LR.append(float((np.mean(ruwe <= 1.10) + 1e-4) / pR0)); LP.append(float(pm.mean() / den))
        LJ.append(float(((ruwe <= 1.10) * pm).mean() / (pR0 * den) + 1e-5)); MR.append(float(np.median(ruwe)))
    out['M1'].append(M1); out['beta'].append(float(beta)); out['eps'].append(float(eps))
    out['L_ruwe'].append(LR); out['L_pm'].append(LP); out['L_joint'].append(LJ); out['med_ruwe'].append(MR)
    print(f'q {q:.2f}: M1 {M1:.3f} eps_G {eps:.3f} beta {beta:.4f}; L_joint at a = 0.03/0.1/0.3/1/3/10 au: ' +
          ' '.join(f'{np.interp(np.log10(x), np.log10(ag), LJ):.3f}' for x in (0.03, 0.1, 0.3, 1, 3, 10)), flush=True)
out['pR0'] = float(pR0); out['sig_tr'] = float(sig_tr)
json.dump(out, open(f'{B}/out/h09_twin_astrometry.json', 'w'), indent=1)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4.5))
L = np.array(out['L_joint'])
im = ax.pcolormesh(ag, qg, np.log10(np.clip(L, 1e-4, 10)), vmin=-3, vmax=0.5, cmap='RdBu', shading='nearest')
ax.contour(ag, qg, L, levels=[0.05, 0.3], colors=['k', 'gray'])
ax.set_xscale('log'); ax.set_xlabel('semi-major axis (au)'); ax.set_ylabel('q'); ax.set_title('Gaia RUWE + DR3-DR2 PM: log10 L(twin)/L(single)')
fig.colorbar(im); fig.savefig(f'{B}/figs/h09_twin_astrometry.png', dpi=110, bbox_inches='tight')
