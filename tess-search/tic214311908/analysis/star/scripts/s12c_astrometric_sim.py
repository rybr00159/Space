"""Limits on bound companions to TIC 214311908 at 0.01-30 au from Gaia astrometry:
(1) DR3 RUWE = 1.076 and (2) the DR3 - DR2 proper-motion difference (0.396 +- 0.118, -0.236 +- 0.132 mas/yr;
calibrated on 1236 comparison M dwarfs: z-scatter x1.50 (RA), x1.60 (Dec), median offset +0.033 / +0.008 mas/yr).
Simulation: Gaia-like epochs (36 field-of-view transits in 18 visibility periods over 2014.65-2017.40, random scan
angle per visibility period, pairs of transits 106.5 min apart); along-scan (AL) photocentre displacements of a
Keplerian orbit (random orientation, phase; e uniform 0-0.6, e = 0 for P < 10 d) scaled by
beta = q/(1+q) - eps/(1+eps) (eps = companion/primary flux ratio in G from the Mamajek main-sequence table;
eps = 0 for M2 < 0.075 Msun); 5-parameter AL fit per data release (DR2: transits before 2016.38).
The per-transit AL error is set so that the simulated parallax error equals the DR3 value (0.0251 mas).
RUWE_sim^2 = RUWE_0^2 + dchi2/(N_CCD - 5), RUWE_0 drawn from comparison M dwarfs (G 13.5-15.5, RUWE < 1.25).
Likelihoods relative to no companion: L_RUWE = P(RUWE_sim <= 1.10)/P(RUWE_0 <= 1.10);
L_PM = mean over sims of N(dmu_obs | dmu_sim, sigma_cal) / N(dmu_obs | 0, sigma_cal). Excluded where L < 0.05.
Output: out/astrometric_limits.json, figs/astrometric_limits.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
B = 'tic214311908/analysis/star'
rng = np.random.default_rng(11)
ra, dec = np.radians(260.969674), np.radians(-51.323821)
plx = 15.5494; M1 = 0.37; NCCD = 328
# empirical RUWE_0
m = pd.read_csv(f'{B}/data/mdwarf_sample.csv')
r0 = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.25) & (m.ipd_frac_multi_peak <= 2)].ruwe.values
print('RUWE_0 sample', len(r0), 'median', np.median(r0), 'P(<=1.10)', np.mean(r0 <= 1.10), 'pct of 1.076:', np.mean(r0 <= 1.076))
# Mamajek MS: mass -> M_G
tab = []
for line in open(f'{B}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] in 'KM' and p[0].endswith('V'):
        try:
            mg = float(p[13]); ms = float(p[30]); tab.append((ms, mg))
        except Exception:
            pass
tab = np.array(sorted(tab)); print('Mamajek mass-MG rows:', len(tab), tab[:3].tolist(), tab[-3:].tolist())
MG1 = 10.233
def epsG(M2):
    if M2 < 0.075: return 0.0
    mg = np.interp(M2, tab[:, 0], tab[:, 1]); return 10 ** (-0.4 * (mg - MG1))
# Earth barycentric positions
def earth(tyr):
    t = Time(tyr, format='jyear', scale='tcb')
    p = get_body_barycentric('earth', t)
    return p.x.to('au').value, p.y.to('au').value, p.z.to('au').value
def epochs():
    t0, t1 = 2014.65, 2017.40
    vp = np.sort(rng.uniform(t0, t1, 18))
    t = np.concatenate([vp, vp + 106.5 / 60 / 24 / 365.25]); psi = np.concatenate([np.full(18, 1.0), np.full(18, 1.0)])
    ang = rng.uniform(0, 2 * np.pi, 18); psi = np.concatenate([ang, ang + np.radians(rng.normal(0, 2, 18))])
    o = np.argsort(t); return t[o], psi[o]
def design(t, psi):
    X, Y, Z = earth(t)
    pa = X * np.sin(ra) - Y * np.cos(ra); pd_ = X * np.cos(ra) * np.sin(dec) + Y * np.sin(ra) * np.sin(dec) - Z * np.cos(dec)
    s, c = np.sin(psi), np.cos(psi); dt = t - 2016.0
    return np.column_stack([s, c, s * dt, c * dt, pa * s + pd_ * c])   # (da*, dd, mua*, mud, plx)
# calibrate per-transit sigma so that sigma_plx = 0.0251
ss = []
for k in range(50):
    t, psi = epochs(); A = design(t, psi); C = np.linalg.inv(A.T @ A); ss.append(np.sqrt(C[4, 4]))
sig_tr = 0.0251 / np.median(ss)
print(f'per-transit AL sigma {sig_tr:.3f} mas (per CCD ~{sig_tr*3:.3f} mas)')
dmu_obs = np.array([0.396 - 0.033, -0.236 - 0.008]); sig_cal = np.array([0.118 * 1.50, 0.132 * 1.60])
def kepler(Mn, e):
    E = Mn.copy()
    for _ in range(30): E = E - (E - e * np.sin(E) - Mn) / (1 - e * np.cos(E))
    return E
M2g = np.array([0.0095, 0.02, 0.04, 0.06, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.34])
ag = np.logspace(-2, 1.5, 29)
NS = 400
LR = np.zeros((len(M2g), len(ag))); LP = np.zeros_like(LR); Lc = np.zeros_like(LR); medR = np.zeros_like(LR); meddmu = np.zeros_like(LR)
r0s = rng.choice(r0, NS)
pR0 = np.mean(r0 <= 1.10)
EP = [epochs() for _ in range(20)]
DES = []
for t, psi in EP:
    A3 = design(t, psi); m2 = t < 2016.38; A2 = A3[m2]
    DES.append((t, psi, A3, np.linalg.pinv(A3), m2, np.linalg.pinv(A2)))
for i, M2 in enumerate(M2g):
    eps = epsG(M2); q = M2 / M1; beta = q / (1 + q) - eps / (1 + eps)
    for j, a in enumerate(ag):
        P = np.sqrt(a ** 3 / (M1 + M2))
        alpha = beta * a * plx      # mas
        dchi = np.zeros(NS); dmu = np.zeros((NS, 2))
        for k in range(NS):
            t, psi, A3, A3p, m2, A2p = DES[k % 20]
            e = 0.0 if P * 365.25 < 10 else rng.uniform(0, 0.6)
            inc = np.arccos(rng.uniform(-1, 1)); om = rng.uniform(0, 2 * np.pi); Om = rng.uniform(0, 2 * np.pi); tp = rng.uniform(0, P)
            Mn = 2 * np.pi * (t - tp) / P
            E = kepler(Mn, e); X = np.cos(E) - e; Yv = np.sqrt(1 - e ** 2) * np.sin(E)
            # Thiele-Innes
            Ati = np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc); Bti = np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc)
            Fti = -np.sin(om) * np.cos(Om) - np.cos(om) * np.sin(Om) * np.cos(inc); Gti = -np.sin(om) * np.sin(Om) + np.cos(om) * np.cos(Om) * np.cos(inc)
            dd = alpha * (Ati * X + Fti * Yv); dra = alpha * (Bti * X + Gti * Yv)
            al = dra * np.sin(psi) + dd * np.cos(psi)
            p3 = A3p @ al; res = al - A3 @ p3
            dchi[k] = np.sum(res ** 2) / sig_tr ** 2
            p2 = A2p @ al[m2]
            dmu[k] = p3[2:4] - p2[2:4]
        ruwe = np.sqrt(r0s ** 2 + dchi / (NCCD - 5))
        LR[i, j] = (np.mean(ruwe <= 1.10) + 1e-4) / pR0
        num = np.exp(-0.5 * np.sum(((dmu_obs - dmu) / sig_cal) ** 2, axis=1)).mean()
        den = np.exp(-0.5 * np.sum((dmu_obs / sig_cal) ** 2))
        LP[i, j] = num / den
        # combined: per-sim joint
        jl = (ruwe <= 1.10) * np.exp(-0.5 * np.sum(((dmu_obs - dmu) / sig_cal) ** 2, axis=1))
        Lc[i, j] = (jl.mean() + 1e-6) / (pR0 * den)
        medR[i, j] = np.median(ruwe); meddmu[i, j] = np.median(np.hypot(*dmu.T))
    print(f'M2 {M2:.3f} (eps_G {eps:.3f}, beta {beta:.3f}): excluded a (combined L<0.05): ' +
          (f'{ag[Lc[i] < 0.05].min():.3f}-{ag[Lc[i] < 0.05].max():.2f} au' if (Lc[i] < 0.05).any() else 'none'), flush=True)
out = dict(M2=M2g.tolist(), a_au=ag.tolist(), L_ruwe=LR.tolist(), L_pm=LP.tolist(), L_combined=Lc.tolist(), median_ruwe=medR.tolist(),
           median_dmu=meddmu.tolist(), sig_tr=sig_tr, notes='L = likelihood relative to no companion; < 0.05 excluded')
json.dump(out, open(f'{B}/out/astrometric_limits.json', 'w'), indent=1)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
for ax, L, lab in zip(axs, [LR, LP, Lc], ['RUWE', 'DR3-DR2 PM', 'combined']):
    im = ax.pcolormesh(ag, M2g, np.log10(np.clip(L, 1e-4, 10)), vmin=-3, vmax=0.5, cmap='RdBu', shading='nearest')
    ax.contour(ag, M2g, L, levels=[0.05], colors='k')
    ax.set_xscale('log'); ax.set_xlabel('semi-major axis (au)'); ax.set_title(f'{lab}: log10 L (black: L = 0.05)')
axs[0].set_ylabel('companion mass (Msun)'); fig.colorbar(im, ax=axs, shrink=0.8)
fig.savefig(f'{B}/figs/astrometric_limits.png', dpi=110, bbox_inches='tight')
