"""b11: what follow-up closes each residual channel (quantitative design numbers).
(1) BD on the M3V (or on the primary of a binary): RV epochs. Posterior draws of (e, w, K) of the BD-on-M3V channel
    (b09, REF baseline) plus the BD prior (b03) for completeness. Epochs at fixed orbital phases relative to the transit
    (ephemeris known to ~2 min: phases are exact). Detection criterion: chi2 of a constant fit to the N RVs
    >= the 5-sigma (p = 5.7e-7) threshold for N - 1 dof, with per-epoch error sigma (instrumental + jitter 0.05 km/s).
    Reported: fraction of BD draws detected for sigma = 0.1, 0.3, 1, 3 km/s and phase sets of 2, 3, 4 epochs.
    Planet signal for comparison: K ~ 7 m/s per 10 M_E (referee_star) -> any variation > 0.3 km/s excludes a planet.
(2) BD on a bound companion: the companion carries a fraction f_c of the light (TESS 0.1-0.2) and moves with
    K2 = 15-35 km/s; the M3V is quiet. Requirement: detect a second, moving CCF component of relative strength
    ~f_c(band) -> high-resolution (R > 30,000) NIR spectra where f_c is larger (J/H: ~0.15-0.3); or two-band transit
    photometry (the depth ratio between blue and red bands differs from 1 when the host is a cooler star).
(3) PTP / BD-on-companion: high-resolution imaging (b09 what-ifs) and chromaticity (V/I depth ratio, b09 PTP_props).
Output: out/b11_followup.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from scipy.stats import chi2 as CHI2
rng = np.random.default_rng(5)
Z = np.load(f'{BD}/data/b09_bdm3v_post_draws.npz')
A = Z['A']; e, w, K = A[:, 0], A[:, 1], A[:, 6]
def rv(Kk, ee, ww, ph):
    # phase measured from mid-transit; transit at true anomaly f_tr = pi/2 - w
    f_tr = np.pi/2 - ww
    Etr = 2*np.arctan2(np.sqrt(1 - ee)*np.sin(f_tr/2), np.sqrt(1 + ee)*np.cos(f_tr/2)); Mtr = Etr - ee*np.sin(Etr)
    M = Mtr + 2*np.pi*ph
    E = M.copy()
    for _ in range(40): E = E - (E - ee*np.sin(E) - M)/(1 - ee*np.cos(E))
    f = 2*np.arctan2(np.sqrt(1 + ee)*np.sin(E/2), np.sqrt(1 - ee)*np.cos(E/2))
    return -Kk*(np.cos(f + ww) + ee*np.cos(ww))        # stellar RV (sign irrelevant for detection)
SETS = {'2 epochs (0.25, 0.75)': [0.25, 0.75], '2 epochs (0.15, 0.6)': [0.15, 0.6], '3 epochs (0.2, 0.5, 0.8)': [0.2, 0.5, 0.8],
        '4 epochs (0.15, 0.4, 0.65, 0.9)': [0.15, 0.4, 0.65, 0.9]}
out = {'BD_on_M3V_posterior': dict(K_kms_p05=float(np.percentile(K, 5)), K_kms_p50=float(np.median(K)), K_kms_p95=float(np.percentile(K, 95)),
                                     e_p05=float(np.percentile(e, 5)), e_p50=float(np.median(e)), e_p95=float(np.percentile(e, 95)))}
res = {}
for nm, phs in SETS.items():
    V = np.array([rv(K, e, w, p) for p in phs])           # [N, draws]
    row = {}
    for sig in (0.1, 0.3, 1.0, 3.0):
        s_ = np.sqrt(sig**2 + 0.05**2)
        Vn = V + rng.normal(0, s_, V.shape)
        c2 = np.sum((Vn - Vn.mean(0))**2, axis=0)/s_**2
        thr = CHI2.isf(5.7e-7, len(phs) - 1)
        row[f'sigma_{sig}_kms'] = float(np.mean(c2 > thr))
    res[nm] = row
    print(nm, row, flush=True)
out['rv_detection_fraction'] = res
# minimal K over the BD prior (13 M_J, e = 0) for reference
out['K_min_prior_kms'] = float((2*np.pi*G/(P*DAY))**(1/3.)*13*MJUP/((0.370 + 13/1047.57)*MSUN)**(2/3.)/1e5)
out['K_max_prior_kms'] = float((2*np.pi*G/(P*DAY))**(1/3.)*78.5*MJUP/((0.370 + 78.5/1047.57)*MSUN)**(2/3.)/1e5)
# companion-hosted BD: K2 for m2 = 0.1-0.3 Msun, M_bd 13-78 M_J
k2 = {}
for m2 in (0.1, 0.15, 0.2, 0.25, 0.3):
    k2[str(m2)] = {str(Mb): float((2*np.pi*G/(P*DAY))**(1/3.)*Mb*MJUP/((m2 + Mb/1047.57)*MSUN)**(2/3.)/1e5) for Mb in (13, 40, 70)}
out['K2_companion_kms'] = k2
print(out['K_min_prior_kms'], out['K_max_prior_kms'], k2, flush=True)
jdump(dict(method=__doc__, **out), f'{BD}/out/b11_followup.json')
