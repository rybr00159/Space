"""Activity, rotation and flares of TIC 214311908.
(1) X-ray: eRASS1 Main source 1eRASS J172352.9-511928 (Merloni+24; J/A+A/682/A34): flux, L_X, L_X/L_bol, positional
    association likelihoods with the M3V, the M-dwarf pair, the background star 3.1" NE and the Mira
    (per-axis sigma = POS_ERR = 5.38"; stars at the eRASS1 epoch 2020.2); chance-coincidence probability;
    Rossby number implied by the Wright+11/18 activity-rotation relation (tau from Wright+18 mass relation).
    GALEX: no source within 20" in GUVcat/AIS (coverage checked by counting GUVcat sources within 15').
(2) Rotation: the deblended M3V light curve of the photometry investigation (../phot/lc_M3V_all.csv,
    flux_raw_norm = 1 + dF/935, model B, NOT detrended), 2-h bins, transits masked. Per-sector Lomb-Scargle; the
    pair (model A raw) for comparison; sinusoid injection-recovery (P = 1-12 d) for the amplitude detectable in
    >= 90% of trials (peak within 10% of P_inj and above the highest peak of the real data); amplitude of a
    sinusoid at P_orb = 5.6915 d (synchronised spots) per sector.
(3) Flares: candidates (>= 3 consecutive points > 3 sigma, or 1 point > 6 sigma in 200-s data) in the detrended
    light curve; detection limit vs. amplitude.
Output: out/activity.json, figs/rotation.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.timeseries import LombScargle
B = 'tic214311908/analysis/star'
PH = 'tic214311908/analysis/phot'
sp = json.load(open(f'{B}/star_params.json'))
out = {}
# ---------------- X-ray
FX, eFX_rel = 5.4014e-14, 3.25 / 7.09   # erg/s/cm2 (0.2-2.3 keV), count-based relative error
d = sp['distance_pc']['median'] * 3.0857e18
LX = 4 * np.pi * d ** 2 * FX; Lbol = sp['luminosity_lsun']['median'] * 3.828e33
xr = dict(name='1eRASS J172352.9-511928', ra=260.9708, dec=-51.324522, pos_err_1sigma_per_axis=5.381, det_like=6.244, counts=[7.09, 3.25],
          rate_cts=[0.0580, 0.0266], flux_0p2_2p3keV=[FX, FX * eFX_rel], LX_erg_s=[LX, LX * eFX_rel], log_LX=float(np.log10(LX)),
          LX_over_Lbol_if_M3V=[LX / Lbol, LX / Lbol * eFX_rel], MJD=58934.0)
ep = 2000.0 + (58934.0 + 2400000.5 - 2451545.0) / 365.25
cands = {'M3V TIC 214311908': (260.969674, -51.323821, -113.671152, 85.54355), 'pair TIC 214311910 A(+B)': (260.974873, -51.324055, -112.796, 87.205),
         'background star 3.1" NE (G 17.5, 1.5 kpc)': (260.970239, -51.323025, 2.56, -0.98), 'Mira TIC 214311906': (260.9776, -51.32372, -5.96, -7.06)}
sig = 5.381; like = {}
for k, (ra, de, pa, pd_) in cands.items():
    dt = ep - 2016.0; ra_e = ra + pa * dt / 3.6e6 / np.cos(np.radians(de)); de_e = de + pd_ * dt / 3.6e6
    r = np.hypot((xr['ra'] - ra_e) * np.cos(np.radians(de)), xr['dec'] - de_e) * 3600
    like[k] = dict(sep_arcsec=float(r), r_over_sigma=float(r / sig), rel_likelihood=float(np.exp(-0.5 * (r / sig) ** 2)))
tot = sum(v['rel_likelihood'] for v in like.values())
for k in like: like[k]['posterior_equal_priors'] = like[k]['rel_likelihood'] / tot
xr['association'] = like
xr['chance_coincidence_within_10arcsec'] = float(9.3e5 / 20600 * np.pi * (10 / 3600) ** 2)
# activity-rotation (Wright+11: Rsat = 10^-3.13, Ro_sat = 0.13, beta = -2.18; tau Wright+18)
M = sp['mass_msun']['median']; tau = 10 ** (2.33 - 1.50 * M + 0.31 * M ** 2)
for lab, frac in (('all_from_M3V', 1.0), ('half_from_M3V', 0.5)):
    Rx = LX * frac / Lbol
    Ro = 0.13 * (Rx / 10 ** -3.13) ** (-1 / 2.18) if Rx < 10 ** -3.13 else 0.13
    xr[f'Ro_{lab}'] = float(Ro); xr[f'Prot_days_{lab}'] = float(Ro * tau)
xr['tau_conv_days'] = float(tau)
out['xray'] = xr
print(json.dumps({k: v for k, v in xr.items() if k != 'association'}, indent=0, default=float)); print(json.dumps(like, indent=0))
# GALEX coverage
try:
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    c0 = SkyCoord(260.969674, -51.323821, unit='deg')
    n15 = {}
    for cat in ('II/335/gufuv', 'II/312/ais'):
        r = Vizier(row_limit=5000, columns=['_r']).query_region(c0, radius=15 * u.arcmin, catalog=cat)
        n15[cat] = int(len(r[0])) if len(r) else 0
    out['galex_sources_within_15arcmin'] = n15; print('GALEX within 15 arcmin:', n15)
except Exception as e:
    out['galex_sources_within_15arcmin'] = str(e)
# ---------------- rotation
lc = pd.read_csv(f'{PH}/lc_M3V_all.csv'); pr = pd.read_csv(f'{PH}/lc_pair_all.csv')
rng = np.random.default_rng(4)
P_orb, T0 = 5.6914968, 3840.9153
def prep(df, col):
    df = df[(df.quality_ok == 1) & np.isfinite(df[col])].copy()
    ph = ((df.time_btjd_tdb - T0) / P_orb + 0.5) % 1 - 0.5
    df = df[np.abs(ph * P_orb) > 0.08]
    return df
rot = {}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 3, figsize=(14, 6.5)); axs = axs.ravel()
for ia, sec in enumerate(sorted(lc.sector.unique())):
    a = prep(lc[lc.sector == sec], 'flux_raw_norm'); b = prep(pr[pr.sector == sec], 'flux_raw_norm_modelA')
    def bin2h(df, col):
        t = df.time_btjd_tdb.values; f = df[col].values; e = np.floor((t - t.min()) / (2 / 24)).astype(int)
        g = pd.DataFrame(dict(e=e, t=t, f=f)).groupby('e').agg(t=('t', 'mean'), f=('f', 'median'), n=('f', 'size'))
        g = g[g.n >= 3]; return g.t.values, g.f.values
    t, f = bin2h(a, 'flux_raw_norm'); tb, fb = bin2h(b, 'flux_raw_norm_modelA')
    f = f - np.median(f); fb = fb - np.median(fb)
    span = t.max() - t.min(); fmin = 1 / min(15, span * 0.8); freq = np.linspace(fmin, 1 / 0.5, 3000)
    pw = LombScargle(t, f).power(freq); pwb = LombScargle(tb, fb).power(freq)
    best = 1 / freq[np.argmax(pw)]; bestb = 1 / freq[np.argmax(pwb)]
    # amplitude of the best sinusoid
    def sinamp(t, f, P):
        A = np.column_stack([np.ones(len(t)), np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P), t - t.mean()])
        c, *_ = np.linalg.lstsq(A, f, rcond=None); r = f - A @ c
        cov = np.linalg.inv(A.T @ A) * np.var(r)
        amp = np.hypot(c[1], c[2]); e = np.sqrt((c[1] ** 2 * cov[1, 1] + c[2] ** 2 * cov[2, 2]) / amp ** 2)
        return float(amp), float(e)
    a_best = sinamp(t, f, best); a_orb = sinamp(t, f, P_orb); a_orb2 = sinamp(t, f, P_orb / 2)
    rms = float(1.4826 * np.median(np.abs(f - np.median(f))))
    # injection-recovery
    thr = pw.max(); det = {}
    for Pinj in (1.0, 2.0, 3.0, 4.0, 5.6915, 7.0, 9.0, 12.0):
        if Pinj > 0.8 * span: continue
        lim = None
        for A in (0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05):
            ok = 0
            for k in range(20):
                ph0 = rng.uniform(0, 2 * np.pi); fi = f + A * np.sin(2 * np.pi * t / Pinj + ph0)
                pi = LombScargle(t, fi).power(freq); Pb = 1 / freq[np.argmax(pi)]
                ok += (abs(Pb - Pinj) / Pinj < 0.1) and (pi.max() > thr)
            if ok >= 18: lim = A; break
        det[f'{Pinj}'] = lim
    rot[int(sec)] = dict(span_d=float(span), rms_2h=rms, best_period_M3V=float(best), best_power=float(pw.max()), amp_best=a_best,
                         best_period_pair=float(bestb), amp_at_Porb=a_orb, amp_at_Porb_half=a_orb2, inj_semi_amp_90pct=det)
    print(f'S{sec}: rms(2h) {rms*1e3:.1f} ppt; best P {best:.2f} d (amp {a_best[0]*1e3:.1f} ppt), pair best {bestb:.2f} d; amp at P_orb {a_orb[0]*1e3:.2f}+-{a_orb[1]*1e3:.2f} ppt; 90% limits {det}')
    ax = axs[ia]; ax.plot(1 / freq, pw, 'C0', lw=0.8, label='M3V'); ax.plot(1 / freq, pwb, 'C1', lw=0.6, alpha=0.7, label='pair (model A)')
    ax.axvline(P_orb, color='C3', ls=':'); ax.set_xscale('log'); ax.set_title(f'S{sec}: LS of deblended raw flux (2-h bins)', fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{B}/figs/rotation.png', dpi=110)
out['rotation'] = rot
# combined amplitude at P_orb (weighted)
am = np.array([[v['amp_at_Porb'][0], v['amp_at_Porb'][1]] for v in rot.values()])
out['rotation_amp_at_Porb_weighted'] = [float(np.sum(am[:, 0] / am[:, 1] ** 2) / np.sum(1 / am[:, 1] ** 2)), float(1 / np.sqrt(np.sum(1 / am[:, 1] ** 2)))]
# ---------------- flares
fl = {}
for sec in sorted(lc.sector.unique()):
    a = lc[(lc.sector == sec) & (lc.quality_ok == 1) & np.isfinite(lc.flux_norm)].sort_values('time_btjd_tdb')
    f = a.flux_norm.values - 1; s = 1.4826 * np.median(np.abs(f - np.median(f))); z = f / s
    run = np.convolve((z > 3).astype(int), np.ones(3, int), 'same') >= 3
    cand3 = int(np.sum(np.diff(np.concatenate([[0], run.astype(int)])) == 1))
    cand6 = int(np.sum(z > 6))
    dur = float(len(a) * a.cadence_s.median() / 86400)
    fl[int(sec)] = dict(days=dur, sigma_ppt=float(s * 1e3), n_3x3sigma_runs=cand3, n_single_gt6sigma=cand6, cadence_s=float(a.cadence_s.median()))
out['flares'] = fl
tot_days = sum(v['days'] for v in fl.values())
out['flare_summary'] = dict(total_days=tot_days, n_candidates=sum(v['n_3x3sigma_runs'] for v in fl.values()),
                            note='phot/ found 5 unconvincing single-point candidates; per-cadence 5-sigma thresholds 40-90 ppt')
print(json.dumps(fl, indent=0)); print('weighted amplitude at P_orb:', out['rotation_amp_at_Porb_weighted'])
json.dump(out, open(f'{B}/out/activity.json', 'w'), indent=1, default=float)
