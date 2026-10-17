"""f03: the Mira (TIC 214311906, Gaia DR3 5925373268528684288): brightness range per band and pulsation phase at the
2027-2028 transit epochs.
Data (read-only): Gaia DR3 epoch photometry G/BP/RP (39 FoV transits 2014.6-2017.4; ctx/lc/gaia_epoch), Gaia DR3
vari_summary/vari_long_period_variable (data/ from f02a), TESS-band Mira flux and slope at each transit (clean/out/
per_transit_correlates.csv, 21 epochs, S12-S104 = 2019-2026), SkyMapper DR4 max-min per band, DENIS 1998 vs 2MASS 1999,
VSX (V 12.91-17.01, P 294 d).
Phase model: a sinusoid with common period P and epoch of maximum t0 fitted jointly to Gaia G magnitudes and TESS
fluxes (+ TESS in-sector slopes as derivative constraints), separate mean/amplitude per band; Mira light curves are
not sinusoidal, so an extra scatter term (0.15 mag in G, 8% in TESS flux, 30% on slopes) is included and the period
posterior comes from the profile chi2 over a 230-380 d grid. Phase predictions at the observable transits of
out/transits_obs_all_sites.csv carry the P/t0 uncertainty plus an ASSUMED random cycle-to-cycle phase jitter of
0.05 per cycle (random walk) - typical Mira period noise, not measured here.
Outputs: out/mira.json, out/mira_phase_at_transits.csv, figs/mira_fit.png"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json
from scipy.optimize import least_squares
g = pd.read_csv(f'{B}/ctx/lc/gaia_epoch/mira_dr3_epoch_G.csv')
g = g[~g.rejected_by_photometry & np.isfinite(g.g_transit_mag)]
tg = g.g_transit_time.values + 2455197.5          # Gaia epoch times: TCB BJD - 2455197.5
mg = g.g_transit_mag.values
gb = g[np.isfinite(g.bp_mag) & np.isfinite(g.rp_mag)]
c = pd.read_csv(f'{B}/clean/out/per_transit_correlates.csv'); c = c[c.pipeline == 'ref'].drop_duplicates('epoch')
tt = c.tc.values + 2457000.0; ft = c.mira_flux.values; st = c.mira_slope_per_d.values
def model(p, P, t_g, t_t):
    t0, aG, AG, F0, F1 = p
    ph_g = 2*np.pi*(t_g - t0)/P; ph_t = 2*np.pi*(t_t - t0)/P
    return aG - AG*np.cos(ph_g), F0 + F1*np.cos(ph_t), -F1*2*np.pi/P*np.sin(ph_t)
sG, sF, sS = 0.15, 0.08, 0.30
def resid(p, P):
    mG, fT, dT = model(p, P, tg, tt)
    return np.concatenate([(mg - mG)/sG, (ft - fT)/(sF*ft), (st - dT)/(sS*np.abs(st) + 20.0)])
Ps = np.arange(230, 380.01, 0.25)
best = []
for P in Ps:
    bestc = None
    for t0s in np.linspace(0, P, 12, endpoint=False):
        p0 = [2457000 + t0s, mg.mean(), 0.8, ft.mean(), 3000.]
        r = least_squares(resid, p0, args=(P,), bounds=([2456500, 5, 0, 1000, 0], [2458500, 20, 3, 1e5, 2e4]))
        if bestc is None or r.cost < bestc.cost: bestc = r
    best.append((P, 2*bestc.cost, bestc.x))
chi = np.array([b[1] for b in best]); i0 = int(np.argmin(chi)); Pbest, pbest = best[i0][0], best[i0][2]
w = np.exp(-0.5*(chi - chi.min())); w /= w.sum()
Pmean = float(np.sum(w*Ps)); Psd = float(np.sqrt(np.sum(w*(Ps - Pmean)**2)))
peaks = sorted([(float(Ps[i]), float(chi[i] - chi.min())) for i in range(1, len(Ps) - 1) if chi[i] < chi[i - 1] and chi[i] < chi[i + 1]], key=lambda x: x[1])[:5]
print('best P', Pbest, 'chi2', chi.min(), 'n', len(tg) + 2*len(tt), 'P mean/sd', Pmean, Psd, 'local minima', peaks)
# t0 uncertainty at best P via Jacobian
r = least_squares(resid, pbest, args=(Pbest,))
J = r.jac; cov = np.linalg.inv(J.T @ J)*max(1.0, 2*r.cost/(len(r.fun) - 5))
t0, aG, AG, F0, F1 = r.x; t0_sd = float(np.sqrt(cov[0, 0]))
# tie t0 to the latest data: epoch of maximum nearest to 2026.3
nlast = np.round((tt.max() - t0)/Pbest); tmax_last = t0 + nlast*Pbest
# ---- brightness ranges by band
lpv = pd.read_csv(f'{F}/data/gaia_mira_lpv.csv'); vs = pd.read_csv(f'{F}/data/gaia_mira_vari_summary.csv')
Tmag = lambda f: 20.44 - 2.5*np.log10(f)         # TESS zero point (TESS Instrument Handbook: 15000 e-/s for T=10)
TmaxF = F0 + F1; TminF = max(F0 - F1, 100.)
ranges = {
 'V (Vega)': dict(bright=12.91, faint=17.01, source='VSX (ASAS-SN) max/min'),
 'BP (Vega)': dict(bright=float(vs.min_mag_bp[0]), faint=float(vs.max_mag_bp[0]), source='Gaia DR3 vari_summary, 39 epochs 2014-17'),
 'G (Vega)': dict(bright=float(vs.min_mag_g_fov[0]), faint=float(vs.max_mag_g_fov[0]), source='Gaia DR3 vari_summary'),
 'RP (Vega)': dict(bright=float(vs.min_mag_rp[0]), faint=float(vs.max_mag_rp[0]), source='Gaia DR3 vari_summary'),
 'T (TESS)': dict(bright=float(Tmag(TmaxF)), faint=float(Tmag(TminF)), source=f'sinusoid fit to 21 TESS epochs (observed {Tmag(ft.max()):.2f}-{Tmag(ft.min()):.2f}); zp 20.44'),
 'g (AB)': dict(bright=float(vs.min_mag_bp[0]) - 0.13, faint=float(vs.max_mag_bp[0]) - 0.13, source='BP range shifted by SkyMapper g mean 15.12 - Gaia BP mean 15.25; SkyMapper g max-min >= 2.98 (8 epochs)'),
 'r (AB)': dict(bright=12.6, faint=15.6, source='SkyMapper r mean 13.95 +- 1.5 (amplitude interpolated between BP 3.7 and G 2.0 mag); assumption'),
 'i (AB)': dict(bright=float(vs.min_mag_rp[0]) + 0.34, faint=float(vs.max_mag_rp[0]) + 0.34, source='RP range shifted by SkyMapper i mean 11.14 - RP mean 10.80; SkyMapper i max-min >= 1.14 (7 epochs)'),
 'z (AB)': dict(bright=9.3, faint=10.6, source='SkyMapper z mean 9.93; amplitude ~1.3 assumed (between RP 1.74 and J >=0.5; SkyMapper z max-min >= 0.66, 4 epochs)'),
 'Ic (Vega)': dict(bright=float(vs.min_mag_rp[0]) + 0.34 - 0.45, faint=float(vs.max_mag_rp[0]) + 0.34 - 0.45, source='i(AB) range - 0.45 (AB->Vega for I); +-0.3 systematic'),
 'J (Vega)': dict(bright=7.6, faint=8.8, source='2MASS 1999 8.08, DENIS 1998 8.56 (>= 0.48 apart); amplitude ~1.2 assumed'),
 'Ks (Vega)': dict(bright=6.3, faint=7.3, source='2MASS 1999 6.68, DENIS 1998 7.21 (>= 0.53 apart); amplitude ~1.0 assumed'),
}
# ---- phase at the observable transits
O = pd.read_csv(f'{OUT}/transits_obs_all_sites.csv')
O = O[O.full_window_sun18 | O.good_partial].drop_duplicates('epoch')
rng = np.random.default_rng(3)
Pdraw = rng.choice(Ps, size=4000, p=w)
out_rows = []
for row in O.itertuples():
    tb = row.bjd_tdb
    ncyc = (tb - tmax_last)/Pbest
    phs = []
    for Pd in Pdraw[:2000]:
        # anchor on the last observed maximum (t0 error) and propagate with the drawn period plus phase jitter
        tl = tmax_last + rng.normal(0, t0_sd)
        ph = ((tb - tl)/Pd + rng.normal(0, 0.05*np.sqrt(abs((tb - tl)/Pd)))) % 1
        phs.append(ph)
    phs = np.array(phs)
    # circular mean / sd
    z = np.mean(np.exp(2j*np.pi*phs)); phm = (np.angle(z)/(2*np.pi)) % 1; phsd = float(np.sqrt(-2*np.log(abs(z)))/(2*np.pi))
    Gp = aG - AG*np.cos(2*np.pi*phs); Fp = F0 + F1*np.cos(2*np.pi*phs)
    out_rows.append(dict(epoch=row.epoch, utc_mid=row.utc_mid, mira_phase_from_max=round(float(phm), 3), phase_sd=round(phsd, 3),
                         G_pred_median=round(float(np.median(Gp)), 2), G_pred_16_84=[round(float(x), 2) for x in np.percentile(Gp, [16, 84])],
                         T_pred_median=round(float(Tmag(np.median(Fp))), 2), frac_of_range_bright=round(float(np.mean(np.cos(2*np.pi*phs) > 0.5)), 2)))
R = pd.DataFrame(out_rows); R.to_csv(f'{OUT}/mira_phase_at_transits.csv', index=False)
print(R.to_string())
jdump(dict(period_fit=dict(P_best_d=Pbest, P_mean_d=Pmean, P_sd_d=Psd, local_minima_P_dchi2=peaks, chi2_min=float(chi.min()),
                           n_data=int(len(tg) + 2*len(tt)), t_max_epoch0_BJD=float(t0), t0_sd_d=t0_sd, last_max_before_2026_5_BJD=float(tmax_last),
                           G_mean=float(aG), G_semi_amp=float(AG), T_flux_mean=float(F0), T_flux_semi_amp=float(F1),
                           note='sinusoid; Mira light curves are asymmetric; extra scatter 0.15 mag (G), 8% (TESS flux), 30% (slopes)'),
           gaia_lpv=dict(frequency=float(lpv.frequency[0]), period_d=float(1/lpv.frequency[0]), period_err_d=float(lpv.frequency_error[0]/lpv.frequency[0]**2),
                         amplitude=float(lpv.amplitude[0])),
           vsx=dict(V_max=12.91, V_min=17.01, P_d=294.0),
           ranges=ranges, phase_jitter_assumption='0.05 per cycle random walk (assumed)',
           tess_observed_flux_range=[float(ft.min()), float(ft.max())]), f'{OUT}/mira.json')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
tm = np.linspace(min(tg.min(), tt.min()) - 100, 2462000, 3000)
mG, fT, _ = model(r.x, Pbest, tm, tm)
ax[0].plot(tg - 2450000, mg, 'k.', label='Gaia DR3 G epochs'); ax[0].plot(tm - 2450000, mG, 'C0', lw=1); ax[0].invert_yaxis(); ax[0].set_ylabel('G')
ax[0].legend(fontsize=8)
ax[1].plot(tt - 2450000, ft, 'r.', label='TESS Mira flux at transits'); ax[1].plot(tm - 2450000, fT, 'C1', lw=1); ax[1].set_ylabel('TESS flux (e-/s)')
for row in R.itertuples():
    ax[1].axvline(O[O.epoch == row.epoch].bjd_tdb.iloc[0] - 2450000, color='0.7', lw=0.5)
ax[1].set_xlabel('BJD - 2450000'); ax[1].legend(fontsize=8); ax[0].set_title(f'Mira sinusoid fit: P = {Pbest:.1f} d (profile mean {Pmean:.1f} +- {Psd:.1f})')
plt.tight_layout(); plt.savefig(f'{F}/figs/mira_fit.png', dpi=110)
