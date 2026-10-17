"""f09: prioritised observing plan, FPP ladder and figures; collects the key numbers of f01-f08.
* PTP rejection power of the two-band test per number of events: per-draw ratio (f07b method) vs per-event sigma_R
  from out/depth_ratio_tests.csv; fraction of the PTP posterior rejected at > 3 sigma.
* FPP ladder: out/fpp_after_followup.csv; for two-band steps that include g or V the WD-occultation term is removed
  (g/TESS depth ratio ~6, f04), which r05's 'chroma' option does not do.
Outputs: out/plan.json, out/fpp_ladder.csv, out/two_band_ptp_rejection.csv, figs/observability.png, figs/hri_ptp.png, figs/fpp_ladder.png"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S1, S2, S3, S4, INK, INK2, SURF = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#0b0b0b', '#52514e', '#fcfcfb'
plt.rcParams.update({'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF, 'axes.facecolor': SURF})
# ---------------- two-band PTP rejection ----------------
W = np.load(f'{OUT}/ptp_draw_weights.npz'); C = pd.read_csv(f'{OUT}/chromatic_ratio_vs_dT.csv').sort_values('dT')
T = pd.read_csv(f'{OUT}/depth_ratio_tests.csv')
dT = W['dT'].astype(float)
rows = []
for case in ('ref_beta__signfixed', 'inv_beta__signfixed', 'ref_beta__as_published'):
    w = W[case].astype(float); ok = (w > 0) & (dT >= C.dT.min()) & (dT <= C.dT.max())
    for t in T.itertuples():
        a, b = t.pair.split('/'); a = {'z': 'z_s'}.get(a, a); b = {'z': 'z_s'}.get(b, b)
        r = np.interp(dT[ok], C.dT, C[f'ratio_{a}'])/np.interp(dT[ok], C.dT, C[f'ratio_{b}'])
        for n in (1, 2, 3, 4, 6, 8):
            nev = n                                    # events per band (simultaneous: n events; alternating: 2n events)
            s = t.sigma_ratio_per_event/np.sqrt(nev)
            rej = float(np.sum(w[ok]*(np.abs(r - 1) > 3*s))/w[ok].sum())
            rows.append(dict(case=case, telescope=t.telescope, pair=t.pair, simultaneous=t.simultaneous, events_per_band=n,
                             total_events=n if t.simultaneous else 2*n, sigma_R=round(s, 4), frac_PTP_rejected_3sigma=round(rej, 3)))
TB = pd.DataFrame(rows); TB.to_csv(f'{OUT}/two_band_ptp_rejection.csv', index=False)
print(TB[(TB.case == 'ref_beta__signfixed') & TB.events_per_band.isin([1, 2, 4])].pivot_table(index=['telescope', 'pair'], columns='events_per_band', values='frac_PTP_rejected_3sigma').to_string())
# ---------------- FPP ladder ----------------
FP = pd.read_csv(f'{OUT}/fpp_after_followup.csv')
FP['FPP_adj'] = np.where(FP.step.str.contains('two-band') & ~FP.step.str.contains('RV'), FP.FPP - FP.WD, FP.FPP)
FP['note'] = np.where(FP.step.str.contains('two-band') & ~FP.step.str.contains('RV'), 'WD term removed (g or V band in the two-band test)', '')
FP.to_csv(f'{OUT}/fpp_ladder.csv', index=False)
lad = FP.pivot_table(index='step', columns=['sign', 'lc'], values='FPP_adj')
print(lad.map(lambda x: f'{x:.1e}').to_string())
# ---------------- figures ----------------
# (1) observability timeline
O = pd.read_csv(f'{OUT}/transits_obs_all_sites.csv')
sites = list(SITES.keys())
fig, ax = plt.subplots(figsize=(10, 3.8))
for i, s in enumerate(sites):
    o = O[O.site == s]
    full = o[o.full_window_sun18]; part = o[o.good_partial]
    tf = pd.to_datetime(full.utc_mid); tp = pd.to_datetime(part.utc_mid)
    ax.scatter(tf, [i]*len(tf), s=42, color=S1, edgecolor=SURF, linewidth=1.5, zorder=3, label='full window (Sun < -18, airmass < 2)' if i == 0 else None)
    ax.scatter(tp, [i]*len(tp), s=42, facecolor='none', edgecolor=S1, linewidth=1.2, zorder=3, label='transit + >= 45 min baseline each side' if i == 0 else None)
ax.set_yticks(range(len(sites))); ax.set_yticklabels([s.split(' (')[0] for s in sites])
ax.grid(axis='x', color='#e6e5e0', lw=0.6); ax.set_axisbelow(True)
ax.set_title('Observable transits of TIC 214311908.01, 2026-11 to 2028-12 (window 4.9-5.0 h)', color=INK, loc='left')
ax.legend(loc='upper left', fontsize=8, frameon=False, ncol=2, bbox_to_anchor=(0, -0.08))
plt.tight_layout(); plt.savefig(f'{F}/figs/observability.png', dpi=130); plt.close()
# (2) HRI: PTP companions vs contrast curves (I band and Ks band, separate panels)
H = json.load(open(f'{OUT}/hri.json'))['curves']
w = W['ref_beta__signfixed'].astype(float); sep = W['sep'].astype(float)
fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
for ax, col, band, names in ((axs[0], 'dI', 'I', ['Zorro 832 nm (nominal)', 'SOAR HRCam I (16th percentile)', 'SOAR HRCam I (median of 2684 TOIs)']),
                             (axs[1], 'dKs', 'Ks', ['VLT/ERIS-NIX Ks LGS (nominal)', 'VLT/ERIS-NIX Ks LGS (pessimistic)', 'VLT/SPHERE IRDIS (degraded, G=14.3)'])):
    y = W[col].astype(float); m = (w > 0) & (sep > 0.005) & (sep < 5)
    hb = ax.hexbin(np.log10(sep[m]), y[m], C=w[m], reduce_C_function=np.sum, gridsize=40, cmap='Greys', mincnt=1, extent=(-2, 0.7, 0, 5), bins='log')
    for nm, cc in zip(names, (S1, S2, S3)):
        pts = np.array(H[nm]['points']); ax.plot(np.log10(pts[:, 0]), pts[:, 1], color=cc, lw=2, label=nm)
    ax.set_xlabel('log10 projected separation (arcsec)'); ax.set_ylabel(f'companion contrast d{band} (mag)')
    ax.set_ylim(5, 0); ax.set_xlim(-2, 0.7); ax.legend(fontsize=7, frameon=False, loc='lower right')
axs[0].set_title('Planet-on-companion (PTP) posterior vs 5-sigma curves, I band', loc='left', color=INK, fontsize=9)
axs[1].set_title('same, Ks band (AO)', loc='left', color=INK, fontsize=9)
plt.tight_layout(); plt.savefig(f'{F}/figs/hri_ptp.png', dpi=130); plt.close()
# (3) FPP ladder
steps = ['0 none', '1a speckle Zorro 832 (nominal)', '1b speckle SOAR HRCam (16th pct)', '1c AO ERIS Ks (nominal)', '2b two-band 5% + pair resolved',
         '3 Zorro + two-band 5%', '4 Zorro + two-band 5% + low-precision RV (BD/WD channels closed)', '6 + precise RV mass on the M3V (K detected at >= 5 sigma)']
fig, ax = plt.subplots(figsize=(9, 4))
for (lc, sg), cc, mk, lab in ((('ref_beta', 'fixed'), S1, 'o', 'REF light curve, sign-fixed'), (('inv_beta', 'fixed'), S2, 's', 'INV light curve, sign-fixed'),
                               (('ref_beta', 'as_published'), S3, '^', 'REF, as published (referee2)')):
    v = [float(FP[(FP.lc == lc) & (FP.sign == sg) & (FP.step == st)].FPP_adj.iloc[0]) for st in steps]
    ax.plot(v, range(len(steps)), marker=mk, ms=8, color=cc, lw=0, label=lab, markeredgecolor=SURF)
ax.axvline(0.015, color=INK2, ls=':', lw=1); ax.text(0.0105, -0.35, 'FPP 1.5%', color=INK2, fontsize=8, ha='right')
ax.set_xscale('log'); ax.set_yticks(range(len(steps))); ax.set_yticklabels([s[:60] for s in steps], fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel('false-positive probability after the step (r05 framework; steps 1a-1c are alternatives)'); ax.legend(fontsize=7.5, frameon=False, loc='upper left')
ax.grid(axis='x', color='#e6e5e0', lw=0.6); ax.set_axisbelow(True)
ax.set_title('Expected FPP after each follow-up step', loc='left', color=INK)
plt.tight_layout(); plt.savefig(f'{F}/figs/fpp_ladder.png', dpi=130); plt.close()
# ---------------- plan summary ----------------
f01 = json.load(open(f'{OUT}/f01_transits.json')); prec = pd.read_csv(f'{OUT}/precision.csv'); hx = pd.read_csv(f'{OUT}/hri_exclusion.csv')
rv = json.load(open(f'{OUT}/f08_rv_atmo.json')); con = json.load(open(f'{OUT}/contamination_summary.json')); mira = json.load(open(f'{OUT}/mira.json'))
def fpp(step, lc='ref_beta', sg='fixed'): return float(FP[(FP.lc == lc) & (FP.sign == sg) & (FP.step == step)].FPP_adj.iloc[0])
plan = [
 dict(priority=1, step='Seeing-limited multi-band transit on the M3V (MuSCAT4 g r i zs at FTS; or LCOGT 1m g/i/zs at CTIO/SAAO/SSO, SPECULOOS/TRAPPIST I+z)',
      tests='signal is real and on the M3V (not a Mira/scene systematic); pair A, 3.1" and 5.5" stars resolved; achromaticity (PTP, WD, background EB)',
      events='1-2 (MuSCAT4), 3-5 (1m alternating bands)', FPP_after_REF_INV_signfixed=[fpp('2b two-band 5% + pair resolved'), fpp('2b two-band 5% + pair resolved', 'inv_beta')]),
 dict(priority=2, step='Speckle imaging Gemini-S/Zorro 832+562 nm (or SOAR/HRCam as a cheaper first look)',
      tests='bound companions >= 0.02" (PTP 97% excluded with Zorro; 63% with HRCam)', events='1 visit',
      FPP_after_REF_INV_signfixed=[fpp('1a speckle Zorro 832 (nominal)'), fpp('1a speckle Zorro 832 (nominal)', 'inv_beta')],
      FPP_with_step1=[fpp('3 Zorro + two-band 5%'), fpp('3 Zorro + two-band 5%', 'inv_beta')]),
 dict(priority=3, step='2-4 low-precision RV epochs near quadrature (sigma <= 0.5 km/s; FEROS, GHOST, SALT-HRS, HARPS/NIRPS)',
      tests='BD (K 2.8-15.6 km/s) and WD (K ~73 km/s) channels; SB2', events='2-4 spectra',
      FPP_after_REF_INV_signfixed=[fpp('4 Zorro + two-band 5% + low-precision RV (BD/WD channels closed)'), fpp('4 Zorro + two-band 5% + low-precision RV (BD/WD channels closed)', 'inv_beta')]),
 dict(priority=4, step='Precise RV mass (ESPRESSO; or NIRPS+HARPS in Jun-Nov when |Vsys-BERV| > 10 km/s)',
      tests='planet mass, K ~15 m/s; confirms the host', events='9-26 ESPRESSO epochs (jitter 5-10 m/s); 16-33 NIRPS+HARPS',
      FPP_after_REF_INV_signfixed=[fpp('6 + precise RV mass on the M3V (K detected at >= 5 sigma)'), fpp('6 + precise RV mass on the M3V (K detected at >= 5 sigma)', 'inv_beta')]),
]
jdump(dict(plan=plan, fpp_now=dict(REF_signfixed=fpp('0 none'), INV_signfixed=fpp('0 none', 'inv_beta'), REF_published=fpp('0 none', 'ref_beta', 'as_published'),
                                   INV_published=fpp('0 none', 'inv_beta', 'as_published')),
           two_band_rejection=TB[TB.case == 'ref_beta__signfixed'].to_dict('records')), f'{OUT}/plan.json')
for p in plan: print(p)
