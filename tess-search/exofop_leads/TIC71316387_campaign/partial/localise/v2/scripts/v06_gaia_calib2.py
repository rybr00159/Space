"""Extended Gaia DR3 IPD calibration with the SOAR speckle sample + field-star base rates.

(a) IPDfmp in bins of (sep, dmag_I) for Gaia-unresolved speckle binaries, G 9-14;
(b) BP blended-transit fraction (NBPblend / o_BPmag) vs separation, binaries vs singles;
(c) ipd_gof_harmonic_phase vs speckle PA (mod 180) for binaries with IPDgofha > 0.02 and 0.2 < sep < 2";
(d) RUWE vs sep;
(e) systems resembling the target: G 10.5-13, RUWE < 1.4, IPDfmp 25-60, blend fraction >= 0.8;
(f) field stars (3 deg, G 10.75-12.75): fraction with IPDfmp >= 30 and with blend fraction >= 0.9, split by
    nearest Gaia neighbour distance.
"""
import numpy as np, pandas as pd, json
from astropy.table import Table
from scipy.spatial import cKDTree
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/speckle'
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
S = pd.read_csv(f'{O}/calib_sample.csv')
S['bfrac'] = S.NBPblend / S.o_BPmag
S['rfrac'] = S.NRPblend / S.o_RPmag
S['dm'] = np.where(np.isfinite(S.dmI), S.dmI, S.dmag)
out = {}
B = S[(S.kind == 'binary') & (S.Gmag > 9) & (S.Gmag < 14) & ~S.gaia_resolved.astype(bool)].copy()
C = S[(S.kind == 'single') & (S.Gmag > 9) & (S.Gmag < 14) & (S.dm1 >= 3.5)].copy()
print(f'G 9-14: Gaia-unresolved binaries {len(B)}, single controls {len(C)}')
# (a)
sepb = [0.03, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.3, 1.8, 3.0]
dmb = [0, 1.5, 2.5, 3.5, 5, 8]
print('\n(a) median IPDfmp [n] by sep (rows) and dm_I (cols)', dmb)
tabA = []
for i in range(len(sepb) - 1):
    line = f'  sep {sepb[i]:.2f}-{sepb[i+1]:.2f}: '
    for j in range(len(dmb) - 1):
        s = B[(B.sep >= sepb[i]) & (B.sep < sepb[i + 1]) & (B.dm >= dmb[j]) & (B.dm < dmb[j + 1])]
        if len(s):
            line += f' {np.median(s.IPDfmp):5.0f}[{len(s):3d}]'
            tabA.append(dict(sep=[sepb[i], sepb[i + 1]], dm=[dmb[j], dmb[j + 1]], n=int(len(s)), fmp_med=float(np.median(s.IPDfmp)),
                             fmp_q16=float(np.percentile(s.IPDfmp, 16)), fmp_q84=float(np.percentile(s.IPDfmp, 84)),
                             ruwe_med=float(np.median(s.RUWE)), bfrac_med=float(np.nanmedian(s.bfrac))))
        else:
            line += '     -     '
    print(line)
out['fmp_grid'] = tabA
# (b)
print('\n(b) BP blended fraction vs sep (binaries, G 9-14):')
tabB = []
for i in range(len(sepb) - 1):
    s = B[(B.sep >= sepb[i]) & (B.sep < sepb[i + 1])]
    if len(s):
        f90 = float((s.bfrac >= 0.9).mean())
        print(f'  sep {sepb[i]:.2f}-{sepb[i+1]:.2f}: n={len(s):3d} median bfrac {np.nanmedian(s.bfrac):.2f}  frac(bfrac>=0.9) {f90:.2f}  median IPDfmp {np.median(s.IPDfmp):.0f} median RUWE {np.median(s.RUWE):.2f}')
        tabB.append(dict(sep=[sepb[i], sepb[i + 1]], n=int(len(s)), bfrac_med=float(np.nanmedian(s.bfrac)), frac_bfrac_ge0p9=f90,
                         ruwe_med=float(np.median(s.RUWE)), fmp_med=float(np.median(s.IPDfmp))))
print(f'  singles: n={len(C)} median bfrac {np.nanmedian(C.bfrac):.2f} frac(bfrac>=0.9) {(C.bfrac >= 0.9).mean():.3f}  frac(IPDfmp>=30) {(C.IPDfmp >= 30).mean():.4f}')
out['blend_vs_sep'] = tabB
out['singles_G9_14'] = dict(n=int(len(C)), frac_bfrac_ge0p9=float((C.bfrac >= 0.9).mean()), frac_fmp_ge30=float((C.IPDfmp >= 30).mean()),
                            n_fmp_ge30=int((C.IPDfmp >= 30).sum()), n_both=int(((C.IPDfmp >= 30) & (C.bfrac >= 0.9)).sum()))
# (c) harmonic phase vs PA
H = B[(B.IPDgofha > 0.02) & (B.sep > 0.2) & (B.sep < 2.0) & np.isfinite(B.pa)].copy()
d = ((H.IPDgofhp - (H.pa % 180)) + 90) % 180 - 90
H['dphi'] = d
print(f'\n(c) harmonic phase - PA (mod 180) for {len(H)} binaries with IPDgofha>0.02: median {np.median(d):+.1f}, '
      f'circular-ish spread (MAD) {1.4826*np.median(np.abs(d-np.median(d))):.1f} deg')
hist, edges = np.histogram(d, bins=np.arange(-90, 91, 15))
print('   histogram (15 deg bins from -90):', hist.tolist())
d2 = ((H.IPDgofhp - (H.pa % 180) - 90) + 90) % 180 - 90
out['phase_minus_pa'] = dict(n=int(len(H)), median=float(np.median(d)), mad_sd=float(1.4826 * np.median(np.abs(d - np.median(d)))),
                             hist=hist.tolist(), frac_within_30=float((np.abs(d) < 30).mean()), frac_within_30_of_90=float((np.abs(d2) < 30).mean()))
# also for higher amplitude
for thr in [0.05, 0.1]:
    Hh = H[H.IPDgofha > thr]
    if len(Hh):
        dd = Hh.dphi
        print(f'   IPDgofha>{thr}: n={len(Hh)} median {np.median(dd):+.1f} MAD-sd {1.4826*np.median(np.abs(dd-np.median(dd))):.1f}; frac |d|<30: {(np.abs(dd)<30).mean():.2f}')
# phase distribution for low amplitude like target (0.02-0.04)
Hl = H[(H.IPDgofha < 0.04)]
dd = Hl.dphi
print(f'   0.02<IPDgofha<0.04: n={len(Hl)} median {np.median(dd):+.1f} MAD-sd {1.4826*np.median(np.abs(dd-np.median(dd))):.1f}; frac |d|<30: {(np.abs(dd)<30).mean():.2f}')
out['phase_lowamp'] = dict(n=int(len(Hl)), median=float(np.median(dd)), mad_sd=float(1.4826 * np.median(np.abs(dd - np.median(dd)))),
                           frac_within_30=float((np.abs(dd) < 30).mean()))
# singles: harmonic amplitude distribution
print(f'   singles IPDgofha: median {np.median(C.IPDgofha):.3f}, 84% {np.percentile(C.IPDgofha,84):.3f}, 95% {np.percentile(C.IPDgofha,95):.3f}; target 0.029')
out['singles_gofha_q'] = [float(np.percentile(C.IPDgofha, q)) for q in (50, 84, 95, 99)]
# (e) target-like
T = B[(B.Gmag > 10.5) & (B.Gmag < 13.0) & (B.RUWE < 1.4) & (B.IPDfmp >= 25) & (B.IPDfmp <= 60)]
print(f'\n(e) target-like (G 10.5-13, RUWE<1.4, IPDfmp 25-60, Gaia-unresolved): n={len(T)}')
print(T[['sep', 'dm', 'filt', 'Gmag', 'BP-RP', 'RUWE', 'IPDfmp', 'IPDgofha', 'IPDgofhp', 'pa', 'bfrac', 'epsi', 'cat']].sort_values('sep').to_string())
T2 = T[T.bfrac >= 0.8]
print(f'   ... with BP blended fraction >= 0.8: n={len(T2)} seps {np.round(T2.sep.values,2).tolist()} dm {np.round(T2.dm.values,2).tolist()}')
out['target_like'] = dict(n=int(len(T)), sep=T.sep.round(3).tolist(), dm=T.dm.round(2).tolist(), bfrac=T.bfrac.round(2).tolist(),
                          n_bfrac08=int(len(T2)), sep_bfrac08=T2.sep.round(3).tolist(), dm_bfrac08=T2.dm.round(2).tolist())
# broader: all binaries with bfrac>=0.9 and fmp 20-60, any G 9-14, RUWE<1.4
T3 = B[(B.RUWE < 1.4) & (B.IPDfmp >= 20) & (B.IPDfmp <= 60) & (B.bfrac >= 0.9)]
print(f'   G 9-14, RUWE<1.4, IPDfmp 20-60, bfrac>=0.9: n={len(T3)}; sep pct 5/16/50/84/95 '
      f'{np.round(np.percentile(T3.sep,[5,16,50,84,95]),2).tolist() if len(T3) else None}; dm pct {np.round(np.percentile(T3.dm,[5,16,50,84,95]),2).tolist() if len(T3) else None}')
out['fmp20_60_bfrac09'] = dict(n=int(len(T3)), sep=T3.sep.round(3).tolist(), dm=T3.dm.round(2).tolist(), G=T3.Gmag.round(2).tolist())
# (f) field base rates
F = pd.read_csv(f'{R}/gaia/field3deg_G10.75-12.75.csv')
A = Table.read(f'{R}/gaia/cone3deg_all.fits').to_pandas()
def xyz(ra, de):
    ra = np.radians(ra); de = np.radians(de)
    return np.c_[np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)]
tree = cKDTree(xyz(A.ra.values, A.dec.values))
dd, ii = tree.query(xyz(F.ra.values, F.dec.values), k=2)
F['nn'] = np.degrees(dd[:, 1]) * 3600
F['nn_g'] = A.phot_g_mean_mag.values[ii[:, 1]]
F['bfrac'] = F.phot_bp_n_blended_transits / F.phot_bp_n_obs
iso = F[F.nn > 5.0]
print(f'\n(f) field (3 deg, G 10.75-12.75): n={len(F)}, isolated (no Gaia source within 5"): {len(iso)}')
for name, s in [('all', F), ('isolated>5"', iso), ('nn<5"', F[F.nn <= 5])]:
    print(f'   {name:12s}: n={len(s):5d} frac IPDfmp>=30 {(s.ipd_frac_multi_peak>=30).mean():.4f} (n={(s.ipd_frac_multi_peak>=30).sum()}), '
          f'>=41 {(s.ipd_frac_multi_peak>=41).mean():.4f}; frac bfrac>=0.9 {(s.bfrac>=0.9).mean():.3f}; both {((s.ipd_frac_multi_peak>=30)&(s.bfrac>=0.9)).mean():.4f} (n={((s.ipd_frac_multi_peak>=30)&(s.bfrac>=0.9)).sum()})')
out['field'] = dict(n=int(len(F)), n_iso=int(len(iso)), iso_frac_fmp_ge30=float((iso.ipd_frac_multi_peak >= 30).mean()),
                    iso_n_fmp_ge30=int((iso.ipd_frac_multi_peak >= 30).sum()), iso_frac_bfrac_ge09=float((iso.bfrac >= 0.9).mean()),
                    iso_n_both=int(((iso.ipd_frac_multi_peak >= 30) & (iso.bfrac >= 0.9)).sum()),
                    iso_frac_fmp_ge41=float((iso.ipd_frac_multi_peak >= 41).mean()))
# among isolated field stars with fmp>=30: blend fraction, ruwe, harmonic amplitude
s = iso[iso.ipd_frac_multi_peak >= 30]
print('   isolated with IPDfmp>=30: bfrac quartiles', np.round(np.percentile(s.bfrac, [25, 50, 75]), 2).tolist() if len(s) else None,
      ' RUWE median', round(float(np.median(s.ruwe)), 2) if len(s) else None)
# blend fraction distribution in isolated field stars by fmp
for lo, hi in [(0, 1), (1, 10), (10, 30), (30, 101)]:
    q = iso[(iso.ipd_frac_multi_peak >= lo) & (iso.ipd_frac_multi_peak < hi)]
    print(f'   isolated IPDfmp {lo}-{hi}: n={len(q)} frac bfrac>=0.9 {(q.bfrac>=0.9).mean():.3f}  median bfrac {np.nanmedian(q.bfrac):.2f}')
json.dump(out, open(f'{R}/results/v06_gaia_calib.json', 'w'), indent=1)
S.to_csv(f'{O}/calib_sample_v06.csv', index=False)
