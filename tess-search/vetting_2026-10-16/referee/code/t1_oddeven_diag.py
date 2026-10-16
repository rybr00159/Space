"""TIC 90086266 S104: is the odd/even depth alternation instrumental?
Per-event diagnostics (background, centroids, pointing, quality flags, flux
level), pixel-by-pixel odd/even ratio, and a 2P-period consistency check."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.io import fits
from reflib import *
from diffimg import *

TIC = 90086266
EV = json.load(open(f'{BASE}/results/t1_events.json'))
T14 = 3.0 / 24
OUT = {}

# all-quality cutout to see flags near events
cutA = load_cut(TIC, 104, qual0=False)
t_all, q_all = cutA['time'], cutA['quality']
qlp = fits.open(lc_path(TIC, 104, 'QLP'))[1].data
tq, qq = np.array(qlp['TIME']), np.array(qlp['QUALITY'])
t, f, e, ex = aperture_lc(TIC, 104, r=2.0)
evs = EV['S104_CUT:r2.0']['events']
rows = []
for v in evs:
    tc = v['tc']
    w_in = np.abs(t - tc) < 0.4 * T14
    w_out = (np.abs(t - tc) > 0.6 * T14) & (np.abs(t - tc) < 0.6 * T14 + 0.2)
    win = np.abs(t_all - tc) < 0.3
    flags = np.unique(q_all[win][q_all[win] > 0])
    nflag = int(np.sum(q_all[win] > 0))
    wq = np.abs(tq - (tc - 448.6 / 86400)) < 0.3
    qflags = np.unique(qq[wq][qq[wq] > 0])
    row = dict(E=v['E'], tc=tc, depth=v['depth'], depth_err=v['depth_err'],
               bkg_in=float(np.median(ex['SAP_BKG'][w_in])), bkg_out=float(np.median(ex['SAP_BKG'][w_out])),
               cx_in_minus_out=float(np.median(ex['MOM_CENTR1'][w_in]) - np.median(ex['MOM_CENTR1'][w_out])),
               cy_in_minus_out=float(np.median(ex['MOM_CENTR2'][w_in]) - np.median(ex['MOM_CENTR2'][w_out])),
               cx_out=float(np.median(ex['MOM_CENTR1'][w_out])), cy_out=float(np.median(ex['MOM_CENTR2'][w_out])),
               pc1=float(np.median(ex['POS_CORR1'][w_out])), pc2=float(np.median(ex['POS_CORR2'][w_out])),
               flux_out=float(np.median(ex['raw'][w_out])), n_flagged_cut_pm0p3d=nflag, cut_flags=flags.tolist(), qlp_flags=qflags.tolist(),
               n_in=int(w_in.sum()), n_out=int(w_out.sum()))
    rows.append(row)
print(' E  depth    err    bkg_in  bkg_out  dcx(mpx) dcy(mpx)  cx_out  cy_out  flux_out  nflag  cut_flags  qlp_flags  n_in n_out')
for r in rows:
    print(f"{r['E']:2d} {r['depth']:.5f} {r['depth_err']:.5f} {r['bkg_in']:8.1f} {r['bkg_out']:8.1f} {r['cx_in_minus_out']*1e3:+7.2f} {r['cy_in_minus_out']*1e3:+7.2f} "
          f"{r['cx_out']:7.3f} {r['cy_out']:7.3f} {r['flux_out']:9.0f} {r['n_flagged_cut_pm0p3d']:5d} {str(r['cut_flags']):12s} {str(r['qlp_flags']):10s} {r['n_in']:4d} {r['n_out']:4d}")
odd = np.array([r['E'] % 2 == 1 for r in rows])
for k in ['bkg_in', 'bkg_out', 'cx_in_minus_out', 'cy_in_minus_out', 'cx_out', 'cy_out', 'flux_out', 'n_flagged_cut_pm0p3d']:
    a = np.array([r[k] for r in rows], float)
    print(f'  {k:22s} even mean {a[~odd].mean():12.5f}  odd mean {a[odd].mean():12.5f}  (diff/sd {(a[~odd].mean()-a[odd].mean())/ (a.std()+1e-12):+.2f})')
OUT['per_event'] = rows

# depth vs background correlation
dep = np.array([r['depth'] for r in rows]); bk = np.array([r['bkg_in'] for r in rows])
OUT['corr_depth_bkg'] = float(np.corrcoef(dep, bk)[0, 1])
print('corr(depth, bkg_in) =', OUT['corr_depth_bkg'])

# pixel-level: odd and even difference images, ratio per pixel
cut = load_cut(TIC, 104)
tcs = np.array([v['tc'] for v in evs])
E = np.array([v['E'] for v in evs])
evd, direct = event_diff_images(cut, tcs, T14, fin=0.3, gap=0.6, wout=0.2)
idx_e = [i for i in range(len(evd)) if E[np.argmin(np.abs(tcs - evd[i][0]))] % 2 == 0]
idx_o = [i for i in range(len(evd)) if E[np.argmin(np.abs(tcs - evd[i][0]))] % 2 == 1]
De, Ee = combine(evd, idx_e)
Do, Eo = combine(evd, idx_o)
g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
xg, yg = [float(v) for v in cut['wcs'].world_to_pixel_values(g.iloc[0]['ra_ep'], g.iloc[0]['dec_ep'])]
xc, yc = int(round(xg)), int(round(yg))
sl = (slice(yc - 2, yc + 3), slice(xc - 2, xc + 3))
ratio = (Do / De)[sl]
rerr = np.abs(ratio) * np.sqrt((Eo / Do) ** 2 + (Ee / De) ** 2)[sl]
snr = (De / Ee)[sl]
print('odd/even deficit ratio per pixel (5x5, only pixels with even S/N>8 shown):')
for j in range(4, -1, -1):
    print('   ' + ' '.join(('%5.3f+-%5.3f' % (ratio[j, i], rerr[j, i])) if snr[j, i] > 8 else '     ---     ' for i in range(5)))
m = snr > 8
w = 1 / rerr[m] ** 2
rbar = np.sum(ratio[m] * w) / np.sum(w)
chi2 = np.sum(((ratio[m] - rbar) / rerr[m]) ** 2)
OUT['pixel_ratio'] = dict(mean=rbar, chi2=chi2, n=int(m.sum()), ratio=ratio.tolist(), err=rerr.tolist())
print(f'  weighted mean odd/even ratio = {rbar:.4f}; chi2 about mean = {chi2:.1f} for {m.sum()-1} dof')
# total-flux sum of odd and even deficits in a 5x5 box (aperture-like)
so, se_ = Do[sl].sum(), De[sl].sum()
print(f'  5x5 summed deficit odd/even = {so/se_:.4f}')
OUT['box_ratio'] = so / se_
save_json(OUT, f'{BASE}/results/t1_oddeven_diag.json')
