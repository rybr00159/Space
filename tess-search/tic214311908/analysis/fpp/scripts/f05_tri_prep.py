"""f05: inputs for TRICERATOPS.
1. Folded, binned light curves (time = days from mid-transit), from the deblended, detrended M3V photometry:
   * inv200 : INV (phot model B, RRab leak removed: search/data/m3v_clean.npz det_clean), 200-s sectors 66/93/103/104,
              binned in phase at 200 s (TRICERATOPS exptime = 200 s).
   * inv600 : INV, sectors 39/66/93/103/104, binned at 600 s (exptime = 600 s; S39 is 600-s data; the 200-s data
              binned to 600 s are smeared like a 600-s exposure).
   * ref200 : REF (referee local-window photometry, shape/data/ref_win.npz), 200-s sectors, binned at 200 s.
   Per-point errors include the per-sector white scale and red-noise beta (shape/out/noise.json); bin errors are
   1/sqrt(sum w). |dt| < 0.25 d. S12 (30-min) is not used: TRICERATOPS takes a single exposure time.
2. Contrast curves (data/tri/cc_cons.txt baseline, cc_opt.txt variant): NOT a high-resolution imaging contrast curve. It is the envelope of
   (a) the Gaia IPD/RUWE calibration (star/out/ipd_limits.json; contrast at which companions are 'excluded', LR<~0.1, in
   I band, used as a proxy for TESS) and (b) the DECaPS2 residual-image 5-sigma limits as CORRECTED by the referee
   (referee_star/out/r04_decaps_inject.json: i 4.24/5.07/5.80/5.58 and z 3.96/4.53/4.96/5.24 at 1/1.25/1.5/2"; the
   investigator's curve beyond 2"). TRICERATOPS rules out any simulated companion brighter than the curve.
3. Writes data/tri/inputs.json with depth, ephemeris and stellar overrides."""
import numpy as np, json, sys
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
os.makedirs(f'{F}/data/tri', exist_ok=True)

def fold_bin(d, width_s, halfwin=0.25):
    dt = d['dt']; m = np.abs(dt) < halfwin
    dt, f, s = dt[m], d['f'][m], d['sig'][m]
    w = width_s/86400.
    edges = np.arange(-halfwin, halfwin + w/2, w)
    idx = np.digitize(dt, edges)
    tb, fb, eb, nb = [], [], [], []
    for k in np.unique(idx):
        q = idx == k
        if q.sum() < 3: continue
        ww = 1/s[q]**2
        tb.append(np.sum(dt[q]*ww)/ww.sum()); fb.append(np.sum(f[q]*ww)/ww.sum()); eb.append(1/np.sqrt(ww.sum())); nb.append(q.sum())
    tb, fb, eb, nb = map(np.array, (tb, fb, eb, nb))
    # normalise the out-of-transit level to 1 (|dt| > 0.1 d)
    oot = np.abs(tb) > 0.10
    c = np.sum(fb[oot]/eb[oot]**2)/np.sum(1/eb[oot]**2)
    fb = fb/c
    return dict(t=tb, f=fb, e=eb, n=nb, oot_level=c)

out = {}
d = inv_clean_windows(halfwin=0.30, sectors=[66, 93, 103, 104])
out['inv200'] = fold_bin(d, 200)
d = inv_clean_windows(halfwin=0.30, sectors=[39, 66, 93, 103, 104])
out['inv600'] = fold_bin(d, 600)
d = ref_windows(sectors=[66, 93, 103, 104])
out['ref200'] = fold_bin(d, 200)
info = {}
for k, v in out.items():
    np.savez(f'{F}/data/tri/lc_{k}.npz', **v)
    intr = np.abs(v['t']) < 0.025
    oot = np.abs(v['t']) > 0.10
    info[k] = dict(nbins=int(len(v['t'])), median_err=float(np.median(v['e'])), rms_oot=float(np.std(v['f'][oot])),
                   chi2_oot_per_dof=float(np.mean(((v['f'][oot]-1)/v['e'][oot])**2)),
                   depth_core=float(1 - np.average(v['f'][intr], weights=1/v['e'][intr]**2)), oot_level=float(v['oot_level']))
    print(k, info[k])

# contrast curve: envelope of IPD 'excluded' contrasts and the corrected DECaPS2 limits (Delta mag, TESS proxy)
ipd = jload('star/out/ipd_limits.json')['cells']
seps_ipd = sorted(set(tuple(c['sep']) for c in ipd))
cell_dm = {}
for s in seps_ipd:
    cells = sorted([c for c in ipd if tuple(c['sep']) == s], key=lambda c: c['dm'][0])
    dm_ex = 0.0
    for c in cells:   # contiguous run of cells with LR_hi84 < 0.2 (strongly disfavoured) from dm = 0 upward
        if c['n'] >= 5 and c['LR_hi84'] == c['LR_hi84'] and c['LR_hi84'] < 0.2: dm_ex = c['dm'][1]
        else: break
    cell_dm[s] = dm_ex
cc_ipd = [(s[0], v) for s, v in cell_dm.items()] + [(s[1], v) for s, v in cell_dm.items()]
r04 = jload('referee_star/out/r04_decaps_inject.json')
dec = {1.0: max(4.24, 3.96), 1.25: max(5.07, 4.53), 1.5: max(5.80, 4.96), 2.0: max(5.58, 5.24)}
inv = jload('star/out/decaps_contrast.json')
for row in inv['i']['curve']:
    if row['sep'] >= 2.5: dec[row['sep']] = row['delta_mag_5sig']
seps = np.array(sorted(set([s for s, _ in cc_ipd] + list(dec.keys()))))
def ipd_at(s):
    v = [dm for (s0, s1), dm in cell_dm.items() if s0 <= s <= s1]
    return max(v) if v else 0.0
def dec_at(s):
    ks = np.array(sorted(dec.keys())); vs = np.array([dec[k] for k in ks])
    if s < ks[0]: return 0.0
    return float(np.interp(s, ks, vs))
rows = []
for s in np.unique(np.concatenate([seps, [0.08, 0.1, 0.5, 0.75, 3.0, 4.0, 5.0, 6.0], [0.5*(a+b) for a, b in seps_ipd]])):
    if s < 0.08: continue
    rows.append((float(s), max(ipd_at(s), dec_at(s))))
# TRICERATOPS interpolates separation at contrast, so the curve must be monotonic. Two versions:
#  cc_opt : running maximum from small separations (optimistic: a contrast reached at small sep is assumed at larger sep)
#  cc_cons: running minimum from large separations (conservative: baseline)
cc_opt = []; mx = 0
for s_, dm in rows:
    mx = max(mx, dm); cc_opt.append((s_, mx))
cc_cons = []; mn = 99
for s_, dm in rows[::-1]:
    mn = min(mn, dm); cc_cons.append((s_, mn))
cc_cons = cc_cons[::-1]
cc_opt = [(a, b) for a, b in cc_opt if b > 0]
cc_cons = [(a, b if b > 0 else 0.01) for a, b in cc_cons]
np.savetxt(f'{F}/data/tri/cc_opt.txt', np.array(cc_opt), fmt='%.3f', delimiter=',')
np.savetxt(f'{F}/data/tri/cc_cons.txt', np.array(cc_cons), fmt='%.3f', delimiter=',')
cc = dict(opt=cc_opt, cons=cc_cons)
print('contrast curve opt (sep arcsec, dmag):', cc_opt)
print('contrast curve cons (sep arcsec, dmag):', cc_cons)
inputs = dict(T0=T0, P=P, tdepth=0.0157, M=M3V['M'], R=M3V['R'], Teff=M3V['Teff'], plx=15.549,
              lc=info, contrast_curve=cc,
              notes=['tdepth = on-target depth of the deblended light curve (15.3-15.4 e-/s of 979 e-/s; LD fit 15.8-15.9 ppt)',
                     'contrast curve is NOT high-resolution imaging; see docstring'])
jdump(inputs, f'{F}/data/tri/inputs.json')
