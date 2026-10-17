"""r06b: referee-2 TRICERATOPS re-run (venv_tri) with my own inputs, to check the analyst's main configuration.
Differences from fpp/f05-f06:
  * light curve: REF photometry (referee_phot/lc_M3V_referee_windows.csv), my own per-window linear baseline fitted to
    |dt| > 0.06 d and divided out, 200-s sectors 66/93/103/104, folded on T0 = 3840.9141, binned at 200 s, |dt| < 0.25 d;
    flux_err = rms of the out-of-transit bins (|dt| > 0.08 d)  [variant 'prop': median propagated bin error x beta]
  * contrast curve (NOT high-resolution imaging): Gaia-IPD 'red_variant' cells (calibrated on red stars) with LR_hi84 < 0.2
    contiguous from dm = 0, and the referee-corrected DECaPS2 5-sigma limits; conservative running minimum from large sep
  * apertures: pixels whose centres are within 2.5 px of the M3V
  * neighbours: only the target and the two catalogued stars inside the localisation region (TIC 1515709327, 1515709330)
Usage: r06b_tri_run.py <seed> <N> <err: emp|prop> [ebcirc] [cccons = use the analyst's conservative curve fpp/data/tri/cc_cons.txt]   -> out/tri/r2_<err>_s<seed>.json"""
import sys, os, json, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from astropy.io import fits
R2 = 'tic214311908/analysis/referee2_fpp'
B = 'tic214311908/analysis'
CUBES = '/home/user/work/tesscut5/214311906'
seed, N, errmode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
ebcirc = 'ebcirc' in sys.argv[4:]
ccsel = 'cons' if 'cccons' in sys.argv[4:] else 'r2'
tag = f'r2_{errmode}{"_circ" if ebcirc else ""}{"_cccons" if ccsel == "cons" else ""}_s{seed}'
os.makedirs(f'{R2}/out/tri', exist_ok=True); os.makedirs(f'{R2}/data/tri', exist_ok=True)
P, T0 = 5.6915087, 3840.9141
# ---- light curve
d = pd.read_csv(f'{B}/referee_phot/lc_M3V_referee_windows.csv')
d = d[np.isfinite(d.flux_norm) & d.sector.isin([66, 93, 103, 104])]
ep = np.round((d.time_btjd_tdb.values - T0)/P).astype(int); dt = d.time_btjd_tdb.values - (T0 + ep*P)
f = d.flux_norm.values.copy(); e = d.flux_err.values
beta = pd.read_csv(f'{B}/referee_phot/out/per_sector.csv').set_index('sector')
bmult = np.array([beta.loc[s, 'err_adopt']/beta.loc[s, 'err_white'] for s in d.sector.values])
for w_ in np.unique(d.window_id.values):
    m = d.window_id.values == w_; o = m & (np.abs(dt) > 0.06)
    c = np.polyfit(dt[o], f[o], 1, w=1/e[o]); f[m] = f[m]/np.polyval(c, dt[m])
m = np.abs(dt) < 0.25
dt, f, e, bmult = dt[m], f[m], e[m], bmult[m]
wd = 200/86400.; edges = np.arange(-0.25, 0.25 + wd/2, wd); idx = np.digitize(dt, edges)
tb, fb, eb = [], [], []
for k in np.unique(idx):
    q = idx == k
    if q.sum() < 3: continue
    ww = 1/e[q]**2
    tb.append(np.sum(dt[q]*ww)/ww.sum()); fb.append(np.sum(f[q]*ww)/ww.sum()); eb.append(np.sqrt(np.sum((ww*e[q]*bmult[q])**2))/ww.sum())
tb, fb, eb = map(np.array, (tb, fb, eb))
oot = np.abs(tb) > 0.08
fb = fb/np.average(fb[oot], weights=1/eb[oot]**2)
sig = float(np.std(fb[oot])) if errmode == 'emp' else float(np.median(eb))
np.savez(f'{R2}/data/tri/lc_ref200_r2.npz', t=tb, f=fb, e=eb)
depth = float(1 - np.average(fb[np.abs(tb) < 0.02], weights=1/eb[np.abs(tb) < 0.02]**2))
# ---- contrast curve (proxy)
ipd = json.load(open(f'{B}/star/out/ipd_limits.json'))['red_variant']['cells']
seps = sorted(set(tuple(c['sep']) for c in ipd))
cell = {}
for s in seps:
    dm = 0.0
    for c in sorted([c for c in ipd if tuple(c['sep']) == s], key=lambda c: c['dm'][0]):
        if c['n'] >= 5 and c['LR_hi84'] == c['LR_hi84'] and c['LR_hi84'] < 0.2: dm = c['dm'][1]
        else: break
    cell[s] = dm
dec = {1.0: 4.24, 1.25: 5.07, 1.5: 5.80, 2.0: 5.58, 2.5: 6.17, 3.0: 6.89, 4.0: 7.53, 5.0: 7.88, 6.0: 8.59}
grid = np.unique(np.concatenate([[0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9], list(dec.keys()), [a for a, b in seps], [b for a, b in seps]]))
rows = []
for s in grid:
    v_ipd = max([dm for (a, b), dm in cell.items() if a <= s <= b] or [0.0])
    ks = np.array(sorted(dec)); v_dec = float(np.interp(s, ks, [dec[k] for k in ks])) if s >= 1.0 else 0.0
    rows.append((float(s), max(v_ipd, v_dec)))
cc = []; mn = 99
for s, v in rows[::-1]:
    mn = min(mn, v); cc.append((s, max(mn, 0.01)))
cc = cc[::-1]
np.savetxt(f'{R2}/data/tri/cc_r2.txt', np.array(cc), fmt='%.3f', delimiter=',')
print('depth', round(depth, 5), 'sig', round(sig, 5), 'nbins', len(tb), 'cc', cc, flush=True)
# ---- TRICERATOPS
os.chdir(f'{R2}/data/tri')
np.random.seed(seed)
import triceratops.triceratops as tr
class _Cut:
    def __init__(self, path): self.hdu = fits.open(path, memmap=True)
class _Search:
    def __init__(self, sector): self.sector = sector
    def download_all(self, cutout_size=None):
        ff = [x for x in os.listdir(CUBES) if x.startswith(f'tess-s{self.sector:04d}-')]
        return [_Cut(f'{CUBES}/{ff[0]}')]
tr.lightkurve.search_tesscut = lambda target=None, sector=None, **kw: _Search(int(sector))
if ebcirc:
    import triceratops.marginal_likelihoods as _ml
    _o = _ml.sample_ecc
    _ml.sample_ecc = lambda x, planet, P_orb: _o(x, planet, P_orb) if planet else np.zeros_like(x)
t0 = time.time()
SECT = [66, 93, 103, 104]
target = tr.target(ID=214311908, sectors=np.array(SECT), search_radius=10)
target.stars.loc[0, 'mass'] = 0.370; target.stars.loc[0, 'rad'] = 0.381; target.stars.loc[0, 'Teff'] = 3350.; target.stars.loc[0, 'plx'] = 15.549
aps = []
for k in range(len(target.pix_coords)):
    x0, y0 = target.pix_coords[k][0]
    xs = np.arange(np.floor(x0) - 4, np.ceil(x0) + 5); ys = np.arange(np.floor(y0) - 4, np.ceil(y0) + 5)
    XX, YY = np.meshgrid(xs, ys); mm = np.hypot(XX - x0, YY - y0) <= 2.5
    aps.append(np.array([XX[mm], YY[mm]]).T)
target.calc_depths(tdepth=depth, all_ap_pixels=aps)
keep = ['214311908', '1515709327', '1515709330']
target.remove_star([i for i in target.stars.ID.astype(str).values if i not in keep])
target.calc_probs(time=tb, flux_0=fb, flux_err_0=sig, P_orb=P, parallel=True, N=N, verbose=0, exptime=200/86400., nsamples=7,
                  contrast_curve_file=(f'{R2}/data/tri/cc_r2.txt' if ccsel == 'r2' else f'{B}/fpp/data/tri/cc_cons.txt'), filt='TESS')
df = target.probs.copy()
res = dict(FPP=float(target.FPP), NFPP=float(target.NFPP), seed=seed, N=N, errmode=errmode, flux_err_0=sig, tdepth=depth,
           eb_circular=ebcirc, lc_tag='ref200_r2', cc=ccsel, nearby='loc', scenarios=df.astype({'ID': str}).to_dict('records'),
           runtime_s=time.time() - t0)
json.dump(res, open(f'{R2}/out/tri/{tag}.json', 'w'), indent=1, default=str)
print(df.sort_values('prob', ascending=False).head(8).to_string())
print('FPP', target.FPP, 'NFPP', target.NFPP, 'runtime', time.time() - t0)
