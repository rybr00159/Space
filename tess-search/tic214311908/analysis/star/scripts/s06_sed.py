"""SED fit of TIC 214311908 with BT-Settl CIFIST2011 models (logg 5.0, [M/H] = 0; logg 4.5/5.5 as systematics).
Photometry: Gaia DR3 G/BP/RP, 2MASS JHKs, SkyMapper DR4 griz (AB), unWISE W1/W2 (AllWISE W1/W2 are halo-flagged
by the Mira 17" away and W3/W4 are dominated by it; they are excluded). Synthetic photometry: photon-counting
<F_lambda> over SVO filter curves; observed magnitudes converted with SVO zero points (Vega, or 3631 Jy for AB).
Free parameters: Teff (linear interpolation of log flux between 100 K grid points) and the dilution factor
(R/d)^2. An error floor of 0.02 mag (optical/NIR) and 0.05 mag (WISE) is added in quadrature.
Monte Carlo over photometric errors (2000 draws) and over the distance gives Teff_SED, Fbol, L, R_SED.
Output: out/sed.json, figs/sed.png"""
import json, glob, os, warnings
warnings.filterwarnings('ignore')
import numpy as np
from astroquery.svo_fps import SvoFps
B = 'tic214311908/analysis/star'
M = f'{B}/data/models'
c_A = 2.99792458e18  # A/s
phot = {  # band: (mag, err, system, filterID)
    'G': (14.281043, 0.0006, 'Vega', 'GAIA/GAIA3.G'), 'BP': (15.830741, 0.0035, 'Vega', 'GAIA/GAIA3.Gbp'),
    'RP': (13.074592, 0.0013, 'Vega', 'GAIA/GAIA3.Grp'),
    'J': (11.363, 0.028, 'Vega', '2MASS/2MASS.J'), 'H': (10.740, 0.031, 'Vega', '2MASS/2MASS.H'), 'Ks': (10.482, 0.026, 'Vega', '2MASS/2MASS.Ks'),
    'SM_g': (15.8171, 0.0299, 'AB', 'SkyMapper/SkyMapper.g'),
    'SM_i': (13.3869, 0.0100, 'AB', 'SkyMapper/SkyMapper.i'), 'SM_z': (12.7683, 0.0098, 'AB', 'SkyMapper/SkyMapper.z'),
    'W1': (22.5 - 2.5 * np.log10(73687.0), 0.0007, 'Vega', 'WISE/WISE.W1'), 'W2': (22.5 - 2.5 * np.log10(80969.0), 0.0015, 'Vega', 'WISE/WISE.W2'),
}
check = {'SM_r': (14.9081, 0.0239, 'AB', 'SkyMapper/SkyMapper.r'), 'W1_allwise': (10.112, 0.022, 'Vega', 'WISE/WISE.W1'), 'W2_allwise': (9.884, 0.020, 'Vega', 'WISE/WISE.W2'),
         'DENIS_I': (13.094, 0.04, 'Vega', 'DENIS/DENIS.I'), 'SM_v': (18.192, 0.059, 'AB', 'SkyMapper/SkyMapper.v'),
         'DECaPS_g': (16.177, 0.02, 'AB', 'CTIO/DECam.g'), 'DECaPS_r': (14.737, 0.02, 'AB', 'CTIO/DECam.r'),
         'DECaPS_i': (13.355, 0.02, 'AB', 'CTIO/DECam.i'), 'DECaPS_z': (12.761, 0.02, 'AB', 'CTIO/DECam.z'),
         'DECaPS_Y': (12.563, 0.02, 'AB', 'CTIO/DECam.Y'), 'TESS': (np.nan, 0, 'Vega', 'TESS/TESS.Red')}
floor = {k: (0.05 if k.startswith('W') else 0.02) for k in phot}
fcache = f'{M}/filters_all.npz'
allf = {v[3] for v in list(phot.values()) + list(check.values())}
if os.path.exists(fcache):
    FT = dict(np.load(fcache, allow_pickle=True)); FT = {k: FT[k] for k in FT}
else:
    FT = {}
info_fn = f'{M}/filters_zp.json'
INFO = json.load(open(info_fn)) if os.path.exists(info_fn) else {}
for f in allf:
    key = f.replace('/', '__')
    if key not in FT:
        t = SvoFps.get_transmission_data(f)
        FT[key] = np.vstack([np.array(t['Wavelength'], float), np.array(t['Transmission'], float)])
    if f not in INFO:
        import requests, re
        txt = requests.get(f'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={f}', timeout=60).text
        kv = dict(re.findall(r'name="(ZeroPoint|WavelengthPivot|ZeroPointType|MagSys)"[^>]*value="([^"]+)"', txt))
        INFO[f] = dict(zp_vega_jy=float(kv['ZeroPoint']), lpivot=float(kv['WavelengthPivot']), zptype=kv['ZeroPointType'], magsys=kv['MagSys'])
np.savez(fcache, **FT); json.dump(INFO, open(info_fn, 'w'), indent=1)

def obs_flam(mag, sysm, f):
    zp = 3631.0 if sysm == 'AB' else INFO[f]['zp_vega_jy']
    fnu = zp * 10 ** (-0.4 * mag) * 1e-23  # erg/s/cm2/Hz
    lp = INFO[f]['lpivot']
    return fnu * c_A / lp ** 2

def load(teff, logg=5.0):
    a = np.load(f'{M}/btsettl_{teff}_{logg:.1f}.npy').astype(float)
    return a[0], a[1]

def synth(lam, fl, f):
    w, T = FT[f.replace('/', '__')]
    Ti = np.interp(lam, w, T, left=0, right=0)
    m = Ti > 0
    return np.trapezoid(fl[m] * Ti[m] * lam[m], lam[m]) / np.trapezoid(Ti[m] * lam[m], lam[m])

def model_fbol(teff, logg=5.0):
    lam, fl = load(teff, logg)
    fn = f'{M}/btsettl_{teff}_{logg:.1f}_fbol.npy'
    tot = float(np.load(fn)[0]) if os.path.exists(fn) else np.trapezoid(fl, lam) + fl[-1] * lam[-1] / 3
    return tot

grid = {}
for logg in (5.0, 4.5, 5.5):
    for fn in sorted(glob.glob(f'{M}/btsettl_*_{logg:.1f}.npy')):
        teff = int(os.path.basename(fn).split('_')[1])
        lam, fl = load(teff, logg)
        grid[(teff, logg)] = dict(s={k: synth(lam, fl, v[3]) for k, v in {**phot, **check}.items()}, fbol=model_fbol(teff, logg))
bands = list(phot)
yobs = np.array([obs_flam(phot[b][0], phot[b][2], phot[b][3]) for b in bands])
emag = np.array([np.hypot(phot[b][1], floor[b]) for b in bands])

def interp_model(teff, logg=5.0):
    ts = sorted(t for (t, g) in grid if g == logg)
    t0 = max([t for t in ts if t <= teff], default=ts[0]); t1 = min([t for t in ts if t > teff], default=ts[-1])
    if t1 == t0: w = 0.0
    else: w = (teff - t0) / (t1 - t0)
    s0 = grid[(t0, logg)]['s']; s1 = grid[(t1, logg)]['s']
    out = {k: np.exp((1 - w) * np.log(s0[k]) + w * np.log(s1[k])) for k in s0}
    fb = np.exp((1 - w) * np.log(grid[(t0, logg)]['fbol']) + w * np.log(grid[(t1, logg)]['fbol']))
    return out, fb

def fit(y, e_mag, logg=5.0, tgrid=np.arange(2800, 4000, 5.0), use=None):
    use = np.ones(len(bands), bool) if use is None else use
    best = None
    for T in tgrid:
        s, fb = interp_model(T, logg)
        m = np.array([s[b] for b in bands])
        # fit in magnitudes: dm = -2.5 log10(y/m) - logscale
        d = -2.5 * np.log10(y / m)
        w = 1 / e_mag ** 2
        mu = np.sum((w * d)[use]) / np.sum(w[use])
        chi2 = np.sum((w * (d - mu) ** 2)[use])
        if best is None or chi2 < best[0]:
            best = (chi2, T, 10 ** (-0.4 * mu), fb, d - mu)
    return best

chi2, T, scale, fb, resid = fit(yobs, emag)
print(f'best fit logg 5.0: Teff {T:.0f} K, chi2 {chi2:.1f} for {len(bands)-2} dof; scale (R/d)^2 = {scale:.4e}')
for b, r, e in zip(bands, resid, emag): print(f'   {b:6s} resid {r*-1:+.3f} mag (obs - model; err {e:.3f})')
# Monte Carlo
plx = 15.5494; eplx = 0.036
rng = np.random.default_rng(1)
res = []
for k in range(1500):
    ypert = yobs * 10 ** (-0.4 * rng.normal(0, emag))
    c2, Tk, sk, fbk, _ = fit(ypert, emag, tgrid=np.arange(3000, 3800, 5.0))
    dk = 1000 / rng.normal(plx, eplx)  # pc
    Fbol = sk * fbk
    Lk = 4 * np.pi * (dk * 3.0857e18) ** 2 * Fbol / 3.828e33
    Rk = np.sqrt(sk) * dk * 3.0857e18 / 6.957e10
    res.append((Tk, Fbol, Lk, Rk, c2))
res = np.array(res)
q = lambda x: [float(np.percentile(x, 50)), float(np.percentile(x, 50) - np.percentile(x, 16)), float(np.percentile(x, 84) - np.percentile(x, 50))]
out = dict(best=dict(teff=T, chi2=chi2, dof=len(bands) - 2, scale=scale, resid_obs_minus_model=dict(zip(bands, (-resid).round(4).tolist()))),
           mc=dict(teff=q(res[:, 0]), fbol_cgs=q(res[:, 1]), L_Lsun=q(res[:, 2]), R_Rsun=q(res[:, 3])))
# logg systematics and leave-band-group-out
for logg in (4.5, 5.5):
    try:
        c2, Tg, sg, fbg, _ = fit(yobs, emag, logg=logg, tgrid=np.arange(3000, 3700, 5.0))
        out[f'logg{logg}'] = dict(teff=Tg, chi2=c2, fbol=sg * fbg, R=float(np.sqrt(sg) * 1000 / plx * 3.0857e18 / 6.957e10))
        print(f'logg {logg}: Teff {Tg:.0f} chi2 {c2:.1f} Fbol {sg*fbg:.4e}')
    except Exception as e:
        print('logg', logg, e)
for lab, drop in [('no_SkyMapper', ['SM_g', 'SM_i', 'SM_z']), ('no_WISE', ['W1', 'W2']), ('no_BP', ['BP']),
                  ('only_Gaia_2MASS', ['SM_g', 'SM_i', 'SM_z', 'W1', 'W2'])]:
    use = np.array([b not in drop for b in bands])
    c2, Tg, sg, fbg, _ = fit(yobs, emag, use=use, tgrid=np.arange(3000, 3800, 5.0))
    out[lab] = dict(teff=Tg, chi2=c2, dof=int(use.sum() - 2), fbol=sg * fbg)
    print(f'{lab}: Teff {Tg:.0f} chi2 {c2:.1f}/{use.sum()-2} Fbol {sg*fbg:.4e}')
# predicted magnitudes in check bands at best fit (including TESS)
s, fbb = interp_model(T)
pred = {}
for k, v in check.items():
    mflam = s[k] * scale
    zp = 3631.0 if v[2] == 'AB' else INFO[v[3]]['zp_vega_jy']
    lp = INFO[v[3]]['lpivot']
    mag = -2.5 * np.log10(mflam * lp ** 2 / c_A / 1e-23 / zp)
    pred[k] = dict(model=float(mag), obs=(None if not np.isfinite(v[0]) else v[0]))
    print(f'  check {k:11s} model {mag:.3f} obs {v[0]}')
out['check_bands'] = pred
print(json.dumps(out['mc'], indent=0))
json.dump(out, open(f'{B}/out/sed.json', 'w'), indent=1, default=float)
# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lam, fl = load(int(round(T / 100) * 100))
fig, ax = plt.subplots(figsize=(7, 4))
k = (lam > 3000) & (lam < 60000)
lamb = lam[k][::20]; flb = np.convolve(fl[k], np.ones(20) / 20, 'same')[::20]
ax.plot(lamb / 1e4, lamb * flb * scale, color='0.6', lw=0.6, label=f'BT-Settl {int(round(T/100)*100)} K (scaled)')
lp = np.array([INFO[phot[b][3]]['lpivot'] for b in bands])
ax.errorbar(lp / 1e4, lp * yobs, yerr=lp * yobs * emag / 1.086, fmt='o', color='C0', label='observed')
s_b = np.array([s[b] for b in bands])
ax.plot(lp / 1e4, lp * s_b * scale, 's', mfc='none', color='C3', label='model synthetic')
for kk in ['W1_allwise', 'W2_allwise']:
    v = check[kk]; l0 = INFO[v[3]]['lpivot']; ax.plot(l0 / 1e4, l0 * obs_flam(v[0], v[2], v[3]), 'x', color='C1')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('wavelength (micron)'); ax.set_ylabel(r'$\lambda F_\lambda$ (erg s$^{-1}$ cm$^{-2}$)')
ax.set_title(f'TIC 214311908 SED: Teff {T:.0f} K, chi2 {chi2:.1f}/{len(bands)-2}; x = AllWISE (halo-flagged)')
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f'{B}/figs/sed.png', dpi=120)
