"""h10b: [Fe/H] prior for the M3V given its space motion. Uses data/gspspec_dwarfs.csv (h10a; best-quality GSP-Spec
dwarfs within 200 pc). UVW (heliocentric, astropy Galactic frame) for those stars and for the M3V (same code);
LSR correction Schoenrich+10 (11.1, 12.24, 7.25) km/s. Kinematic analogues: |W_LSR - 47.6| < 10 km/s and
|v_tot,LSR - 58.4| < 15 km/s (and a looser set: 30 < |W_LSR| < 65). Also Bensby+14 thin/thick probabilities.
[Fe/H] = mh_gspspec + fem_gspspec (fem missing -> mh). No GSP-Spec calibration applied (the Recio-Blanco+23 logg
polynomial is ~0.02 dex for dwarfs). Fractions with [Fe/H] above +0.1/+0.2/+0.3/+0.4, KDE density, for all local
dwarfs and for kinematic analogues; a Gaussian-mixture fit used as the prior in h15.
Output: out/h10b_feh_prior.json, figs/h10b_feh_prior.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
d = pd.read_csv(f'{B}/data/gspspec_dwarfs.csv', dtype={'flags_gspspec': str, 'flag13': str})
d = d[d.flags_gspspec.astype(str).str[:13] == '0000000000000'].copy()
d['feh'] = d.mh_gspspec + d.fem_gspspec.fillna(0.0)
def uvw(ra, dec, plx, pmra, pmdec, rv):
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=(1000 / plx) * u.pc, pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                 radial_velocity=rv * u.km / u.s, frame='icrs').galactic
    v = c.velocity
    return v.d_x.to(u.km / u.s).value, v.d_y.to(u.km / u.s).value, v.d_z.to(u.km / u.s).value
U, V, W = uvw(d.ra.values, d.dec.values, d.parallax.values, d.pmra.values, d.pmdec.values, d.radial_velocity.values)
lsr = np.array([11.1, 12.24, 7.25])
d['U'], d['V'], d['W'] = U + lsr[0], V + lsr[1], W + lsr[2]
d['vtot'] = np.sqrt(d.U ** 2 + d.V ** 2 + d.W ** 2)
tu, tv, tw = uvw(np.array([260.969674]), np.array([-51.323821]), np.array([15.549]), np.array([-113.671152]), np.array([85.54355]), np.array([16.2837]))
T = np.array([tu[0], tv[0], tw[0]]) + lsr; Tv = float(np.sqrt(np.sum(T ** 2)))
print('M3V UVW_LSR', T.round(2), 'vtot', round(Tv, 1))
# Bensby+14 TD/D
def gauss3(U, V, W, sU, sV, sW, Vasym):
    k = 1 / ((2 * np.pi) ** 1.5 * sU * sV * sW)
    return k * np.exp(-U ** 2 / (2 * sU ** 2) - (V - Vasym) ** 2 / (2 * sV ** 2) - W ** 2 / (2 * sW ** 2))
def pops(U, V, W):
    fD = 0.94 * gauss3(U, V, W, 35, 20, 16, -15); fTD = 0.06 * gauss3(U, V, W, 67, 38, 35, -46); fH = 0.0015 * gauss3(U, V, W, 160, 90, 90, -220)
    return fD, fTD, fH
fD, fTD, fH = pops(d.U, d.V, d.W); d['TDD'] = fTD / fD
tD, tTD, tH = pops(*T); print('M3V TD/D', tTD / tD)
sel = {'all_dwarfs': np.ones(len(d), bool),
       'analogues_W_vtot': (np.abs(d.W - T[2]) < 10) & (np.abs(d.vtot - Tv) < 15),
       'analogues_absW_30_65': (d.W.abs() > 30) & (d.W.abs() < 65),
       'analogues_TDD_0.2_1': (d.TDD > 0.2) & (d.TDD < 1.0)}
out = dict(target_UVW_LSR=T.tolist(), target_vtot=Tv, target_TDD=float(tTD / tD), n_total=int(len(d)))
xs = np.linspace(-1.2, 0.8, 401)
fig, ax = plt.subplots(figsize=(7, 4.5))
for lab, k in sel.items():
    f = d.feh[k].values
    res = dict(n=int(k.sum()), median=float(np.median(f)), mean=float(np.mean(f)), sd=float(np.std(f)),
               **{f'frac_gt_{t:+.1f}': float(np.mean(f > t)) for t in (0.0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)})
    # 2-component Gaussian mixture by EM
    w_, m_, s_ = np.array([0.7, 0.3]), np.array([0.0, -0.4]), np.array([0.2, 0.3])
    for _ in range(300):
        pdf = w_ * np.exp(-0.5 * (f[:, None] - m_) ** 2 / s_ ** 2) / s_
        r = pdf / pdf.sum(1, keepdims=True); w_ = r.mean(0); m_ = (r * f[:, None]).sum(0) / r.sum(0); s_ = np.sqrt((r * (f[:, None] - m_) ** 2).sum(0) / r.sum(0))
    res['gmm'] = dict(w=w_.tolist(), mu=m_.tolist(), sd=s_.tolist())
    out[lab] = res
    print(lab, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in res.items() if kk != 'gmm'}, res['gmm'])
    if k.sum() > 30:
        ax.hist(f, bins=np.arange(-1.2, 0.8, 0.05), histtype='step', density=True, label=f'{lab} (n={k.sum()})')
ax.set_xlabel('[Fe/H] (GSP-Spec, dwarfs < 200 pc)'); ax.legend(fontsize=7); ax.axvline(0.3, c='k', ls=':')
plt.tight_layout(); plt.savefig(f'{B}/figs/h10b_feh_prior.png', dpi=110)
json.dump(out, open(f'{B}/out/h10b_feh_prior.json', 'w'), indent=1)
