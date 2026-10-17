"""h06b: which near-IR photometry of the M3V is right? Local (8') empirical checks:
 (1) DENIS - 2MASS for the same stars (all, and in the DENIS image 1724791 that measured the M3V), vs magnitude;
 (2) unWISE - AllWISE W1/W2 for isolated stars (puts unWISE on the AllWISE scale used by the comparison samples);
 (3) VHS DR5 saturated-star behaviour: VHS aperture mags (ap3, ap6, pmag) vs 2MASS for stars of similar 2MASS
     brightness; inverse-calibrate the M3V's saturated VHS ap3/ap6 to an independent J and Ks at the VHS epoch (2014);
 (4) DENIS I vs Gaia-predicted I for local stars.
Output: out/h06b_nir_local.json, figs/h06b_nir_local.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
B = 'tic214311908/analysis/host'
D = pd.read_csv(f'{B}/data/local_denis.csv'); T = pd.read_csv(f'{B}/data/local_2mass.csv')
UW = pd.read_csv(f'{B}/data/local_unwise.csv'); AW = pd.read_csv(f'{B}/data/local_allwise.csv')
V = pd.read_csv(f'{B}/data/local_vhs.csv'); G = pd.read_csv(f'{B}/data/local_gaia.csv')
tgt = SkyCoord(260.970532, -51.324203, unit='deg')  # 2MASS position of the M3V (epoch 1999.5)
mira = SkyCoord(260.977681, -51.323727, unit='deg')
def xm(a, b, ra1, de1, ra2, de2, r=1.5):
    c1 = SkyCoord(a[ra1].values, a[de1].values, unit='deg'); c2 = SkyCoord(b[ra2].values, b[de2].values, unit='deg')
    i, d, _ = c1.match_to_catalog_sky(c2); k = d.arcsec < r
    return a[k].reset_index(drop=True), b.iloc[i[k]].reset_index(drop=True), d.arcsec[k]
out = {}
# ---- (1) DENIS vs 2MASS
d, t, sep = xm(D, T, 'RAJ2000', 'DEJ2000', 'RAJ2000', 'DEJ2000', 1.5)
cD = SkyCoord(d.RAJ2000.values, d.DEJ2000.values, unit='deg')
d['r_tgt'] = cD.separation(tgt).arcsec; d['r_mira'] = cD.separation(mira).arcsec
good = t.Qflg.astype(str).str[0].isin(['A']) & t.Qflg.astype(str).str[2].isin(['A']) & (t.Cflg.astype(str) == '000') & (d.r_tgt > 3)
res = {}
for band, dj, tj in (('J', 'Jmag', 'Jmag'), ('K', 'Kmag', 'Kmag')):
    dd = (d[dj] - t[tj]).values
    for lab, sel in (('all', good), ('same_image', good & (d.Image == 1724791)), ('within_120', good & (d.r_tgt < 120))):
        for lo, hi in ((9.5, 12.0), (10.0, 11.5), (12.0, 13.5)):
            k = sel & (t[tj] > lo) & (t[tj] < hi) & np.isfinite(dd)
            if k.sum() < 3: continue
            res[f'{band}_{lab}_{lo}-{hi}'] = dict(n=int(k.sum()), median=float(np.median(dd[k])), robust_sd=float(1.4826 * np.median(np.abs(dd[k] - np.median(dd[k])))),
                                                frac_ge_0p30=float(np.mean(dd[k] >= 0.30)))
    # stars near the Mira (within 30") other than the target
    k = (d.r_mira < 40) & (d.r_tgt > 3) & np.isfinite(dd)
    res[f'{band}_near_mira_lt40'] = [dict(r_mira=float(a), dmag=float(b), m2mass=float(c)) for a, b, c in zip(d.r_mira[k], dd[k], t[tj][k])]
it = np.argmin(d.r_tgt.values)
res['target'] = dict(dJ=float(d.Jmag[it] - t.Jmag[it]), dK=float(d.Kmag[it] - t.Kmag[it]), image=int(d.Image[it]), r=float(d.r_tgt[it]))
out['denis_minus_2mass'] = res
print('DENIS-2MASS:', json.dumps({k: v for k, v in res.items() if 'near_mira' not in k}, indent=0)[:1800])
print('near Mira:', res['J_near_mira_lt40'], res['K_near_mira_lt40'])
# ---- (2) unWISE vs AllWISE
UW['W1u'] = 22.5 - 2.5 * np.log10(UW.FW1.clip(lower=1e-3)); UW['W2u'] = 22.5 - 2.5 * np.log10(UW.FW2.clip(lower=1e-3))
u_, a_, sep = xm(UW, AW, 'RAJ2000', 'DEJ2000', 'RAJ2000', 'DEJ2000', 1.5)
cU = SkyCoord(u_.RAJ2000.values, u_.DEJ2000.values, unit='deg')
u_['r_tgt'] = cU.separation(tgt).arcsec; u_['r_mira'] = cU.separation(mira).arcsec
ok = (a_.ccf.astype(str) == '0000') & (a_.ex == 0) & (u_.fFW1 > 0.95) & (u_.r_mira > 60) & (u_.FlagsW1 == 0)
rw = {}
for b, ua, aa in (('W1', 'W1u', 'W1mag'), ('W2', 'W2u', 'W2mag')):
    dd = (u_[ua] - a_[aa]).values
    for lo, hi in ((8.5, 11.5), (9.5, 11.0), (11.5, 13.0)):
        k = ok & (a_[aa] > lo) & (a_[aa] < hi) & np.isfinite(dd)
        if k.sum() >= 3:
            rw[f'{b}_{lo}-{hi}'] = dict(n=int(k.sum()), median=float(np.median(dd[k])), robust_sd=float(1.4826 * np.median(np.abs(dd[k] - np.median(dd[k])))))
out['unwise_minus_allwise'] = rw
print('unWISE-AllWISE:', rw)
# unWISE fracflux for stars near bright stars: how often is fFW1 < 0.9
# ---- (3) VHS saturation calibration against 2MASS
v_, t2, sep = xm(V, T, 'RAJ2000', 'DEJ2000', 'RAJ2000', 'DEJ2000', 1.0)
cV = SkyCoord(v_.RAJ2000.values, v_.DEJ2000.values, unit='deg')
v_['r_tgt'] = cV.separation(tgt).arcsec; v_['r_mira'] = cV.separation(mira).arcsec
okv = (t2.Qflg.astype(str).str[0] == 'A') & (t2.Qflg.astype(str).str[2] == 'A') & (t2.Cflg.astype(str) == '000') & (v_.r_tgt > 3) & (v_.r_mira > 40) & (v_.PriOrSec.fillna(0) == 0) & (t2.prox > 6)
# the M3V moved ~2.1" between 2MASS (1999.5) and VHS (2014.6): take its VHS row by the 2014.6 position
cVall = SkyCoord(V.RAJ2000.values, V.DEJ2000.values, unit='deg')
ra14 = 260.969674 + (-113.671152) * (2014.58 - 2016.0) / 3.6e6 / np.cos(np.radians(-51.323821)); de14 = -51.323821 + 85.54355 * (2014.58 - 2016.0) / 3.6e6
itv = np.argmin(cVall.separation(SkyCoord(ra14, de14, unit='deg')).arcsec)
vt = V.iloc[itv].copy(); vt['r_tgt'] = cVall.separation(SkyCoord(ra14, de14, unit='deg')).arcsec[itv]
print('target VHS row: r=%.2f Jap3 %.3f Jap6 %.3f Jpmag %.3f Ksap3 %.3f Ksap6 %.3f Kspmag %.3f' % (vt.r_tgt, vt.Jap3, vt.Jap6, vt.Jpmag, vt.Ksap3, vt.Ksap6, vt.Kspmag))
vres = {}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for j, (band, tm, cols) in enumerate((('J', 'Jmag', ['Jap3', 'Jap6', 'Jpmag']), ('Ks', 'Kmag', ['Ksap3', 'Ksap6', 'Kspmag']))):
    for col in cols:
        x = t2[tm].values; y = v_[col].values; k = okv & np.isfinite(x) & np.isfinite(y) & (x > 9.0) & (x < 13.5)
        # inverse: for stars whose VHS value is within +-0.15 of the target's, distribution of their 2MASS mag
        yt = vt[col]
        nb = k & (np.abs(y - yt) < 0.15)
        # also forward: at 2MASS mag within +-0.1 of the target's 2MASS value, distribution of VHS value
        m2t = {'J': 11.363, 'Ks': 10.482}[band]
        fw = k & (np.abs(x - m2t) < 0.12)
        # local linear fit VHS = a + b * 2MASS over 2MASS within +-0.6 mag of the target, then invert
        loc = k & (np.abs(x - m2t) < 0.6)
        cf = np.polyfit(x[loc], y[loc], 1); rr = y[loc] - np.polyval(cf, x[loc]); sdl = 1.4826 * np.median(np.abs(rr - np.median(rr)))
        inv = (yt - cf[1]) / cf[0]; inv_err = sdl / abs(cf[0])
        vres[f'{band}_{col}_inverse_fit'] = dict(n=int(loc.sum()), slope=float(cf[0]), scatter=float(sdl), true_mag=float(inv), err=float(inv_err))
        vres[f'{band}_{col}'] = dict(target_vhs=float(yt), n_inverse=int(nb.sum()),
                                     inverse_2mass_median=float(np.median(x[nb])) if nb.sum() else None,
                                     inverse_2mass_robust_sd=float(1.4826 * np.median(np.abs(x[nb] - np.median(x[nb])))) if nb.sum() > 2 else None,
                                     n_forward=int(fw.sum()), forward_vhs_median=float(np.median(y[fw])) if fw.sum() else None,
                                     forward_vhs_robust_sd=float(1.4826 * np.median(np.abs(y[fw] - np.median(y[fw])))) if fw.sum() > 2 else None,
                                     forward_target_pct=float(np.mean(y[fw] <= yt)) if fw.sum() else None)
        if col.endswith('ap3'):
            ax[j].scatter(x[k], y[k] - x[k], s=4, alpha=0.5, label='local stars')
            ax[j].scatter([m2t], [yt - m2t], c='r', s=60, marker='*', label='M3V (2MASS)')
            ax[j].scatter([{'J': 11.774, 'Ks': 10.823}[band]], [yt - {'J': 11.774, 'Ks': 10.823}[band]], c='m', s=60, marker='s', label='M3V (DENIS)')
            ax[j].set_xlabel(f'2MASS {band}'); ax[j].set_ylabel(f'VHS {col} - 2MASS'); ax[j].set_ylim(-0.6, 1.6); ax[j].legend(fontsize=7)
out['vhs_saturation'] = vres
for k_, v__ in vres.items():
    print(k_, v__)
# ---- (4) DENIS I vs Gaia-predicted I_C (Riello+21 G-I_C polynomial)
g_, d2, sep = xm(G, D, 'RAJ2000', 'DEJ2000', 'RAJ2000', 'DEJ2000', 1.5)
x = g_['BP-RP'].values
GI = 0.02085 + 0.7419 * x - 0.09631 * x ** 2
dI = d2.Imag.values - (g_.Gmag.values - GI)
k = np.isfinite(dI) & (x > 1.5) & (x < 2.75) & (d2.Imag > 10) & (d2.Imag < 13.5)
cG = SkyCoord(g_.RAJ2000.values, g_.DEJ2000.values, unit='deg')
rt = cG.separation(tgt).arcsec
ii = np.argmin(rt)
out['denis_I_minus_gaiaIc'] = dict(n=int((k & (rt > 3)).sum()), median=float(np.median(dI[k & (rt > 3)])), robust_sd=float(1.4826 * np.median(np.abs(dI[k & (rt > 3)] - np.median(dI[k & (rt > 3)])))),
                                   target=float(dI[ii]), note='Riello+21 G-I_C relation valid to BP-RP 2.75; local stars 1.5<BP-RP<2.75, 10<I<13.5')
print('DENIS I - Gaia I_C:', out['denis_I_minus_gaiaIc'])
ax[2].scatter(x[k], dI[k], s=4, alpha=0.5); ax[2].scatter([x[ii]], [dI[ii]], c='r', marker='*', s=80); ax[2].set_xlabel('BP-RP'); ax[2].set_ylabel('DENIS I - Gaia-predicted I_C')
plt.tight_layout(); plt.savefig(f'{B}/figs/h06b_nir_local.png', dpi=110)
json.dump(out, open(f'{B}/out/h06b_nir_local.json', 'w'), indent=1)
