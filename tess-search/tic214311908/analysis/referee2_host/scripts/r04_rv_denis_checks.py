"""r04: quick independent checks of (a) the Gaia DR3 rv_error percentile of the M3V among comparable faint-star
(rv_method_used = 2) M dwarfs, with r01's local ridge and with/without restricting to the target colour and to
the fraction of deblended transits; (b) the DENIS-2MASS offset of the field around the M3V (host/data/local_*.csv
raw VizieR downloads), recomputed with an independent cross-match; (c) the h15 single-star CMD model's predicted
fraction of single stars at least as over-luminous as the M3V, vs the empirical fraction among all clean field stars.
Output: out/r04_checks.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import norm
from astropy.coordinates import SkyCoord
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
out = {}
# ---- (a) RV
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
F['dMK'] = F.ks_m + 5 * np.log10(F.parallax / 100) - (R['MK'][0] * (F.bp_rp - C0) + R['MK'][1])
T = dict(err=2.21737, grvs=12.5484, sn=14.6958, ntr=21, ndb=12)
c = (F.rv_method_used == 2) & np.isfinite(F.radial_velocity_error) & ((F.grvs_mag - T['grvs']).abs() < 0.25) & ((F.rv_expected_sig_to_noise / T['sn'] - 1).abs() < 0.25) & \
    F.rv_nb_transits.between(15, 30) & (F.ph_qual == 'AAA')
rv = {}
for lab, sel in (('normal RUWE<1.2 |dMK|<0.25', c & (F.ruwe < 1.2) & (F.dMK.abs() < 0.25)),
                 ('normal, |dBPRP|<0.15', c & (F.ruwe < 1.2) & (F.dMK.abs() < 0.25) & ((F.bp_rp - C0).abs() < 0.15)),
                 ('normal, deblended frac>0.4', c & (F.ruwe < 1.2) & (F.dMK.abs() < 0.25) & (F.rv_nb_deblended_transits / F.rv_nb_transits > 0.4)),
                 ('overlum RUWE<1.1 -1.0<dMK<-0.5', c & (F.ruwe < 1.1) & F.dMK.between(-1.0, -0.5)),
                 ('RUWE>2', c & (F.ruwe > 2))):
    e = F.radial_velocity_error[sel].values
    rv[lab] = dict(n=int(len(e)), median=float(np.median(e)), p95=float(np.percentile(e, 95)), target_pct=float(np.mean(e <= T['err'])),
                   frac_gt_p95_normal=None)
    print('RV', lab, rv[lab])
p95n = rv['normal RUWE<1.2 |dMK|<0.25']['p95']
for lab in rv:
    sel_e = None
out['rv'] = rv
# ---- (b) DENIS - 2MASS
D = pd.read_csv(f'{H}/data/local_denis.csv'); Tm = pd.read_csv(f'{H}/data/local_2mass.csv')
cD = SkyCoord(D.RAJ2000.values, D.DEJ2000.values, unit='deg'); cT = SkyCoord(Tm.RAJ2000.values, Tm.DEJ2000.values, unit='deg')
i, d2, _ = cD.match_to_catalog_sky(cT); k = d2.arcsec < 1.0
Dm = D[k].reset_index(drop=True); Tt = Tm.iloc[i[k]].reset_index(drop=True)
tgt = SkyCoord(260.970532, -51.324203, unit='deg'); rt = SkyCoord(Dm.RAJ2000.values, Dm.DEJ2000.values, unit='deg').separation(tgt).arcsec
good = Tt.Qflg.astype(str).str[0].eq('A') & Tt.Qflg.astype(str).str[2].eq('A') & Tt.Cflg.astype(str).eq('000') & (rt > 3)
den = {}
it = int(np.argmin(rt))
den['target'] = dict(sep=float(rt[it]), image=int(Dm.Image[it]), dJ=float(Dm.Jmag[it] - Tt.Jmag[it]), dK=float(Dm.Kmag[it] - Tt.Kmag[it]), DENIS_J=float(Dm.Jmag[it]), DENIS_K=float(Dm.Kmag[it]))
for band in ('J', 'K'):
    dd = (Dm[band + 'mag'] - Tt[band + 'mag']).values
    for lab, sel in (('same image', good & (Dm.Image == den['target']['image'])), ('other images', good & (Dm.Image != den['target']['image']))):
        kk = sel & np.isfinite(dd) & Tt[band + 'mag'].between(9.5, 12.0)
        if kk.sum() >= 3:
            den[f'{band} {lab}'] = dict(n=int(kk.sum()), median=float(np.median(dd[kk])), robust_sd=float(1.4826 * np.median(np.abs(dd[kk] - np.median(dd[kk])))),
                                        images=sorted(set(Dm.Image[kk].astype(int)))[:10])
            print('DENIS', band, lab, den[f'{band} {lab}'])
print('DENIS target', den['target'])
out['denis'] = den
# ---- (c) h15 single-star CMD model tail vs empirical
h10 = json.load(open(f'{H}/out/h10b_feh_prior.json')); g = h10['analogues_absW_30_65']['gmm']
fg = np.linspace(-1.2, 0.8, 2001); pf = sum(w * norm.pdf(fg, m, s) for w, m, s in zip(g['w'], g['mu'], g['sd'])); pf /= np.trapezoid(pf, fg)
a, b, s = -0.2818, -1.8804, np.hypot(0.2905, 0.035)
tx = -0.7606
P_tail = float(np.trapezoid(pf * norm.cdf(tx, a + b * fg, s), fg)); dens = float(np.trapezoid(pf * norm.pdf(tx, a + b * fg, s), fg))
dens_T = float(np.trapezoid(pf * norm.pdf(tx, a + b * fg - 0.7526, s), fg))
# empirical: all clean stars at the target colour
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & ((F.bp_rp - C0).abs() < 0.08)
x = F.dMK[clean].values; txr = r01['target']['cmd']['raw_plx']['dMK']
emp_tail = float(np.mean(x <= txr)); emp_dens = float(np.mean(np.abs(x - txr) < 0.1) / 0.2); emp_dens0 = float(np.mean(np.abs(x - (txr + 0.7526)) < 0.1) / 0.2)
out['h15_cmd_check'] = dict(model_single_frac_as_bright=P_tail, model_single_density_at_target=dens, model_twin_component_density=dens_T,
                            empirical_all_frac_as_bright=emp_tail, empirical_all_density_at_target=emp_dens, empirical_density_at_twin_component_position=emp_dens0,
                            ratio_model_single_to_empirical_all=dens / emp_dens,
                            note='h15 single-star dMKs model = N(a + b f, 0.29) integrated over the kinematic [Fe/H] prior; empirical = all clean field stars (singles + binaries) at the M3V colour')
print('h15 CMD check', out['h15_cmd_check'])
# CMD-only P(twin) with empirical densities: prior twin (q>=0.8, surviving) 0.0234 of stars; twin density = 0.0234 * empirical density of components at x+0.753
pT = 0.0234
num = pT * emp_dens0; out['cmd_only_P_twin_empirical'] = float(num / emp_dens)
out['cmd_only_P_twin_h15_style'] = float(pT * dens_T / (pT * dens_T + (1 - 0.268 * 0.39) * dens))
print('CMD-only P(twin): empirical-density version', out['cmd_only_P_twin_empirical'], ' h15-style', out['cmd_only_P_twin_h15_style'])
json.dump(out, open(f'{B}/out/r04_checks.json', 'w'), indent=1)
