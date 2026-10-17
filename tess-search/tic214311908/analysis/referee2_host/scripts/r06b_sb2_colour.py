"""r06b: colour residual r(G-Ks) and CMD offset dMKs of Gaia DR3 double-lined spectroscopic binaries (SB2, SB2C) and
of SB1/astrometric orbits among nearby M dwarfs, relative to GLOBAL ridges built here from the clean field sample
(binned KDE mode of M_Ks and median G-Ks in 0.05-mag BP-RP bins, BP-RP 2.30-3.25, linear interpolation).
An SB2 with q near 1 is a known unresolved near-equal pair: the twin discriminant predicts r(G-Ks) ~ normal-star value
(~0 relative to the ridge) at dMKs ~ -0.5 to -0.75.
Output: out/r06b_sb2_colour.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import gaussian_kde
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
def Cs(d):
    x = d.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    return (d.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * d.phot_g_mean_mag ** 7.618399)
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv'); F['Cs'] = Cs(F)
F['MK'] = F.ks_m + 5 * np.log10(F.parallax / 100); F['GK'] = F.phot_g_mean_mag - F.ks_m
cl = F[(F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3)]
ce = np.arange(2.30, 3.2501, 0.05); cc, mk, gk, rgk_sd = [], [], [], []
for lo in ce[:-1]:
    d = cl[cl.bp_rp.between(lo, lo + 0.05)]
    if len(d) < 50: continue
    k = gaussian_kde(d.MK, bw_method=0.15); g = np.linspace(d.MK.quantile(0.02), d.MK.quantile(0.98), 1500)
    cc.append(np.median(d.bp_rp)); mk.append(g[np.argmax(k(g))]); gk.append(np.median(d.GK[(d.MK - g[np.argmax(k(g))]).abs() < 0.2]))
cc, mk, gk = map(np.array, (cc, mk, gk))
S = pd.read_csv(f'{B}/data/sb2/nss_mdwarfs.csv'); S['Cs'] = Cs(S)
S = S[S.bp_rp.between(cc[0], cc[-1]) & (S.ph_qual.astype(str).str[2] == 'A')].copy()
S['dMK'] = S.ks_m + 5 * np.log10(S.parallax / 100) - np.interp(S.bp_rp, cc, mk)
S['rGK'] = S.phot_g_mean_mag - S.ks_m - np.interp(S.bp_rp, cc, gk)
# reference: clean stars with the same treatment
cl = cl.copy(); cl['dMK'] = cl.MK - np.interp(cl.bp_rp, cc, mk); cl['rGK'] = cl.GK - np.interp(cl.bp_rp, cc, gk)
out = dict(ridge=dict(bp_rp=cc.tolist(), MK=mk.tolist(), GK=gk.tolist()))
sb2 = S[S.nss_solution_type.isin(['SB2', 'SB2C'])].drop_duplicates('source_id')
rows = []
for _, r in sb2.iterrows():
    k = cl.dMK.between(r.dMK - 0.1, r.dMK + 0.1) & cl.bp_rp.between(r.bp_rp - 0.15, r.bp_rp + 0.15)
    rows.append(dict(source_id=int(r.source_id), type=r.nss_solution_type, P=float(r.period), q=float(r.mass_ratio) if np.isfinite(r.mass_ratio) else None,
                     bp_rp=float(r.bp_rp), dMK=float(r.dMK), rGK=float(r.rGK), Cs=float(r.Cs), ruwe=float(r.ruwe), ks_err=float(r.ks_msigcom),
                     clean_same_dMK_median_rGK=float(np.median(cl.rGK[k])) if k.sum() else None, clean_n=int(k.sum())))
    print(rows[-1])
out['sb2'] = rows
for typ in (['SB1', 'AstroSpectroSB1'], ['Orbital', 'OrbitalTargetedSearch', 'OrbitalTargetedSearchValidated']):
    d = S[S.nss_solution_type.isin(typ)].drop_duplicates('source_id')
    res = {}
    for lo, hi in ((-1.2, -0.55), (-0.55, -0.25), (-0.25, 0.25)):
        k = d.dMK.between(lo, hi); kc = cl.dMK.between(lo, hi)
        res[f'{lo}..{hi}'] = dict(n=int(k.sum()), median_rGK=float(np.median(d.rGK[k])) if k.sum() else None, clean_median_rGK=float(np.median(cl.rGK[kc])))
    out['/'.join(typ)] = res; print(typ, res)
sv = sb2[(sb2.mass_ratio > 0.8)]
out['sb2_q_gt_0.8_summary'] = dict(n=int(len(sv)), median_dMK=float(np.median(sv.dMK)) if len(sv) else None, median_rGK=float(np.median(sv.rGK)) if len(sv) else None)
print(out['sb2_q_gt_0.8_summary'])
json.dump(out, open(f'{B}/out/r06b_sb2_colour.json', 'w'), indent=1)
