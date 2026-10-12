"""Assemble results.json from all intermediate products (no new fitting)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import chi2 as chi2dist, norm
from common import HERE, SECTORS, to_json

L = lambda f: json.load(open(f'{HERE}/{f}'))
reg = L('registration.json')
prf = L('diff_prf_main.json'); gau = L('diff_gauss_main.json')
jp = L('joint_prf_main.json'); jg = L('joint_gauss_main.json')
combos = {k: L(f'combine_{k}.json') for k in ['prf_main__prf_prf', 'gauss_main__gauss_gauss', 'gauss_main__prf_gauss', 'prf_main__gauss_prf']}
inj = {k: L(f'inj_{k}.json') for k in ['prf_prf', 'gauss_gauss', 'prf_gauss', 'gauss_prf']}
jit = L('jitter_pertransit.json'); cen = L('centroid_shift.json'); var = L('variants_table.json')
inj3 = L('inj3_prf.json'); dep = L('depth_consistency.json')

out = dict(task='Independent pixel-level localisation of the P=4.136 d transit near TIC 87055304',
           folder=HERE, ephemeris=dict(P=4.136195, T0=1658.41325, T14=0.1875),
           neighbour='Gaia DR3 6720986208895208960 = TIC 87055315 (T=13.54)',
           target='Gaia DR3 6720986204578421888 = TIC 87055304 (T=11.74)')
# geometry
out['geometry'] = dict(neighbour_offset_from_target_EN_arcsec_2016=[6.594, 6.600], separation_2016_arcsec=9.33,
                       separation_at_S66_epoch_arcsec=float(np.hypot(*prf['66']['neigh_EN'])),
                       note='8.94" in CONTEXT.md is the separation from the TIC (J2000-epoch) position, not star-to-star')
# registration
out['registration'] = {s: dict(prf_local_r4_offset_px=[reg[s]['variants']['prf:local_r4']['dx'], reg[s]['variants']['prf:local_r4']['dy']],
                               prf_local_r4_offset_EN_arcsec=[reg[s]['variants']['prf:local_r4']['dE'], reg[s]['variants']['prf:local_r4']['dN']],
                               spread_over_prf_variants_px=float(np.ptp([[v['dx'], v['dy']] for k, v in reg[s]['variants'].items() if k.startswith('prf')], axis=0).max()),
                               wcs_check_star_rms_EN_arcsec=reg[s].get('wcs_check_rms_EN_arcsec'),
                               gauss_sigma_px=reg[s]['gauss_sigma_xy'],
                               resid_rms_frac_r2_prf=reg[s]['prf_resid_rms_frac_r2'], resid_rms_frac_r2_gauss=reg[s]['gauss_resid_rms_frac_r2'],
                               fitted_N_over_T_flux=reg[s]['variants']['prf:local_r4']['ratio_NT'],
                               catalog_N_over_T=reg[s]['neigh_target_ratio_catalog'],
                               blend_photocentre_s_par_prf=reg[s]['prf_blend_photocentre']['s_par'],
                               blend_photocentre_s_par_expected=reg[s]['prf_blend_photocentre']['s_par_expected_catalog'])
                       for s in map(str, SECTORS)}
# per-sector
per = {}
for s in map(str, SECTORS):
    d = {}
    for name, r in [('prf', prf), ('gauss', gau)]:
        st = r[s]['stacked']; b = r[s]['bootstrap']
        d[name] = dict(n_transits=r[s]['n_transits'], offT_EN=[st['offT_E'], st['offT_N']],
                       offT_EN_boot_sd=[b['offT_E']['std'], b['offT_N']['std']], offN_EN=[st['offN_E'], st['offN_N']],
                       dist_target=st['offT'], dist_neigh=st['offN'],
                       s_par=st['s_par'], s_par_boot_p16_p84=[b['s_par']['p16'], b['s_par']['p84']], perp_arcsec=st['perp_arcsec'],
                       f_N_2src=st['f_N'], f_N_boot_sd=b['f_N']['std'], a_T=st['a_T'], a_N=st['a_N'], sa_T=st['sa_T'], sa_N=st['sa_N'],
                       dchi2_T_minus_N=st['dchi2_TminusN'], dchi2_C3_minus_N=st['dchi2_C3minusN'],
                       TN3_amp=st['TN3_amp'], TN3_amp_err=st['TN3_amp_err'], amp_e_per_s=st['free_amp'],
                       free_chi2_per_dof=st['free_chi2_per_dof'],
                       per_transit_s_par=[p['s_par'] for p in r[s]['per_transit']])
    for k, c in combos.items():
        o = c['sectors'][s]
        d[f'cal_{k}'] = dict(s_par_H_T=[o['s_par']['T_mean'], o['s_par']['T_std']], s_par_H_N=[o['s_par']['N_mean'], o['s_par']['N_std']],
                             z_T=o['s_par']['z_T'], z_N=o['s_par']['z_N'], p_emp_T=o['s_par']['p_emp_T'], p_emp_N=o['s_par']['p_emp_N'],
                             lnLR_TN_spar=o['s_par']['lnLR_TN_gauss'], lnLR_TN_EN2d=o['EN2d']['lnLR_TN'], lnLR_TN_dchi2=o['dchi2']['lnLR_TN_gauss'])
    d['centroid_shift_r1.5'] = cen[s]['r1.5']
    d['three_host_EN_mahal'] = {h: inj3[s][h]['mahal'] for h in ['T', 'N', 'C3']}
    d['pointing'] = dict(leverage_FT_over_A=jit[s]['leverage_FT_over_A'],
                         spurious_s_par_per_transit=[t['spurious_s_par'] for t in jit[s]['transits']],
                         single_transit_P_spar_gt_half=[jit[s]['single_transit_cal']['P_spar_gt_half_T'], jit[s]['single_transit_cal']['P_spar_gt_half_N']])
    d['depth'] = dep[s]
    per[s] = d
out['per_sector'] = per
# combined
comb = {}
for k, c in combos.items():
    cc = c['combined']
    # Fisher combination of empirical per-sector p-values under H_T (floored at 1/(n+1))
    ps = []
    for s in map(str, SECTORS):
        o = c['sectors'][s]['s_par']
        ps.append(max(o['p_emp_T'], 1.0 / (o['nsim'] + 1)))
    X = -2 * np.sum(np.log(ps)); pF = chi2dist.sf(X, 2 * len(ps))
    comb[k] = dict(s_par_real=cc['s_par_real'], H_T=[cc['s_par_T_mean'], cc['s_par_T_std']], H_N=[cc['s_par_N_mean'], cc['s_par_N_std']],
                   z_T=cc['z_T'], z_N=cc['z_N'], frac_T_sims_ge_real=cc['frac_T_sims_ge_real'], frac_N_sims_le_real=cc['frac_N_sims_le_real'],
                   lnLR_TN_combined_spar=cc['lnLR_TN_combined_spar_gauss'], lnLR_TN_sum_EN2d=cc['sum_lnLR_EN'],
                   lnLR_TN_sum_spar=cc['sum_lnLR_spar'], dchi2_sum_real=cc['dchi2_sum_real'], z_dchi2_T=cc['z_dchi2_T'], z_dchi2_N=cc['z_dchi2_N'],
                   fisher_pT_per_sector=ps, fisher_p_T=float(pF), fisher_sigma_T=float(norm.isf(pF / 2)))
out['combined_injection_calibrated'] = comb
# combined f_T (2-source) with calibrated and bootstrap errors
for name, r, ik in [('prf', prf, 'prf_prf'), ('gauss', gau, 'gauss_gauss')]:
    fT = np.array([1 - r[s]['stacked']['f_N'] for s in map(str, SECTORS)])
    sb = np.array([r[s]['bootstrap']['f_N']['std'] for s in map(str, SECTORS)])
    sc = np.array([inj[ik][s]['summary']['resampled_N']['f_N']['std'] for s in map(str, SECTORS)])
    wb = 1 / sb ** 2; wc = 1 / sc ** 2
    out[f'f_target_2source_{name}'] = dict(per_sector=fT.tolist(), boot_sd=sb.tolist(), cal_sd=sc.tolist(),
                                           comb_boot=[float((wb * fT).sum() / wb.sum()), float(1 / np.sqrt(wb.sum()))],
                                           comb_cal=[float((wc * fT).sum() / wc.sum()), float(1 / np.sqrt(wc.sum()))])
out['joint_sky_fit'] = {n: dict(best_EN=j['best_EN'], boot_sd_EN=j['boot_sd_EN'], dist_target=j['offset_from_target'],
                                dist_neigh=j['offset_from_neigh'], s_par=j['s_par'], s_par_boot_sd=j['s_par_boot_sd'],
                                perp=j['perp'], perp_boot_sd=j['perp_boot_sd'], dchi2_target=j['dchi2_target'],
                                dchi2_neigh=j['dchi2_neigh'], dchi2_c3=j['dchi2_c3'], mahal_target=j['mahal_target_boot'],
                                mahal_neigh=j['mahal_neigh_boot'], mahal_c3=j['mahal_c3_boot'])
                        for n, j in [('prf', jp), ('gauss', jg)]}
out['three_host_calibration_prf'] = dict(sum_lnL=inj3['sum_lnL'], lnLR_N_over_T=inj3['lnLR_N_over_T'], lnLR_N_over_C3=inj3['lnLR_N_over_C3'])
out['sign_test_per_transit'] = dict(n_transits=int(sum(jit[s]['n_real'] for s in map(str, SECTORS))),
                                    n_s_par_gt_half=int(sum(jit[s]['n_real_gt_half'] for s in map(str, SECTORS))),
                                    log10_P_H_T=jit['combined']['log10_P_all_T'], log10_P_H_N=jit['combined']['log10_P_all_N'])
# conservative sign test: per-sector 95% upper limit on P(s_par>0.5 | H_T), effective n = n_fake/2
from scipy.stats import beta as betad
lp = 0.0
cons = {}
for s_ in map(str, SECTORS):
    c_ = jit[s_]['single_transit_cal']; neff = c_['n_fake'] / 2.0; keff = c_['P_spar_gt_half_T'] * neff
    pu = float(betad.ppf(0.95, keff + 1, neff - keff))
    cons[s_] = pu
    lp += jit[s_]['n_real_gt_half'] * np.log10(pu)
out['sign_test_per_transit']['conservative_pT_upper95'] = cons
out['sign_test_per_transit']['log10_P_H_T_conservative'] = float(lp)
out['variants'] = {k: dict(desc=v['desc'], s_par_comb=v['s_par_comb'], s_par_comb_err_boot=v['s_par_comb_err'],
                           per_sector=[v['sectors'][str(s)]['s_par'] for s in SECTORS], sum_dchi2_T_minus_N=v['dchi2_sum'])
                   for k, v in var.items()}
to_json(out, f'{HERE}/results.json')
for k, v in comb.items():
    print(k, 'z_T=%.1f z_N=%.1f lnLR=%.1f fisher p_T=%.2e (%.1f sigma)' % (v['z_T'], v['z_N'], v['lnLR_TN_combined_spar'], v['fisher_p_T'], v['fisher_sigma_T']))
for name in ['prf', 'gauss']:
    print('f_T', name, out[f'f_target_2source_{name}']['comb_boot'], out[f'f_target_2source_{name}']['comb_cal'])
print('sign test', out['sign_test_per_transit'])
