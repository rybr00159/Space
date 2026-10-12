"""Check 6: statistical independence and the quoted significances.

(a) My pixel test (r3): per-sector z_T (H_target rejected) for the a/b pair (primary, a priori) and the other pairs;
    Stouffer combination; empirical combined null built by drawing one null stack per sector independently;
    sensitivity of the combined z to a 20-30 % error in the null sigma (finite number of independent fake windows).
(b) Effective number of independent fake / control windows in each analysis (overlapping windows share cadences):
    N_eff per transit epoch ~ 1 + (phase range) / (window length in phase).
(c) The quoted significances of the other analyses, recomputed from their own JSON numbers:
    loc_independent z_T (Gaussian approx.), psf_photometry delta chi2 and f_target with systematics.
"""
import json
import numpy as np
from rcommon import OUT, P, T14, dump

RNG = np.random.default_rng(7)
res = {}
for tag in ('none_plane', 'mad10_plane', 'mad10_median'):
    r = json.load(open(f'{OUT}/r3_pixel_ratio_{tag}.json'))
    nul = np.load(f'{OUT}/r3_nulls_{tag}.npz') if tag != 'mad10_median' else None
    out = {}
    for pair in ('a/b', 'a/t', 'Ga/Gb', 'Ga/t', 'n/t'):
        secs = [s for s in r if pair in r[s]['pairs']]
        z = np.array([r[s]['pairs'][pair]['z_T_models'] for s in secs])
        zm = np.array([r[s]['pairs'][pair]['z_T_models_mad'] for s in secs])
        zst = float(z.sum() / np.sqrt(len(z)))
        d = dict(sectors=secs, z_T_per_sector=z.tolist(), stouffer=zst, stouffer_mad=float(zm.sum() / np.sqrt(len(zm))),
                 stouffer_if_null_sigma_x1p2=float((z / 1.2).sum() / np.sqrt(len(z))),
                 stouffer_if_null_sigma_x1p3=float((z / 1.3).sum() / np.sqrt(len(z))),
                 Q_obs=[r[s]['pairs'][pair]['Q'] for s in secs],
                 Q_T_max=[r[s]['pairs'][pair]['Q_T_max_models'] for s in secs],
                 Q_N_range=[r[s]['pairs'][pair]['Q_N_range_models'] for s in secs],
                 z_N_range=[r[s]['pairs'][pair]['z_N_models'] for s in secs])
        if nul is not None:
            keys = [f"S{s}_{pair.split('/')[0]}_{pair.split('/')[1]}_zT_null" for s in secs]
            arrs = [nul[k] for k in keys]
            comb = np.sum([a[RNG.integers(0, len(a), 200000)] for a in arrs], 0) / np.sqrt(len(arrs))
            d['empirical_combined_null'] = dict(mean=float(comb.mean()), std=float(comb.std()), max=float(comb.max()),
                                                frac_ge_obs=float(np.mean(comb >= zst)), n=int(len(comb)))
        out[pair] = d
        print(f"{tag:13s} {pair:6s} z_T per sector {np.round(z,1)} Stouffer {zst:.1f} (mad {d['stouffer_mad']:.1f}; "
              f"null sigma x1.3 -> {d['stouffer_if_null_sigma_x1p3']:.1f})"
              + (f"; empirical combined null max {d['empirical_combined_null']['max']:.1f}, frac>=obs {d['empirical_combined_null']['frac_ge_obs']:.1e}"
                 if 'empirical_combined_null' in d else ''))
    res[f'pixel_test_{tag}'] = out
# also the model-free (kappa) bound for n/t
r = json.load(open(f'{OUT}/r3_pixel_ratio_none_plane.json'))
zl = [r[s]['pairs']['n/t']['z_L'] for s in r if 'n/t' in r[s]['pairs']]
res['pixel_test_none_plane']['n/t']['model_free_kappa_z'] = zl
res['pixel_test_none_plane']['n/t']['model_free_kappa_stouffer'] = float(np.sum(zl) / np.sqrt(len(zl)))
print('model-free (envelope kappa) n/t z:', np.round(zl, 2), 'Stouffer %.1f' % (np.sum(zl) / np.sqrt(len(zl))))

# (b) effective numbers of independent windows
def neff(phase_lo, phase_hi, win_days, excl=0.0):
    rng_ = (phase_hi - phase_lo) - excl
    return 1 + rng_ / (win_days / P)
neff_tab = {
    'loc_independent fake windows (phi 0.25-0.75 step 0.05, window +-2.5 T14)':
        dict(per_epoch_nominal=11, per_epoch_eff=neff(0.25, 0.75, 5 * T14)),
    'this audit fake windows (phi 0.20-0.80 step 0.04, window +-2.0 T14)':
        dict(per_epoch_nominal=16, per_epoch_eff=neff(0.20, 0.80, 4 * T14)),
    'psf_photometry control phases (step 0.01, window +-0.45 d, all epochs fitted jointly)':
        dict(n_nominal=65, n_eff_joint=neff(0.1435, 0.8565, 0.90, excl=0.07)),
}
for k, v in neff_tab.items():
    print(k, {kk: round(vv, 1) for kk, vv in v.items()})
res['effective_window_counts'] = neff_tab

# (c) quoted significances of the other analyses
li = json.load(open('/home/user/work/tic87/loc_independent/combine_prf_main__prf_prf.json'))['combined']
pp = json.load(open('/home/user/work/tic87/psf_photometry/results.json'))['nominal']['combined']
f = pp['f_target']
f_best, f_lo, f_hi = f['f_best'], f['f_lo'], f['f_hi']
stat = (f_hi - f_lo) / 2; sys_ = 0.07
quoted = dict(
    loc_independent=dict(s_par_real=li['s_par_real'], null_T_mean=li['s_par_T_mean'], null_T_std=li['s_par_T_std'],
                         z_T_gauss=li['z_T'], frac_T_sims_ge_real=li['frac_T_sims_ge_real'],
                         z_T_if_null_sigma_x1p3=(li['s_par_real'] - li['s_par_T_mean']) / (1.3 * li['s_par_T_std']),
                         note='0 of 20000 sims >= real: empirical p < 5e-5 (3.9 sigma) by count; >8 sigma relies on the '
                              'Gaussian shape of the null far beyond the sampled range'),
    psf_photometry=dict(chi2_f1=f['chi2_f1'], chi2_f0=f['chi2_f0'], delta_chi2=pp['delta_chi2_target_minus_neighbour'],
                        f_target=f_best, f_stat=stat, f_sys=sys_,
                        z_from_f1_stat_only=(1 - f_best) / stat, z_from_f1_with_sys=(1 - f_best) / np.hypot(stat, sys_),
                        note='delta chi2 = 129 is a stat-only likelihood ratio (no PRF/registration systematic); '
                             'with the adopted 0.07 systematic f_target = 1 is excluded at ~7.4 sigma'))
for k, v in quoted.items():
    print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != 'note'})
res['quoted_significances_recomputed'] = quoted
dump(res, f'{OUT}/r6_stats.json')
