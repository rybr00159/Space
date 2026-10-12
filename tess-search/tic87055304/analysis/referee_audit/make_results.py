"""Collect the audit numbers into results.json (all values read from the per-check JSON files)."""
import json
import numpy as np
from rcommon import OUT, SECTORS, dump

L = lambda f: json.load(open(f'{OUT}/{f}'))
r1, r2, r4, r5, r6 = L('r1_wcs_check.json'), L('r2_registration.json'), L('r4_timing.json'), L('r5_jitter.json'), L('r6_stats.json')
r3 = L('r3_pixel_ratio_none_plane.json'); r3c = L('r3c_flux_inequality_none.json'); r3cm = L('r3c_flux_inequality_mad10.json')
b1 = L('b1_centroid_audit.json')

S = [str(s) for s in SECTORS]
res = dict(
    folder=OUT,
    verdict_task_A=('No shared systematic found that could make a target-hosted signal look neighbour-hosted. The localisation '
                    'conclusion (host = TIC 87055315 or a source within ~1" of it) survives; a PRF-free pixel test written for '
                    'this audit independently rejects a target host. Verdict scale: supports false positive (for a planet on '
                    'TIC 87055304). Some quoted significances are overstated (see statistics).'),
    verdict_task_B=('src/centroid.py could not exclude the neighbour: it is biased toward the aperture centre (positive-clipped '
                    'flux-weighted centroid in an r = 3 px circle on the target), its S66 stack is dominated by two scattered-light '
                    'windows, and its own sigma includes a 0.15 px floor. The ">6 sigma" exclusion divided the target-neighbour '
                    'distance by the floor-free per-axis bootstrap error (err_px) and assumed an unbiased estimator; both are wrong.'),
    geometry={s: dict(cam_ccd=r1[s]['header_cam_ccd'], det_sign=r1[s]['det_sign'], PA_plus_x_deg=r1[s]['PA_of_+x_deg'],
                      PA_plus_y_deg=r1[s]['PA_of_+y_deg'], nb_vec_px=r1[s]['nb_vec_px'], nb_sep_arcsec=r1[s]['nb_sep_arcsec'],
                      nb_PA_deg=r1[s]['nb_PA_deg'], tess_point_minus_wcs_vec_px=r1[s]['tess_point_minus_wcs_vec_px'],
                      fieldstar_median_dist_px={k: v['median_dist_px'] for k, v in r1[s]['field_star_test'].items()},
                      fieldstar_n=r1[s]['field_star_n']) for s in S},
    registration={s: dict(pair_info_fraction_along_u={k: v['frac_info_pair_along_u'] for k, v in r2[s]['info'].items()},
                          error_needed_px=r2[s]['registration_error_needed_px'],
                          wcs_only_fieldstar_offset_along_u_px=r2[s]['wcs_only_offset_along_u_px']) for s in S},
    jitter={s: dict(leverage_analytic=r5[s]['leverage_analytic'], leverage_numeric_PRF=r5[s]['leverage_numeric_PRF'],
                    s_u_needed_mpx=r5[s]['s_u_needed_mpx'], s_u_per_transit_mpx=[o['s_u_mpx'] for o in r5[s]['real']],
                    s_u_stack_mpx=r5[s]['real_s_u_stack_mpx'], spurious_s_par_stack=r5[s]['real_spurious_s_par_stack'],
                    fake_s_u_mad_mpx=r5[s]['fake_s_u_mad_mpx']) for s in S},
    timing={s: dict(time_keys=r4[s]['tesscut_time_keys'], pixel_dt0_min=1440 * r4[s]['pair_aperture_dt0_d'],
                    pixel_dt0_err_min=1440 * r4[s]['pair_aperture_dt0_err_d'], ephemeris_sigma_min=1440 * r4[s]['ephemeris_sigma_d'],
                    depth_peak_shift_d=r4[s]['depth_vs_shift']['shift_d'][int(np.argmax(r4[s]['depth_vs_shift']['depth_ppt']))],
                    pair_r2_depth_ppt=r4[s]['pair_aperture_depth_ppt'], flag_frac_in_transit=r4[s]['flag_frac_in_transit'],
                    flag_frac_near_oot=r4[s]['flag_frac_near_oot'],
                    tesscut_minus_qlp_nearest_time_s=r4[s]['tesscut_minus_qlp_nearest_s']['median'],
                    centroidpy_oot_frac_inside_transit=r4[s]['centroidpy_oot_frac_inside_true_transit']) for s in S},
    pixel_test=dict(
        definition='delta(p) = D(p)/O(p); Q = delta(x)/delta(y); H_target => Q = fT(x)/fT(y) <= 1/fT(y)',
        per_sector={s: {p: dict(x=v['x'], y=v['y'], dF_x=v['dF_x'], dF_y=v['dF_y'], delta_x_ppt=v['delta_x_ppt'],
                                delta_y_ppt=v['delta_y_ppt'], sigma_delta_x_ppt=v['sigma_delta_x_null_ppt'],
                                sigma_delta_y_ppt=v['sigma_delta_y_null_ppt'], Q=v['Q'],
                                Q_T_scene_range=v['model_context']['scene']['Q_T_range'],
                                Q_N_scene_range=v['model_context']['scene']['Q_N_range'],
                                z_T=v['z_T_models'], z_N_range=v['z_N_models'], kappa_model_free=v['kappa'],
                                z_model_free_kappa=v['z_L'])
                        for p, v in r3[s]['pairs'].items()} for s in S},
        combined={k: v for k, v in r6['pixel_test_none_plane'].items()},
        combined_cleaned={k: dict(stouffer=v['stouffer'], empirical=v.get('empirical_combined_null'))
                          for k, v in r6['pixel_test_mad10_plane'].items()},
        zero_model_flux_inequality=dict(no_cleaning=r3c, cleaned=r3cm)),
    statistics=dict(effective_windows=r6['effective_window_counts'], quoted=r6['quoted_significances_recomputed']),
    centroid_py=dict(
        replication={s: (dict(offset_arcsec=b1[s]['real']['offset_arcsec'], offset_px=b1[s]['real']['offset_px'],
                              err_px=b1[s]['real']['err_px'], offset_sigma=b1[s]['real']['offset_sigma'],
                              s_par=b1[s]['real']['s_par'], z_vs_nb_errpx_only=b1[s]['real']['z_vs_neighbour_errpx_only'],
                              z_vs_nb_with_floor=b1[s]['real']['z_vs_neighbour_with_floor'],
                              per_window_floor=b1[s]['real']['per_window_floor_e_s_px'],
                              estimator_variants_s_par={k: v['s_par'] for k, v in b1[s]['real']['estimator_variants'].items()},
                              noiseless_s_par={k: v['s_par'] for k, v in b1[s]['real']['noiseless_models'].items()})
                         if b1[s]['real'] else None) for s in S},
        injections={}),
)
for host in ('target', 'neighbour'):
    for m in ('prf', 'gauss0.6'):
        v = [o for s in S for o in b1[s]['injections'][f'{m}_{host}'] if o]
        a = lambda k: np.array([o[k] for o in v])
        res['centroid_py']['injections'][f'{m}_{host}'] = dict(
            n=len(v), offset_arcsec_median=float(np.median(a('offset_arcsec'))),
            offset_arcsec_16_84=np.percentile(a('offset_arcsec'), [16, 84]).tolist(),
            offset_sigma_median=float(np.median(a('offset_sigma'))), offset_sigma_max=float(a('offset_sigma').max()),
            frac_offset_sigma_ge_3=float(np.mean(a('offset_sigma') >= 3)),
            s_par_median=float(np.median(a('s_par'))), s_par_16_84=np.percentile(a('s_par'), [16, 84]).tolist(),
            frac_lead_statistic_gt_6=float(np.mean(a('z_vs_neighbour_errpx_only') > 6)))
dump(res, f'{OUT}/results.json')
print(json.dumps(res['centroid_py']['injections'], indent=1))
