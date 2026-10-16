"""Collect the key referee numbers into one JSON."""
import json, glob, os
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee'
S = {
 "qlp_timestamp_issue": {"finding": "QLP S103/S104 TIME lacks barycentric correction: QLP = TESScut TIME - TIMECORR within ~5 s for 5 targets (S104) and 1 (S103); momentum-dump flag at 4185.0384 (QLP) vs 4185.0418 (SPOC 2-min TIC 90087028) vs 4185.0433 (TESScut); QLP-TESScut eclipse times for TIC 90086266 = -7.2 +- 0.3 min (TIMECORR 7.5 min). S13-S94 QLP agree with TESScut to <=20 s.",
                         "scripts": ["code/t1_events.py"], "data": ["data/dl_spoc104/"]},
 "TIC90086266": {"verdict": "CONFIRM (not a planet); weaken precision of claims a and c",
   "cycle_count_N": 1078, "next_alias_sigma": 45,
   "ephem_barycentric": {"T0_BTJD": 4180.47842, "T0_err": 0.00019, "P": 2.3394739, "P_err": 2e-7, "chi2": "25.8/18"},
   "ephem_with_QLP_S104": {"T0_BTJD": 4180.47342, "P": 2.3394693, "chi2": "38.4/17"},
   "oddeven_S13_pct": "+1.0..+1.7 (0.3-1.2 sigma; TESS-SPOC SAP -0.2)", "oddeven_S104_pct": "+5.4..+6.2 (4.3-6.3 sigma)",
   "oddeven_S104_pixel_ratio": "0.910, chi2 10.6/11 (uniform)", "oddeven_instrumental_correlates": "none (bkg, centroid, pointing, flags, flux level)",
   "secondary_phase05_ppm": {"S13_PDCSAP": [358, 165], "S13_SAP": [99, 157], "S13_CUT": [-66, 223], "S104_QLP": [-453, 126], "S104_CUT": [-286, 190]},
   "expected_secondary_Mdwarf_ppm": "1600-2200 (BB 3200-3400 K)", "R2_min_Rsun": 0.26,
   "K_0p27Msun_kms": 31.6, "P_range_le_5kms": {"random5": 3.6e-4, "clustered_3visits": 5.0e-3, "clustered_2visits": 1.0e-2},
   "max_depth_gaia_nbrs_pct": {"G15.83@3.4arcsec": 1.94, "G15.72@6.0arcsec": 2.14}, "observed_raw_depth_pct": "2.7-3.5",
   "diffimg_offset_S13_arcsec": "(+0.2..+0.8, +0.3..+0.8) +-0.5; injection (0,0)", "diffimg_offset_S104_arcsec": "(+0.8,+1.2) [up to (-0.8,+3.1)]; injections (-0.9,-0.7)+-(0.6,0.3) or (+1.4,+0.5) with planar bkg"},
 "TIC61332759": {"verdict": "CONFIRM",
   "harmonic_map_offsets_from_target": {"S13": [-36.0, -43.9], "S27": [-36.0, -42.1], "S67": [-35.2, -43.3]}, "gaia_neighbour_offset": [-36.2, -44.0],
   "dip_phased_diffimg_offsets": {"S27": [-36.0, -42.0], "S67": [-36.7, -43.3], "S94": [-36.8, -43.2], "S104": [-34.4, -42.9]},
   "neighbour_gaia_class": "RS (0.61)", "target_r1_decorrelated_dip_ppm": {"S13": [277, 113], "S27": [55, 104], "S94": [171, 111], "S104": [-182, 139], "S67": [-883, 164]}},
 "TIC214311906": {"verdict": "CONFIRM not a planet of the Mira; OVERTURN 'blended EB within ~10 arcsec, unidentified'",
   "orbit_a_Rsun": "13-19 for 1-3 Msun", "host_density_gcc_b0": 10.7, "per_sector_snr_at_ephem": [4.2, 6.3, 5.9, 7.2, 6.0, 7.9], "combined_snr": 15.5, "trial_phase_max": 6.4,
   "source_offset_arcsec": [-17.2, 3.7], "source_offset_err": [1.7, 1.0],
   "injection_chi2_2dof": {"TIC214311908_A": 1.6, "TIC214311910_B": 28.9, "Mira": 56.0},
   "abs_dip_eps_r0.75": {"A": [5.13, 0.41], "B": [4.35, 0.48], "Mira": [2.55, 0.41], "mirror": [1.28, 0.36]},
   "QLP_depth_ppm_TIC214311908": [12741, 15114, 19712, 15613, 10658, 14630], "QLP_depth_ppm_Mira": [1431, 2844, 3606, 3076, 1227, 2412],
   "deltaF_eps": [16.5, 1.1], "host": "TIC 214311908 = Gaia DR3 5925373272850784768 (M3V, 64 pc, T=13.01, R=0.381, M=0.369; CPM with Gaia DR3 5925373272850782208)",
   "depth_on_host_pct": [1.8, 0.2], "Rp_Rearth": "5-6", "stack_T14_h": 1.65, "stack_T12_over_T14": 0.20, "V_shape_dchi2": 8.0, "oddeven_ppm": [[1983, 273], [1613, 194]],
   "my_tesscut_depth_on_Mira_ppm": "550-1340 vs QLP 1070-3740"},
 "TIC55829034": {"verdict": "CONFIRM", "per_sector_blind_BLS_P": [0.6162, 0.6153, 0.6162, 0.6163], "host_density_gcc_b0": 0.65, "a_Rsun_max": 3.8,
   "source_offset_arcsec": [3.2, -5.5], "err": [0.8, 0.9], "sector_scatter": 3, "catalog_match": "none (nearest Gaia EBs 88-164 arcsec; 2x0.307166=0.61433 != 0.616146)"},
 "TIC61141412": {"verdict": "CONFIRM", "lin_chi2_range": "90-240 (10-13 dof)", "quad_chi2_range": "12-31 (9-12 dof)",
   "OC_min_by_sector": {"S13": "-10..-19", "S27": "-4..+7", "S67": "+17..+27", "S94": "-11..0", "S104": "-4..-12"},
   "S104_QLP_time_shift_min": 6.45, "P": 6.52854, "half_phase_ppm": [53, 72], "odd_even_ppm": [[1371, 106], [1548, 102]],
   "gaia_sb1": {"P": 5.452382, "P_err": 0.00048, "K1": 18.13, "e": 0.011, "n_rv": 10}, "spot_period_d": "5.4-5.8", "vbroad": 11.0}
}
json.dump(S, open(f'{B}/referee_summary.json', 'w'), indent=1)
print('written', f'{B}/referee_summary.json')
