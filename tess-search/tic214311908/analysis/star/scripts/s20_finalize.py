"""Adds the system, kinematic, activity, companion-limit, neighbour and planet summaries to star_params.json
(written by s11) from the outputs of s09, s10, s12b/c, s13b/c, s14, s15, s16, s17, s18. Run last."""
import json, numpy as np
B = 'tic214311908/analysis/star'
P = json.load(open(f'{B}/star_params.json'))
L = lambda f: json.load(open(f'{B}/out/{f}'))
kin, act, met, cmd, ipd, ast, dss, dec, bl, pl, misc = (L('kinematics.json'), L('activity.json'), L('metallicity.json'), L('cmd.json'), L('ipd_limits.json'),
                                                       L('astrometric_limits.json'), L('dss_limits.json'), L('decaps_contrast.json'), L('blend_density.json'),
                                                       L('planet_params.json'), L('gaia_misc_calib.json'))
d = P['distance_pc']['median']
# Mamajek mass from M_G
tab = []
for line in open(f'{B}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] in 'KM' and p[0].endswith('V'):
        try: tab.append((float(p[13]), float(p[30])))
        except Exception: pass
tab = np.array(sorted(tab))
mass_of_MG = lambda mg: float(np.interp(mg, tab[:, 0], tab[:, 1]))
MG_A = 14.654037 + 5 * np.log10(15.5494 / 100); MG_B = 15.740896 + 5 * np.log10(15.5494 / 100)
mA, mB = mass_of_MG(MG_A), mass_of_MG(MG_B)
dpm = np.array([-112.795992 + 113.671152, 87.205288 - 85.54355]); edpm = np.hypot([0.441198, 0.322918], [0.031124, 0.0224])
mtot = P['mass_msun']['median'] + mA + mB
P['wide_system'] = dict(
    members=dict(M3V='TIC 214311908 = Gaia DR3 5925373272850784768', A='TIC 214311910 = Gaia DR3 5925373272850782208 (G 14.654, RUWE 13.7)',
                 B='Gaia DR3 5925373272827253632 (G 15.741, 2-parameter solution, 0.83" from A)'),
    sep_M3V_A_arcsec=11.728, sep_M3V_A_au=float(11.728 * d), sep_A_B_arcsec=0.83, sep_A_B_au=float(0.83 * d),
    parallax_A=[16.4126, 0.3929], parallax_diff_sigma_formal=float((16.4126 - 15.4973) / np.hypot(0.3929, 0.0251)),
    parallax_note='A has RUWE 13.7 (unresolved B in its solution), so its formal parallax error is underestimated; the 2.3-sigma difference is not significant',
    dpm_masyr=dpm.tolist(), dpm_err=edpm.tolist(), dv_tan_kms=float(4.74 * np.hypot(*dpm) / 15.5494),
    v_escape_at_projected_sep_kms=float(29.78 * np.sqrt(2 * mtot / (11.728 * d))),
    MG_A=float(MG_A), MG_B=float(MG_B), mass_A_msun_Mamajek=mA, mass_B_msun_Mamajek=mB,
    elbadry21=dict(in_catalogue=True, sepAU=756.8, R_chance_alignment=2e-5, bin_type='MSMS'),
    B_chance_alignment_within_0p83=float(0.00056 * np.pi * 0.83 ** 2),
    other_comovers_within_2deg='none (s17: 94 Gaia sources with parallax 13-18.5 mas; only A has dv_tan < 3 km/s)',
    hierarchy='(M3V) -- ~760 au -- (A -- ~53 au -- B); a_out/a_in ~ 14, dynamically stable; no Gaia NSS solution, no RV for A/B; '
              'M3V RV 16.28 +- 2.22 km/s (21 transits). If the M3V is itself an unresolved near-equal binary the system is a quadruple.',
    outer_orbit_period_yr_approx=float((11.728 * d) ** 1.5 / np.sqrt(mtot)), inner_AB_period_yr_approx=float((0.83 * d) ** 1.5 / np.sqrt(mA + mB)))
P['kinematics'] = dict(UVW_kms=[kin['U'][0], kin['V'][0], kin['W'][0]], UVW_err=[kin['U'][1], kin['V'][1], kin['W'][1]], v_tot_LSR=kin['vtot_lsr'],
                       W_LSR=kin['W_lsr'], P_thin_Bensby14=kin['P_thin'], P_thick_Bensby14=kin['P_thick'], TD_over_D=kin['TD_over_D'],
                       banyan_sigma='field (P_field = 1.00 with and without RV)',
                       interpretation='W_LSR = +48 km/s: hot kinematics typical of an old (several Gyr) thin-disk or thick-disk star; not a member of any young association')
xr = act['xray']
P['activity'] = dict(xray=dict(source=xr['name'], det_like=xr['det_like'], counts=xr['counts'], flux_0p2_2p3_cgs=xr['flux_0p2_2p3keV'], log_LX=xr['log_LX'],
                               LX_over_Lbol_if_all_from_M3V=xr['LX_over_Lbol_if_M3V'], association_rel_likelihoods={k: v['rel_likelihood'] for k, v in xr['association'].items()},
                               Ro_and_Prot_implied=dict(all_from_M3V=[xr['Ro_all_from_M3V'], xr['Prot_days_all_from_M3V']], half=[xr['Ro_half_from_M3V'], xr['Prot_days_half_from_M3V']]),
                               chance_coincidence_10arcsec=xr['chance_coincidence_within_10arcsec']),
                     galex='no coverage (0 GUVcat/AIS sources within 15\')', rosat_2rxs='none within 20"', xmm_swift='none within 20"',
                     rotation='not measured: deblended TESS photometry detects sinusoids only above 1.5-5% semi-amplitude (90%); highest periodogram peaks (3.6-10 d) vary between sectors and coincide with peaks in the pair light curve',
                     flares='no convincing flare in 111 d (5-sigma per-cadence thresholds 40-90 ppt)',
                     gaia_halpha_espels_ew_nm=[0.0971, 0.0404])
P['photometric_metallicity'] = dict(dMKs_from_ridge=cmd['MKs_vs_BPRP']['target_offset_mag'], dMG_from_ridge=cmd['MG_vs_BPRP']['target_offset_mag'],
                                    frac_comparison_as_overluminous=cmd['MKs_vs_BPRP']['frac_comparison_as_bright_or_brighter'],
                                    feh_if_single=met['feh_phot_if_single'], calibration=dict(n=met['n_cal'], fit=met['fit_feh_vs_dMK'], scatter=met['scatter_feh']),
                                    pairA_deblended=met['pairA_deblended'])
P['companion_limits'] = dict(
    gaia_ipd=dict(clean_single_rate=ipd['single']['L'], red_variant=ipd['red_variant'], summary='companions with Delta(I) < 2 at 0.08-2.5" (5-160 au) excluded (likelihood ratio < 0.1); Delta 2-3.5 disfavoured at 0.08-0.8", excluded at 0.8-2.5"; Delta > 3.5 essentially unconstrained'),
    astrometry=dict(file='out/astrometric_limits.json', summary='RUWE 1.076 + DR3-DR2 PM: bound stellar companions 0.08-0.3 Msun excluded at ~0.1-4 au, brown dwarfs >= 40 MJup at ~0.24-2.4 au; near-equal masses (>= 0.34 Msun, photocentre ~ barycentre) not excluded at any separation'),
    decaps2_residual_5sigma=dict((f'{k}', v) for k, v in json.load(open(f'{B}/out/exclusion_summary.json'))['decaps_dT_5sig'].items()),
    dss_lookback_unbound={k: dict(within2=v['grid_T_lim_min_within2'], within3=v['grid_T_lim_min'], calib_sd_dex=v['calib']['robust_sd_dex']) for k, v in dss.items()},
    gaia_misc=misc)
P['blend_density'] = bl
P['planet_summary'] = dict(file='out/planet_params.json', single=dict(Rp_earth=pl['single']['Rp_b06_depth1p64'], a_au=pl['single']['a_au'], S=pl['single']['insolation_Searth'],
                                                                     Teq=pl['single']['Teq_A0'], TSM=pl['single']['TSM'], ESM=pl['single']['ESM'], K_ChenKipping=pl['single']['rv']['ChenKipping17']),
                           twin=dict(Rp_earth=pl['twin']['Rp_b06_depth1p64'], a_au=pl['twin']['a_au'], Teq=pl['twin']['Teq_A0']))
json.dump(P, open(f'{B}/star_params.json', 'w'), indent=1, default=float)
print(json.dumps(P['wide_system'], indent=1))
