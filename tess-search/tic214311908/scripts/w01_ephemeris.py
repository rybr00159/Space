"""w01: write writeup/ephemeris_adopted.json.
Adopted = the conservative planning ephemeris suggested by referee3_clean (README item 4):
T0 = 2460840.9150 +- 0.0013 BJD_TDB (epoch 0), P = 5.691506 +- 0.000011 d.
The referee did not state a T0-P correlation for these rounded values; we adopt rho = 0.35, the correlation at epoch 0 of the
referee's correlated combined fit (referee3_clean/out/q05_timing.json comb_cor 0.345; comb_ind 0.377), which reproduces the
referee's '1 sigma ~3 min in 2027' (2.9-3.5 min). The file uses the block layout of clean/out/ephemeris_final.json
('adopted' names the block), so followup/scripts/f01_transits.py and referee3_followup/scripts/r01_times.py can read it."""
import json, os, numpy as np
W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the write-up folder (writeup/ or tic214311908/)
B = os.path.dirname(W) if os.path.isdir(os.path.join(os.path.dirname(W), 'clean')) else os.path.join(W, 'analysis')
cl = json.load(open(f'{B}/clean/out/ephemeris_final.json'))
q5 = json.load(open(f'{B}/referee3_clean/out/q05_timing.json'))
T0, sT0, P, sP, rho = 2460840.9150, 0.0013, 5.691506, 0.000011, 0.35
cov = [[sT0**2, rho*sT0*sP], [rho*sT0*sP, sP**2]]
def sig_min(E):
    return float(np.sqrt(cov[0][0] + E*E*cov[1][1] + 2*E*cov[0][1])*1440)
def blk_from_clean(name):
    b = cl[name]
    return dict(T0_epoch0_BJD_TDB=b['T0_epoch0_BJD_TDB'], T0_epoch0_err_d=b['T0_epoch0_err_d'], P_d=b['P_d'],
                P_err_d=b['P_err_d'], reference_epoch=b['reference_epoch'], T0_ref_BJD_TDB=b['T0_ref_BJD_TDB'],
                T0_ref_err_d=b['T0_ref_err_d'], cov_T0ref_P=b['cov_T0ref_P'], chi2=b['chi2'], dof=b['dof'])
def blk_from_q5(name):
    b = q5[name]
    s0, sp, r = b['T0_err'], b['P_err'], b['corr']
    return dict(T0_epoch0_BJD_TDB=b['T0_E0'] + 2457000, T0_epoch0_err_d=s0, P_d=b['P'], P_err_d=sp, corr_epoch0=r,
                reference_epoch=0, T0_ref_BJD_TDB=b['T0_E0'] + 2457000, cov_T0ref_P=[[s0*s0, r*s0*sp], [r*s0*sp, sp*sp]],
                chi2=b['chi2'], dof=b['dof'], sigma_2027_min_E100_130_160=[p['sig_min'] for p in b['predictions']])
out = {
 "adopted": "adopted_planning",
 "adopted_planning": dict(
    description="Conservative planning ephemeris (referee3_clean README item 4). Use +-10 min (about 3 sigma) margins around predicted mid-times in 2027-2028.",
    reference_epoch=0, T0_ref_BJD_TDB=T0, T0_ref_err_d=sT0, T0_ref_BTJD_TDB=round(T0 - 2457000, 4),
    P_d=P, P_err_d=sP, corr_T0ref_P=rho, cov_T0ref_P=cov,
    corr_note="rho = 0.35 assumed: referee3_clean comb_cor has 0.345 at epoch 0 (comb_ind 0.377); with rho = 0 the 2027 1-sigma would be 2.6-3.0 min instead of 2.9-3.5 min",
    sigma_min_at_epochs={str(E): round(sig_min(E), 2) for E in (100, 108, 130, 149, 160, 227)},
    time_system="BJD_TDB, target-barycentred (TESScut TIME - TIMECORR + Roemer delay to the M3V; phot p00 / referee_phot r02 / referee3_clean q01 agree to < 0.6 s)",
    epoch_numbering="epoch 0 = the transit at BTJD 3840.915 (2025-06-14, in Sector 93); phase-1 numbering used in every analysis folder",
    T14_h=1.70, planning_window="T14 + 2 x 3 sigma + 2 x 1.5 h baseline",
    source="referee3_clean/README.md (item 4), referee3_clean/out/q05_timing.json, referee3_clean/out/q05b_eph_variants.json"),
 "why_this_one": [
    "clean/ adopted the combined INV+REF fit T0 = 2460840.91488 +- 0.00093, P = 5.6915020 +- 0.0000078 d (clean/out/ephemeris_final.json).",
    "referee3_clean re-measured the per-transit times (reproduced to 0.1 min rms) but found that P moves by about 1 sigma with the timing-error model: with sector-uniform light-curve-level injection-calibrated errors the combined fit gives P = 5.6915097-5.6915111 and T0 = 3840.91500-3840.91510 BTJD_TDB (q05_timing.json comb_ind / comb_cor); global MCMC fits give P = 5.6915045-5.6915067.",
    "The 2027 predictions of the clean and referee fits differ by +1.5 to +2.2 min, comparable to their quoted 1 sigma (clean 1.9-2.5 min; referee 2.6-3.5 min). referee3_followup adds that the INV, REF and phase-1 ephemerides differ from the clean one by up to 2.2 min over 2026-2028 (pipeline systematic of about 1 min).",
    "The adopted values sit between the clean and referee solutions with errors that cover both; 1 sigma is about 3 min in 2027, so +-10 min windows are about 3 sigma.",
    "Difference adopted minus clean combined: +0.17 min at epoch 0, +0.9 min at epoch 130 (mid-2027), +1.5 min at epoch 227 (end of 2028)."],
 "alternatives": {
    "clean_combined": blk_from_clean('combined'),
    "clean_INV": blk_from_clean('inv'),
    "clean_REF": blk_from_clean('ref'),
    "referee3_comb_cor": blk_from_q5('comb_cor'),
    "referee3_comb_ind": blk_from_q5('comb_ind'),
    "referee3_INV": blk_from_q5('inv'),
    "referee3_REF": blk_from_q5('ref'),
    "phase1_phot_global": dict(T0_epoch0_BJD_TDB=2460840.91414, T0_epoch0_err_d=0.00134, P_d=5.6915087, P_err_d=0.0000090, source="phot/out/transits.json (phot/README.md)"),
    "phase1_phot_timing": dict(T0_epoch0_BJD_TDB=2460840.91272, T0_epoch0_err_d=0.00175, P_d=5.6915092, P_err_d=0.0000122, source="phot/out/timing.json"),
    "phase1_referee_phot": dict(T0_epoch0_BJD_TDB=2460840.91564, T0_epoch0_err_d=0.00118, P_d=5.6915051, P_err_d=0.0000092, source="referee_phot/out/oddeven_timing.json (README)"),
    "phase2_shape": dict(T0_epoch0_BJD_TDB=2460840.9155, T0_epoch0_err_d=0.0015, P_d=5.691508, P_err_d=0.000009, source="shape/out/final_params.json"),
    "repo_vetting_v2_SUPERSEDED": dict(T0_epoch0_BJD_TDB=2460840.9153, T0_epoch0_err_d=0.0019, P_d=5.6914968, P_err_d=0.0000145,
        source="tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json (13 QLP transits, errors scaled by chi2/dof 53.4/11)",
        note="used for the repo's predicted_transits_2027.csv; its P is 0.8 sigma (adopted errors) shorter than adopted; its mid-2027 (epoch 130) mid-time is 1.3 min earlier than adopted"),
    "repo_findings_first_SUPERSEDED": dict(T0_epoch0_BJD_TDB=2460840.9101, T0_epoch0_err_d=0.0022, P_d=5.691471, P_err_d=0.000018,
        source="tess-search/FINDINGS.md line 440 (14 events, uncorrected QLP S103/S104 time stamps)")}
}
# difference adopted - clean combined
c = cl['combined']
for E in (0, 130, 227):
    tc = c['T0_ref_BJD_TDB'] + c['P_d']*(E - c['reference_epoch'])
    ta = T0 + P*E
    print(E, round((ta - tc)*1440, 2), 'min; sigma', round(sig_min(E), 2))
json.dump(out, open(f'{W}/ephemeris_adopted.json', 'w'), indent=1)
print('written')
