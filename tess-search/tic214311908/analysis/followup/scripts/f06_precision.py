"""f06: ground-based photometric precision, transit S/N, saturation and the depth-ratio (achromaticity) test.
Noise model per exposure: photon noise of the star, sky (dark-sky surface brightness, ASSUMED values below), read noise,
and scintillation from the modified Young formula of Osborn et al. 2015 (MNRAS 452, 1707, eq. 7):
  sigma_Y^2 = 10e-6 C_Y^2 D^(-4/3) t^(-1) (cos z)^(-3) exp(-2 h/H), H = 8000 m,
with C_Y = 1.42 (Tololo median; used for CTIO, Pachon, La Silla, Las Campanas), 1.56 (Paranal median), 1.5 elsewhere
(Siding Spring, Sutherland, amateur: not in Osborn's table; ASSUMED) and 1.0 at Dome C (ASSUMED, weak high-altitude
turbulence); quartile variant with C_Y x 1.3. Scintillation of the comparison ensemble: factor 1.2. Airmass 1.3.
Count rates: AB magnitudes of the M3V (XP-synthetic / model-offset values from f05) and the photon zero point
5.48e10 (dlambda/lambda) photons s^-1 m^-2 at AB = 0, times collecting area x total throughput (ASSUMED, +-30%).
Per-event depth error: white part over T_in = 1.45 h in transit and 3 h out of transit (2-min bins), times 1.42 for a
quadratic baseline fitted per event (calibrated by the f10 simulation), plus a red-noise term sigma_red per event
(ASSUMED: 0.5 ppt for 1-2 m, 0.8 ppt for 0.4-0.6 m; variant 2x; f10 shows an OU process with tau 30 min and rms s
adds ~0.9-1.0 s to the per-event depth error).
Depth-ratio test: sigma(R)/R = sqrt((s1/d1)^2 + (s2/d2)^2); N(5%) = (sigma_R/0.05)^2. Hypothesis predictions from
out/chromatic.json (planet: ratio 1 after LD; PTP; WD; background K-dwarf EB).
Saturation: peak-pixel rate of a Moffat (beta 2.5) PSF at FWHM 1.0" (worst case) for the M3V and the Mira at maximum.
Outputs: out/precision.csv, out/depth_ratio_tests.csv, out/saturation.csv, out/precision_summary.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json
con = json.load(open(f'{OUT}/contamination_summary.json'))
chrom = json.load(open(f'{OUT}/chromatic.json'))
m3v = con['m3v_mags']; mira = con['mira_mag_ranges']
TOAB = {'V': 0.02, 'Ic': 0.45, 'J': 0.91}
def ab(b, m): return m + TOAB.get(b, 0.0)
BAND = {  # centre nm, width nm, dark sky AB mag/arcsec^2 (ASSUMED), extinction mag/airmass (ASSUMED)
    'g': (477, 137, 22.0, 0.20), 'r': (623, 137, 21.2, 0.12), 'i': (763, 153, 20.4, 0.08), 'z': (870, 100, 19.7, 0.06),
    'V': (551, 88, 21.7, 0.15), 'Ic': (806, 149, 20.1, 0.06), 'I+z': (900, 250, 19.8, 0.06), 'ASTEP_blue': (555, 293, 21.7, 0.15),
    'ASTEP_red': (850, 276, 20.0, 0.06), 'ExTrA': (1200, 700, 15.8, 0.08), 'J': (1235, 162, 16.9, 0.08)}
TEL = {  # name: D [m], area [m^2], pix ["], RN [e-], full well [e-], readout [s], C_Y, h [m], bands {band: throughput}, sigma_red [ppt], fwhm_typ ["], notes
 'LCOGT 1m Sinistro': dict(D=1.0, A=0.60, pix=0.389, RN=8.0, FW=1.0e5, ro=28.0, CY=1.42, h=2207, thr={'g': 0.35, 'r': 0.42, 'i': 0.38, 'z': 0.22, 'V': 0.38, 'Ic': 0.35}, red=0.5, fwhm=2.0,
                           src='LCO Sinistro page: 0.389"/px, RN 7.5-8 e-, FW 100-120 ke-, 28 s cycle (9 s for 2x2 central window)'),
 'FTS 2m MuSCAT4': dict(D=2.0, A=2.8, pix=0.27, RN=12.0, FW=7.1e4, ro=6.0, CY=1.5, h=1165, thr={'g': 0.35, 'r': 0.40, 'i': 0.35, 'z': 0.25}, red=0.5, fwhm=1.8,
                        src='LCO MuSCAT page: simultaneous g r i zs, 0.27"/px, 9.1\' FOV, RN 12 e- (17 zs) fast, FW >= 71 ke-, 6 s readout; installed Oct 2023 at Siding Spring'),
 'SPECULOOS-South 1m': dict(D=1.0, A=0.724, pix=0.35, RN=6.0, FW=1.0e5, ro=10.0, CY=1.56, h=2635, thr={'I+z': 0.45}, red=0.5, fwhm=1.2,
                            src='Paranal; Andor iKon-L deep-depletion 2k, 0.35"/px, 12\' FOV, I+z (>90% from 750 nm); RN/FW/readout ASSUMED'),
 'TRAPPIST-South 0.6m': dict(D=0.6, A=0.24, pix=0.64, RN=9.0, FW=1.0e5, ro=6.0, CY=1.42, h=2400, thr={'I+z': 0.35, 'z': 0.25, 'V': 0.30}, red=0.8, fwhm=1.5,
                             src='La Silla; FLI ProLine 2k, 0.64"/px, 22\' FOV; RN/FW/readout ASSUMED'),
 'ExTrA 0.6m (8" fibre)': dict(D=0.6, A=0.24, pix=None, fib=8.0, RN=20.0, FW=None, ro=2.0, CY=1.42, h=2400, thr={'ExTrA': 0.10}, red=0.5, fwhm=1.5,
                               src='La Silla; 8" (or 4") fibres, R~20, 0.85-1.55 um, 60-s exposures (ExTrA TOI papers); throughput ASSUMED'),
 'ExTrA 0.6m (4" fibre)': dict(D=0.6, A=0.24, pix=None, fib=4.0, RN=20.0, FW=None, ro=2.0, CY=1.42, h=2400, thr={'ExTrA': 0.09}, red=0.5, fwhm=1.5, src='as above, 4" fibres'),
 'ASTEP+ 0.4m (Dome C)': dict(D=0.4, A=0.10, pix=1.39, RN=5.0, FW=6.0e4, ro=2.0, CY=1.0, h=3233, thr={'ASTEP_blue': 0.35, 'ASTEP_red': 0.20}, red=0.8, fwhm=2.0,
                              src='Dome C; two channels 555/293 and 850/276 nm, 1.39"/px, ~40\' FOV (ASTEP+ papers); RN/FW/throughput ASSUMED (red channel transmits less)'),
 'Amateur 0.4m': dict(D=0.4, A=0.10, pix=0.6, RN=3.0, FW=5.0e4, ro=2.0, CY=1.5, h=500, thr={'Ic': 0.30, 'r': 0.35}, red=1.0, fwhm=2.5, sky_off=-1.0,
                      src='generic CMOS 0.4 m; all ASSUMED'),
}
AIRMASS = 1.3
BASEFIT = 1.42   # f10 Monte Carlo: depth sd / analytic white term = 1.40-1.48 with a quadratic baseline fitted per event
T_IN_H, T_OUT_H = 1.45, 3.0
def rates(tel, b):
    c, w, sky, k = BAND[b]
    thr = TEL[tel]['thr'][b]*10**(-0.4*k*AIRMASS)
    n0 = 5.48e10*(w/c)*TEL[tel]['A']*thr                     # e-/s at AB = 0
    star = n0*10**(-0.4*ab(b, m3v[b]))
    skyab = sky + TEL[tel].get('sky_off', 0.0)
    return n0, star, n0*10**(-0.4*skyab)                    # sky per arcsec^2
def scint(tel, t):
    T = TEL[tel]; return np.sqrt(10e-6*T['CY']**2*T['D']**(-4/3)/t*AIRMASS**3*np.exp(-2*T['h']/8000.))*1.2
rows, sat = [], []
for tel, T in TEL.items():
    for b in T['thr']:
        n0, star, skyas = rates(tel, b)
        fwhm = T['fwhm']
        if T.get('fib'):
            rap = T['fib']/2; area = np.pi*rap**2; enc = moffat_enclosed(rap, fwhm, 2.5); npix = 40.0   # spectrum pixels (ASSUMED)
            texp = 60.0
            mira_peak = None; m3v_peak = None; tsat_m3v = np.inf; tsat_mira = np.inf
        else:
            rap = max(1.5, 1.0*fwhm); area = np.pi*rap**2; enc = moffat_enclosed(rap, fwhm, 2.5); npix = area/T['pix']**2
            pk = moffat_peak_frac(1.0, 2.5, T['pix'])
            m3v_peak = star*pk; mb = mira.get(b if b in mira else {'I+z': 'z', 'ASTEP_blue': 'V', 'ASTEP_red': 'z'}.get(b, b))[0]
            mira_peak = n0*10**(-0.4*ab(b, mb))*pk
            tsat_m3v = 0.7*T['FW']/(m3v_peak + skyas*T['pix']**2); tsat_mira = 0.7*T['FW']/mira_peak
            texp = float(min(120.0, max(5.0, 0.8*tsat_m3v)))        # in focus at FWHM 1.0" worst case; else defocus
        ncyc = 120.0/(texp + T['ro'])
        S = star*enc*texp; Bk = skyas*area*texp; R2 = npix*T['RN']**2
        sig_ph = np.sqrt(S + Bk + R2)/S
        sig_sc = scint(tel, texp)
        sig_exp = np.hypot(sig_ph, sig_sc)
        sig2 = sig_exp/np.sqrt(max(ncyc, 1.0)) if ncyc >= 1 else sig_exp
        nin = T_IN_H*30; nout = T_OUT_H*30
        white = BASEFIT*sig2*np.sqrt(1/nin + 1/nout)   # x1.42: quadratic-baseline penalty calibrated by simulation (f10)
        sdel = np.hypot(white, T['red']*1e-3)
        sdel2 = np.hypot(white, 2*T['red']*1e-3)
        sdel_q = np.hypot(BASEFIT*np.hypot(sig_ph, sig_sc*1.3)/np.sqrt(max(ncyc, 1.0))*np.sqrt(1/nin + 1/nout), T['red']*1e-3)
        rows.append(dict(telescope=tel, band=b, m3v_AB=round(ab(b, m3v[b]), 2), star_e_per_s=round(star*enc), sky_e_per_s_ap=round(skyas*area),
                         r_ap=rap, fwhm=fwhm, t_exp=round(texp, 1), n_exp_per_2min=round(ncyc, 2), sigma_phot_exp_ppt=round(sig_ph*1e3, 3),
                         sigma_scint_exp_ppt=round(sig_sc*1e3, 3), sigma_2min_ppt=round(sig2*1e3, 3), sigma_depth_white_ppt=round(white*1e3, 3),
                         sigma_red_ppt=T['red'], sigma_depth_event_ppt=round(sdel*1e3, 3), SN_event=round(DEPTH/sdel, 1),
                         SN_event_red2x=round(DEPTH/sdel2, 1), SN_event_CYq3=round(DEPTH/sdel_q, 1)))
        if m3v_peak is not None:
            sat.append(dict(telescope=tel, band=b, M3V_peak_e_s_fwhm1=round(m3v_peak), t_sat_M3V_s=round(tsat_m3v, 1),
                            Mira_mag_bright=mb, Mira_peak_e_s_fwhm1=round(mira_peak), t_sat_Mira_s=round(tsat_mira, 2),
                            Mira_saturates_in_texp=bool(tsat_mira < texp), t_exp=round(texp, 1)))
P = pd.DataFrame(rows); P.to_csv(f'{OUT}/precision.csv', index=False)
S = pd.DataFrame(sat); S.to_csv(f'{OUT}/saturation.csv', index=False)
pd.set_option('display.width', 250)
print(P[['telescope', 'band', 'm3v_AB', 't_exp', 'sigma_2min_ppt', 'sigma_scint_exp_ppt', 'sigma_depth_event_ppt', 'SN_event', 'SN_event_red2x']].to_string())
print(S.to_string())
# ---- depth-ratio tests
H2 = chrom['H2_PTP_depth_over_TESS']; H3 = chrom['H3_WD']['depth_over_TESS']; H4 = chrom['H4_bg_Kdwarf_EB']
mapb = {'z': 'z_s'}
def Rpred(H, b1, b2):
    g = lambda b: H[mapb.get(b, b)]
    return g(b1)/g(b2)
PAIRS = {'LCOGT 1m Sinistro': [('g', 'i'), ('g', 'z'), ('V', 'Ic'), ('r', 'z')], 'FTS 2m MuSCAT4': [('g', 'z'), ('g', 'i'), ('r', 'z'), ('r', 'i')],
         'ASTEP+ 0.4m (Dome C)': [('ASTEP_blue', 'ASTEP_red')], 'Amateur 0.4m': [('r', 'Ic')], 'TRAPPIST-South 0.6m': [('V', 'I+z')]}
tests = []
for tel, prs in PAIRS.items():
    for b1, b2 in prs:
        s1 = P[(P.telescope == tel) & (P.band == b1)].iloc[0]; s2 = P[(P.telescope == tel) & (P.band == b2)].iloc[0]
        d1 = DEPTH*chrom['H1_planet_raw_depth_over_TESS'][mapb.get(b1, b1)]; d2 = DEPTH*chrom['H1_planet_raw_depth_over_TESS'][mapb.get(b2, b2)]
        sr = np.hypot(s1.sigma_depth_event_ppt*1e-3/d1, s2.sigma_depth_event_ppt*1e-3/d2)
        sim = tel in ('FTS 2m MuSCAT4', 'ASTEP+ 0.4m (Dome C)')
        row = dict(telescope=tel, pair=f'{b1}/{b2}', simultaneous=sim, sigma_ratio_per_event=round(sr, 3),
                   events_for_5pct=round((sr/0.05)**2*(1 if sim else 2), 1),
                   note='non-simultaneous bands need one event per band (x2 events)' if not sim else 'simultaneous')
        for name, H in (('PTP_dT1.16', H2['dT1.16']), ('PTP_dT0.88', H2['dT0.88']), ('PTP_dT1.45', H2['dT1.45']), ('WD', H3), ('bgK_EBV0.2', H4['E(B-V)=0.2'])):
            R = Rpred(H, b1, b2)
            row[f'R_{name}'] = round(R, 3)
            row[f'events_5sigma_{name}'] = round((5*sr/abs(R - 1))**2*(1 if sim else 2), 1) if abs(R - 1) > 1e-3 else np.inf
        tests.append(row)
Tt = pd.DataFrame(tests); Tt.to_csv(f'{OUT}/depth_ratio_tests.csv', index=False)
print(Tt.to_string())
jdump(dict(assumptions=dict(airmass=AIRMASS, T_in_h=T_IN_H, T_out_h=T_OUT_H, bands=BAND, telescopes={k: {kk: vv for kk, vv in v.items()} for k, v in TEL.items()},
                            note='throughputs, sky brightness, extinction, read noise of SPECULOOS/TRAPPIST/ASTEP/amateur and red-noise levels are ASSUMED; '
                                 'scale S/N by sqrt(throughput ratio) when better values are known'),
           best_SN=P.sort_values('SN_event', ascending=False)[['telescope', 'band', 'SN_event']].head(8).to_dict('records')),
      f'{OUT}/precision_summary.json')
