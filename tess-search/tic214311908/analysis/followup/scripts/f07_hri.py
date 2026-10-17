"""f07: high-resolution imaging design and the FPP after each follow-up step.
* Posterior-weighted bound-companion draws: ../referee2_fpp/out/r04_companions.npz with the PTP (planet on a fainter bound
  companion), BD-on-companion and TP-diluted weights of ../referee2_fpp/scripts/r05_fpp.py (imported read-only, no
  bytecode written), in two versions: as published, and with the photo-eccentric sign fix of ../bd (b10: rho_circ =
  rho_host g^+3, implemented exactly as b10 by returning 1/g from the eccentricity sampler).
* Contrast curves (5 sigma, this target: V 15.6, I 12.9, Ks 10.5): see CURVES; verified anchors:
  - SOAR/HRCam I: six-point median and 16th-percentile limits of 2,684 TESS targets (Ziegler et al. 2026, arXiv:2608.26250,
    Table 2: 0.1" 2.04 [1.90], 0.2" 4.04 [3.43], 1" 4.75 [3.95], 1.5" 5.28 [4.31], 3" 6.19 [4.91]); inside 0.1" linear in
    log(sep) to 0 mag at 0.04" (ASSUMED; 40-mas ACF peak);
  - Gemini-S/Zorro: inner working angles 0.017" (832 nm) and 0.026" (562 nm), contrasts of 4-5 mag near the diffraction
    limit and 5-9 mag out to 1.2", lower at 0.5-1" for V >~ 15 and ~2.5 mag worse when seeing doubles (Howell & Furlan
    2022; Howell et al. 2025); the curve shape for this V = 15.6 target is ASSUMED from these statements;
  - VLT/ERIS-NIX Ks: NGS mode needs R < 11 (target RP 13.07: not possible on-axis); LGS mode with the target as
    tip-tilt star (R 7-18 allowed) (ESO ERIS pages); the contrast curve is ASSUMED (typical LGS-AO K-band shape);
  - VLT/SPHERE: SAXO contrast degrades strongly for G > 12-12.5 (Jones et al. 2022); target G 14.28 -> not recommended;
    a degraded curve is shown for completeness; NaCo is decommissioned (ESO).
* Excluded fraction = posterior-weighted fraction of companions with band contrast below the curve at their projected
  separation (random orbital phase, single epoch). Bands: I for HRCam and Zorro-832, V for Zorro-562, Ks for ERIS/SPHERE.
* FPP after each step: r05 evaluate() with its 'hri' option (contrast array substituted for the right band), 'chroma'
  (two-band depth ratio at 5%/10%) and pair A resolved (lr_pairA 1e-8), and f_BD = 0 for RV; REF and INV light curves,
  as published and sign-fixed.
Outputs: out/hri_exclusion.csv, out/fpp_after_followup.csv, out/hri.json, out/ptp_draw_weights.npz"""
import sys, os
sys.dont_write_bytecode = True
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'): os.environ[k] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json
sys.path.insert(1, f'{B}/referee2_fpp/scripts')
import r05_fpp as R
CP = R.CP
dV = -2.5*np.log10(np.clip(CP['DcV'], 1e-12, None)/np.clip(1 - CP['DcV'], 1e-12, None))
CURVES = {   # name: (band, [(sep", dmag)], note)
 'Zorro 832 nm (nominal)': ('I', [(0.017, 0.0), (0.025, 2.5), (0.05, 3.8), (0.1, 4.3), (0.2, 4.8), (0.5, 5.2), (1.2, 5.5)], 'ASSUMED shape from Howell & Furlan 2022 statements'),
 'Zorro 832 nm (poor seeing)': ('I', [(0.02, 0.0), (0.03, 2.0), (0.05, 3.0), (0.1, 3.5), (0.2, 3.8), (0.5, 4.0), (1.2, 4.2)], 'nominal minus ~1-1.3 mag, larger IWA'),
 'Zorro 562 nm (nominal)': ('V', [(0.026, 0.0), (0.04, 2.0), (0.1, 3.5), (0.2, 4.0), (0.5, 4.5), (1.2, 4.8)], 'ASSUMED; faint 562-nm targets are often not detected'),
 'SOAR HRCam I (median of 2684 TOIs)': ('I', [(0.04, 0.0), (0.1, 2.04), (0.2, 4.04), (1.0, 4.75), (1.5, 5.28), (3.0, 6.19)], 'Ziegler+2026 Table 2 median'),
 'SOAR HRCam I (16th percentile)': ('I', [(0.04, 0.0), (0.1, 1.90), (0.2, 3.43), (1.0, 3.95), (1.5, 4.31), (3.0, 4.91)], 'Ziegler+2026 Table 2 16th pct (fainter target)'),
 'SOAR HRCam I (faint, pessimistic)': ('I', [(0.06, 0.0), (0.1, 1.4), (0.2, 2.9), (1.0, 3.45), (1.5, 3.8), (3.0, 4.4)], '16th pct - 0.5 mag, IWA 0.06"'),
 'VLT/ERIS-NIX Ks LGS (nominal)': ('Ks', [(0.05, 0.0), (0.06, 1.5), (0.1, 2.5), (0.2, 4.0), (0.5, 5.5), (1.0, 6.5), (2.0, 7.0)], 'ASSUMED typical LGS-AO K curve'),
 'VLT/ERIS-NIX Ks LGS (pessimistic)': ('Ks', [(0.07, 0.0), (0.1, 1.5), (0.2, 3.0), (0.5, 4.5), (1.0, 5.5), (2.0, 6.0)], 'ASSUMED'),
 'VLT/SPHERE IRDIS (degraded, G=14.3)': ('Ks', [(0.1, 0.0), (0.15, 2.0), (0.3, 4.0), (0.5, 5.0), (1.0, 6.0)], 'ASSUMED; SAXO degraded for G > 12.5'),
}
def lim(curve, sep):
    s, v = map(np.array, zip(*curve))
    return np.where(sep >= s[0], np.interp(np.log10(sep), np.log10(s), v), -99.0)
bandval = {'I': CP['dI'], 'V': dV, 'Ks': CP['dKs']}
orig_ecc = R.ecc
def ecc_fixed(kind, n=2500):
    ge, geo = orig_ecc(kind, n); return 1/ge, geo
def set_sign(fixed):
    R.ZT.clear(); R.ZCACHE.clear(); R.ecc = ecc_fixed if fixed else orig_ecc; R.rng = np.random.default_rng(11)
def weights(p):
    """per-draw evidence weights of PTP, BD-on-companion and TP-diluted (copied from r05 evaluate / r07)"""
    tab = R.build(p['ds'], p['ecc'], p['alpha'])
    R.evaluate(p)                                   # fills ZCACHE
    Z1, Z2 = R.ZCACHE[(p['ds'], p['ecc'], p['alpha'])]
    Lc = CP['L_ipd']*CP['L_ast']*CP['L_dec']*CP['L_wide']
    Lcmd = R.cmd_like(p, R.DX_C, R.DR_C)
    q = CP['q']; qmin = 0.075/CP['m1']
    wq = ((1 - p['F_twin'])*np.ones_like(q) + p['F_twin']*(q >= 0.95)/0.05)*(1 - qmin)
    Sa = np.where(CP['a'] <= 10, p['S_close'], np.where(CP['a'] >= 200, 1.0, p['S_close'] + (1 - p['S_close'])*np.log(CP['a']/10)/np.log(20)))
    occ = lambda cl, M: p['f_' + cl]/R.LNP*(np.asarray(M)/0.37)**p['beta_' + cl]
    PL = ('sub', 'N', 'S', 'J')
    p2 = sum(occ(cl, CP['m2'])*Z2[cl] for cl in PL)*Sa*CP['stab2']
    p1 = sum(occ(cl, CP['m1'])*Z1[cl] for cl in PL)*Sa*CP['stab1']
    b2 = occ('BD', CP['m2'])*Z2['BD']*CP['stab2']
    base = wq*Lc*Lcmd
    twin = CP['dT'] <= p['dT_twin']
    return dict(PTP=base*p2*~twin, BDHEB=base*b2, TPdil=base*p1, config=base)
out_ex, wsave = [], {}
for fixed in (False, True):
    set_sign(fixed)
    for ds in ('ref_beta', 'inv_beta'):
        p = dict(R.BASE); p['ds'] = ds
        W = weights(p)
        tag = f"{ds}|{'signfixed' if fixed else 'as_published'}"
        wsave[tag.replace('|', '__')] = W['PTP'].astype(np.float32)
        for name, (bnd, curve, note) in CURVES.items():
            det = bandval[bnd] < lim(curve, CP['sep'])
            row = dict(case=tag, instrument=name, band=bnd, note=note)
            for k, w in W.items():
                row[f'excl_{k}'] = float(np.sum(w*det)/np.sum(w))
            row['PTP_sep_below_IWA'] = float(np.sum(W['PTP']*(CP['sep'] < curve[0][0]))/np.sum(W['PTP']))
            out_ex.append(row)
        pct = lambda v, w: [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(w[np.argsort(v)])/w.sum(), np.sort(v))]
        print(tag, 'PTP sep', pct(CP['sep'], W['PTP']), 'dI', pct(CP['dI'], W['PTP']), 'dV', pct(dV, W['PTP']), 'dKs', pct(CP['dKs'], W['PTP']), flush=True)
EX = pd.DataFrame(out_ex); EX.to_csv(f'{OUT}/hri_exclusion.csv', index=False)
np.savez_compressed(f'{OUT}/ptp_draw_weights.npz', **wsave, sep=CP['sep'].astype(np.float32), dT=CP['dT'].astype(np.float32),
                    dI=CP['dI'].astype(np.float32), dV=dV.astype(np.float32), dKs=CP['dKs'].astype(np.float32))
print(EX[['case', 'instrument', 'excl_PTP', 'excl_BDHEB', 'excl_TPdil', 'PTP_sep_below_IWA']].to_string())
# ---------------- FPP after follow-up steps ----------------
def hri_param(name):
    bnd, curve, _ = CURVES[name]; return bnd, tuple(zip(*curve))
def run(p_extra, band='I'):
    saved = CP['dI']
    if band != 'I': CP['dI'] = bandval[band]
    try:
        p = dict(R.BASE); p.update(p_extra); r = R.evaluate(p)
    finally:
        CP['dI'] = saved
    po = r['post']
    return dict(FPP=r['FPP'], FPP_strict=r['FPP_strict'], NFPP=r['NFPP'], PTP=po['PTP_companion'],
                BD=po['BD_on_M3V_single'] + po['BD_on_primary_of_binary'] + po['BDHEB_on_companion'], WD=po['EB_WD_on_M3V'],
                HEB=po['HEB_MSMS'], BG=po['BG_all'], dominant=r['dominant_FP'])
STEPS = [('0 none', {}, 'I'),
         ('1a speckle Zorro 832 (nominal)', dict(hri=hri_param('Zorro 832 nm (nominal)')[1]), 'I'),
         ('1b speckle SOAR HRCam (16th pct)', dict(hri=hri_param('SOAR HRCam I (16th percentile)')[1]), 'I'),
         ('1c AO ERIS Ks (nominal)', dict(hri=hri_param('VLT/ERIS-NIX Ks LGS (nominal)')[1]), 'Ks'),
         ('1d referee2 generic HRI curve', dict(hri=([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0], [2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0])), 'I'),
         ('2a two-band 10% + pair resolved', dict(chroma=0.10, lr_pairA=1e-8), 'I'),
         ('2b two-band 5% + pair resolved', dict(chroma=0.05, lr_pairA=1e-8), 'I'),
         ('3 Zorro + two-band 5%', dict(hri=hri_param('Zorro 832 nm (nominal)')[1], chroma=0.05, lr_pairA=1e-8), 'I'),
         ('3b HRCam(16th) + two-band 5%', dict(hri=hri_param('SOAR HRCam I (16th percentile)')[1], chroma=0.05, lr_pairA=1e-8), 'I'),
         ('4 Zorro + two-band 5% + low-precision RV (BD/WD channels closed)', dict(hri=hri_param('Zorro 832 nm (nominal)')[1], chroma=0.05, lr_pairA=1e-8, f_BD=0.0, r_wd=0.0), 'I'),
         ('4b two-band 5% + low-precision RV (no imaging)', dict(chroma=0.05, lr_pairA=1e-8, f_BD=0.0, r_wd=0.0), 'I'),
         ('5 BD prior from ../bd (f_BD 1.3e-4), no follow-up', dict(f_BD=1.3e-4), 'I'),
         ('5b BD prior 1.3e-4 + Zorro + two-band 5%', dict(f_BD=1.3e-4, hri=hri_param('Zorro 832 nm (nominal)')[1], chroma=0.05, lr_pairA=1e-8), 'I')]
rows = []
for fixed in (False, True):
    set_sign(fixed)
    for ds in ('ref_beta', 'inv_beta'):
        for name, kw, bnd in STEPS:
            kw2 = dict(kw); kw2['ds'] = ds
            r = run(kw2, bnd)
            rows.append(dict(sign='fixed' if fixed else 'as_published', lc=ds, step=name, **r))
            print(f"{'fixed' if fixed else 'publ.'} {ds} {name:65s} FPP {r['FPP']:.2e} strict {r['FPP_strict']:.2e} PTP {r['PTP']:.1e} BD {r['BD']:.1e} WD {r['WD']:.1e} dom {r['dominant']}", flush=True)
FP = pd.DataFrame(rows)
# precise-RV step: a K ~ 10-20 m/s detection on the M3V in phase with the ephemeris removes PTP, BD, WD, HEB and pair-A channels;
# what is left is the background term (carried at 3e-5 posterior) -> FPP <~ 1e-4 (bookkeeping, not an r05 evaluation)
for (sg, ds), s in FP.groupby(['sign', 'lc']):
    base = s[s.step.str.startswith('3 ')].iloc[0]
    rows.append(dict(sign=sg, lc=ds, step='6 + precise RV mass on the M3V (K detected at >= 5 sigma)', FPP=float(base.BG + 1e-5), FPP_strict=np.nan,
                     NFPP=0.0, PTP=0.0, BD=0.0, WD=0.0, HEB=0.0, BG=float(base.BG), dominant='background (carried 3e-5) + residual twin-host ambiguity'))
FP = pd.DataFrame(rows); FP.to_csv(f'{OUT}/fpp_after_followup.csv', index=False)
jdump(dict(curves={k: dict(band=v[0], points=v[1], note=v[2]) for k, v in CURVES.items()},
           exclusion=EX.to_dict('records'), fpp_steps=FP.to_dict('records'),
           notes=['excl_PTP = posterior-weighted fraction of the PTP companions that the curve would detect (single epoch)',
                  'sign-fixed rows follow ../bd/b10; BD channel priors in r05 are 1e-3 per star unless step 5 (1.3e-4 from ../bd)',
                  'RV step 6 is bookkeeping: background scenarios are carried at 3e-5 in r05']),
      f'{OUT}/hri.json')
