"""s16: independent (non-light-curve) tests of the 'white dwarf totally occulted by the M3V' scenario.
(a) Gaia astrometry: RUWE simulation of the referee_star r10 type (41 FoV transits x 9 CCDs, random scan angles,
    5-parameter fit, RUWE0 drawn from comparison M dwarfs G 13.5-15.5 from ../referee_star/data/mdw_cmd.csv) for a
    dark companion of mass M2 at P = 5.6915 d (a from Kepler, M1 = 0.37), photocentre semi-major axis
    alpha = a (q/(1+q) - eps), eps = WD share of the G-band flux (from (b)); sigma_CCD = 0.22 and 0.35 mas.
    Reported: L = P(RUWE_sim <= 1.10)/P(RUWE0 <= 1.10) (target RUWE 1.076).
(b) Blue/UV excess: a WD supplying f_A of the TESS flux at the fitted (M, R) has the Teff that gives this ratio
    (Kurucz logg 5 surface fluxes, specs.py); its SkyMapper v, g and DECam g share of the M3V's observed flux
    -> predicted colour excess in v - g and g - i. Compared with an empirical locus: southern M dwarfs with
    BP-RP 2.66-2.86 from mdw_cmd.csv cross-matched with SkyMapper DR4 (CDS XMatch, II/379/smssdr4, 3").
(c) Radial velocity: K1 of the M3V for each WD mass; Gaia DR3 combined-spectrum RV error of the M3V vs comparison
    stars of the same grvs_mag (percentile).
Output: out/wd_checks.json, data/sm4_xmatch.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import G, MSUN, RSUN, DAY, AU, SD
from specs import S, absmag
rng = np.random.default_rng(3)
RS = 'tic214311908/analysis/referee_star'
P = 5.6915087; M1 = 0.370; R1 = 0.3813; T1 = 3350; d = 64.31; DM = 5*np.log10(d/10)
fA = 0.0155  # WD share of the TESS flux in the s06 fits
def r_nauenberg(M):
    mu = M/1.44; return 0.0112*np.sqrt(mu**(-2/3.) - mu**(2/3.))
# M3V observed fluxes (AB) in v, g (SkyMapper), DECam g, and its TESS-band surface flux from BT-Settl 3350 K
obs = dict(smv=18.192, smg=15.817, g=16.177, smi=13.387)
out = dict(fA=fA, cases={})
SB_M = S(T1, 'T')*R1**2
Tg = np.exp(np.linspace(np.log(4000), np.log(50000), 400))
for M2 in [0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2]:
    fac = 1.8 if M2 < 0.45 else 1.1
    RA = r_nauenberg(M2)*fac
    ratio = fA/(1 - fA)*SB_M/RA**2
    sT = np.array([S(t, 'T', True) for t in Tg])
    TA = float(np.interp(ratio, sT, Tg)) if ratio < sT.max() else np.nan
    if not np.isfinite(TA):
        out['cases'][str(M2)] = dict(note='no Teff <= 50 kK gives the required TESS flux'); continue
    mv = absmag(TA, RA, 'smv', True) + DM; mg = absmag(TA, RA, 'smg', True) + DM; mdg = absmag(TA, RA, 'g', True) + DM
    mG_ = absmag(TA, RA, 'r', True) + DM  # proxy for Gaia G share (DECam r)
    fv = 10**(-0.4*(mv - obs['smv'])); fg = 10**(-0.4*(mg - obs['smg'])); fdg = 10**(-0.4*(mdg - obs['g']))
    eps = 10**(-0.4*(mG_ - 14.737))  # share of the DECaPS2 r flux of the M3V, as the G-band proxy
    a = (G*(M1+M2)*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/AU
    K1 = (2*np.pi*G/(P*DAY))**(1/3.)*M2*MSUN/((M1+M2)*MSUN)**(2/3.)/1e5
    out['cases'][str(M2)] = dict(R_WD=RA, Teff_WD_needed=TA, a_au=a, K1_kms=K1, WD_v_mag=mv, WD_smg_mag=mg, WD_decam_g_mag=mdg,
                                 frac_v=fv/(1+fv), frac_smg=fg/(1+fg), frac_decam_g=fdg/(1+fdg),
                                 dmag_v=-2.5*np.log10(1+fv), dmag_smg=-2.5*np.log10(1+fg), dmag_v_minus_g=-2.5*np.log10((1+fv)/(1+fg)), eps_G=eps)
    print(M2, {k: round(v, 4) for k, v in out['cases'][str(M2)].items()}, flush=True)
# ---------------- (a) astrometry
m = pd.read_csv(f'{RS}/data/mdw_cmd.csv')
r0 = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.25) & (m.ipd_frac_multi_peak <= 2)].ruwe.values
pR0 = np.mean(r0 <= 1.10)
plx = 15.55; NTR = 41; NCCD = 9
def sim(M2, eps, sig, nsim=400):
    q = M2/M1; beta = q/(1+q) - eps/(1+eps)
    a = (M1+M2)**(1/3.)*(P/365.25)**(2/3.)
    alpha = beta*a*plx; Py = P/365.25
    res = []
    for k in range(nsim):
        t = np.sort(rng.uniform(2014.6, 2017.4, NTR)); psi = rng.uniform(0, 2*np.pi, NTR)
        lam = 2*np.pi*(t - 2015.0); pf = np.sin(lam)*np.sin(psi) + 0.6*np.cos(lam)*np.cos(psi)
        A = np.column_stack([np.sin(psi), np.cos(psi), np.sin(psi)*(t-2016), np.cos(psi)*(t-2016), pf])
        inc = np.arccos(rng.uniform(-1, 1)); Om = rng.uniform(0, 2*np.pi); ph = rng.uniform(0, 2*np.pi)
        u = 2*np.pi*t/Py + ph
        x, y = np.cos(u), np.sin(u)*np.cos(inc)
        dE = alpha*(x*np.cos(Om) - y*np.sin(Om)); dN = alpha*(x*np.sin(Om) + y*np.cos(Om))
        al = dE*np.sin(psi) + dN*np.cos(psi)
        rr = al - A @ np.linalg.lstsq(A, al, rcond=None)[0]
        chi = NCCD*np.sum(rr**2)/sig**2
        res.append(np.sqrt(rng.choice(r0)**2 + chi/(NTR*NCCD - 5)))
    res = np.array(res)
    return float((np.mean(res <= 1.10) + 1e-3)/pR0), float(np.median(res)), float(alpha)
out['astrometry'] = dict(n_RUWE0=int(len(r0)), P_RUWE0_le_1p10=float(pR0))
for M2 in [0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2]:
    eps = out['cases'].get(str(M2), {}).get('eps_G', 0.02)
    for sig in (0.22, 0.35):
        L, med, al = sim(M2, eps, sig)
        out['astrometry'][f'M{M2}_sig{sig}'] = dict(L=L, median_ruwe=med, alpha_mas=al)
    print('astrometry', M2, out['astrometry'][f'M{M2}_sig0.22'], out['astrometry'][f'M{M2}_sig0.35'], flush=True)
# ---------------- (c) RV error percentile
sel = m[(m.grvs_mag > 12.3) & (m.grvs_mag < 12.8) & m.radial_velocity_error.notna()]
out['rv'] = dict(target_rv_err=2.217, n_comp=int(len(sel)), percentile=float(np.mean(sel.radial_velocity_error < 2.217)*100),
                 comp_median=float(sel.radial_velocity_error.median()), comp_p90=float(sel.radial_velocity_error.quantile(0.9)),
                 comp_p99=float(sel.radial_velocity_error.quantile(0.99)))
print('rv', out['rv'])
# ---------------- (b) empirical SkyMapper locus
try:
    from astroquery.xmatch import XMatch
    from astropy.table import Table
    import astropy.units as u
    s = m[(m.dec < -5) & (m.bp_rp > 2.66) & (m.bp_rp < 2.86) & (m.ruwe < 1.4) & (m.phot_g_mean_mag > 12.5) & (m.phot_g_mean_mag < 16)]
    s = s.sample(min(len(s), 3000), random_state=1)
    tab = Table.from_pandas(s[['source_id', 'ra', 'dec', 'bp_rp', 'phot_g_mean_mag']])
    xm = XMatch.query(cat1=tab, cat2='vizier:II/379/smssdr4', max_distance=3*u.arcsec, colRA1='ra', colDec1='dec').to_pandas()
    xm.to_csv(f'{SD}/data/sm4_xmatch.csv', index=False)
    cols = {c.lower(): c for c in xm.columns}
    vcol = [c for c in xm.columns if c.lower() in ('vpsf', 'v_psf')][0]; gcol = [c for c in xm.columns if c.lower() in ('gpsf', 'g_psf')][0]
    icol = [c for c in xm.columns if c.lower() in ('ipsf', 'i_psf')][0]
    q = xm[np.isfinite(xm[vcol]) & np.isfinite(xm[gcol]) & np.isfinite(xm[icol])].copy()
    q = q.sort_values('angDist').drop_duplicates('source_id')
    q['vg'] = q[vcol] - q[gcol]; q['gi'] = q[gcol] - q[icol]
    co = np.polyfit(q.bp_rp, q.vg, 1); res = q.vg - np.polyval(co, q.bp_rp)
    sd = 1.4826*np.median(np.abs(res - np.median(res)))
    tv = obs['smv'] - obs['smg']; tres = tv - np.polyval(co, 2.756)
    co2 = np.polyfit(q.bp_rp, q.gi, 1); res2 = q.gi - np.polyval(co2, q.bp_rp); sd2 = 1.4826*np.median(np.abs(res2 - np.median(res2)))
    tres2 = (obs['smg'] - obs['smi']) - np.polyval(co2, 2.756)
    out['skymapper_locus'] = dict(n=int(len(q)), vg_fit=list(map(float, co)), vg_robust_sd=float(sd), target_vg=float(tv), target_vg_resid=float(tres),
                                  target_vg_resid_percentile=float(np.mean(res < tres)*100),
                                  gi_robust_sd=float(sd2), target_gi_resid=float(tres2), target_gi_resid_percentile=float(np.mean(res2 < tres2)*100))
    for M2, c in out['cases'].items():
        if 'dmag_v_minus_g' in c:
            c['predicted_vg_resid'] = float(tres_pred) if False else float(np.polyval(co, 2.756) + c['dmag_v_minus_g'])  # placeholder
    # what residual would an M3V+WD system show? intrinsic M3V colour = target - WD effect; we test the reverse:
    # the observed residual tres must equal (intrinsic M-dwarf scatter) + dmag_v_minus_g(WD)
    for M2, c in out['cases'].items():
        if 'dmag_v_minus_g' in c:
            c['implied_intrinsic_vg_resid'] = float(tres - c['dmag_v_minus_g'])
            c['implied_intrinsic_percentile'] = float(np.mean(res < tres - c['dmag_v_minus_g'])*100)
            c['z_vs_locus'] = float((tres - c['dmag_v_minus_g'])/sd)
    print('locus', out['skymapper_locus'])
except Exception as ex:
    out['skymapper_locus'] = dict(error=str(ex)); print('xmatch failed', ex)
json.dump(out, open(f'{SD}/out/wd_checks.json', 'w'), indent=1, default=float)
