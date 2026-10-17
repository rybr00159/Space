"""f08: planet mass predictions, RV semi-amplitudes, RV feasibility (precise and low-precision), TSM/ESM and ranking.
Mass-radius (verified from the papers): Chen & Kipping 2017 Table 1 (C = 1.008 R_E, Terran S = 0.2790, Neptunian
S = 0.589 with 14.6% radius dispersion, T1 = 2.04 M_E); Otegi et al. 2020 volatile-rich M = (1.74 +- 0.38) R^(1.58 +- 0.10);
Kempton+18 eq. 2 (Louie+18 implementation M = 1.436 R^1.70, used inside TSM); empirical: NASA Exoplanet Archive
transiting planets 4-7 R_E with measured masses around Teff < 4000 K hosts (refreshed here by TAP GET; ctx copy as fallback).
Planet radius 4.90 +- 0.28 R_E (clean/out/final_depth.json, stat+sys), M* 0.370 +- 0.011, R* 0.383 +- 0.012, Teff 3350 +- 55,
P 5.6915020 d, i 88.8 deg.
RV precision anchors (verified text): TOI-756 (V 14.6, J 11.1; Parc et al. 2025, arXiv:2510.14927): HARPS 2400-s photon noise
median 5.5 m/s; NIRPS 3x800 s photon noise 8.4 m/s (15.4 m/s after removing lines hit by OH/telluric overlap because
|Vsys - BERV| <= 10 km/s); fitted jitter ~15 m/s for both. TOI-674 (Murgas et al. 2021, arXiv:2106.01246): HARPS 30-min
exposures, median 7.3 m/s. Scaled to this star photon-limited (10^(0.2 dm)); ESPRESSO scaled from HARPS by collecting
area (8.2/3.57)^2 and an ASSUMED 1.25x throughput ratio.
TSM/ESM: Kempton et al. 2018 eqs. 1, 3, 4 (scale factor 1.15 for 4-10 R_E; thresholds TSM > 90 for 1.5-10 R_E, ESM 7.5).
Outputs: out/rv_feasibility.csv, out/berv_season.csv, out/tsm_esm_rank.csv, out/f08_rv_atmo.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, json, io, urllib.parse, urllib.request
rng = np.random.default_rng(5)
N = 200000
Rp = rng.normal(4.90, 0.28, N); Ms = rng.normal(0.370, 0.011, N); Rs = rng.normal(0.383, 0.012, N); Ts = rng.normal(3350, 55, N)
P = 5.6915020; inc = np.radians(88.8)
G = 6.674e-11; MSUN_ = 1.989e30; MEARTH = 5.972e24; RSUN_ = 6.957e8; AU = 1.496e11
def K_of(Mp_E, Ms_sun, e=0.0):
    return (2*np.pi*G/(P*86400))**(1/3)*Mp_E*MEARTH*np.sin(inc)/((Ms_sun*MSUN_ + Mp_E*MEARTH)**(2/3))/np.sqrt(1 - e**2)
# Chen & Kipping
R_T1 = 1.008*2.04**0.279
M_ck = 2.04*(Rp*np.exp(rng.normal(0, np.log(1.146), N))/R_T1)**(1/0.589)
# Otegi volatile-rich (parameter uncertainties only)
M_ot = rng.normal(1.74, 0.38, N).clip(0.3)*Rp**rng.normal(1.58, 0.10, N)
M_lo = 1.436*Rp**1.70
# empirical from the NASA Exoplanet Archive
q = ("select pl_name,hostname,pl_orbper,pl_rade,pl_radeerr1,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_bmasselim,pl_bmassprov,pl_orbsmax,pl_eqt,"
     "st_teff,st_rad,st_mass,sy_jmag,sy_kmag,sy_dist from pscomppars where tran_flag=1 and pl_rade>3.5 and pl_rade<8 and pl_orbper<15 and st_teff<4000")
url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?' + urllib.parse.urlencode(dict(query=q, format='csv'))
try:
    txt = urllib.request.urlopen(url, timeout=120).read().decode()
    NEA = pd.read_csv(io.StringIO(txt)); nea_src = 'NASA Exoplanet Archive pscomppars TAP (GET), fetched 2026-10-17'
    NEA.to_csv(f'{F}/data/nea_mdwarf_neptunes_refresh.csv', index=False)
except Exception as ex:
    NEA = pd.read_csv(f'{B}/ctx/out/known_mdwarf_neptunes.csv'); nea_src = f'ctx copy (refresh failed: {ex})'
lim_ = NEA['pl_bmasselim'].fillna(0) if 'pl_bmasselim' in NEA else pd.Series(0, index=NEA.index)
relerr = NEA.pl_bmasseerr1.abs()/NEA.pl_bmasse
emp = NEA[(NEA.pl_rade >= 4) & (NEA.pl_rade <= 7) & NEA.pl_bmasse.notna() & (lim_ == 0) & (relerr < 0.5) &
          (NEA.get('pl_bmassprov', pd.Series(['Mass']*len(NEA))) == 'Mass')]
lm = np.log(emp.pl_bmasse.values); lr = np.log(emp.pl_rade.values)
A = np.vstack([np.ones_like(lr), lr - np.log(4.9)]).T
coef, res, *_ = np.linalg.lstsq(A, lm, rcond=None)
sc = float(np.sqrt(np.sum((lm - A @ coef)**2)/max(len(lm) - 2, 1)))
M_emp = np.exp(coef[0] + coef[1]*(np.log(Rp) - np.log(4.9)) + rng.normal(0, sc, N))
def summ(M):
    K = K_of(M, Ms)
    return dict(M_16_50_84=[float(x) for x in np.percentile(M, [16, 50, 84])], K_16_50_84=[float(x) for x in np.percentile(K, [16, 50, 84])],
                K_2p5_97p5=[float(x) for x in np.percentile(K, [2.5, 97.5])])
masses = {'Chen&Kipping17 (Neptunian, 14.6% R scatter)': summ(M_ck), 'Otegi+20 volatile-rich': summ(M_ot),
          'Kempton18/Louie18 M=1.436R^1.70': summ(M_lo),
          f'empirical NEA M-dwarf 4-7 R_E (n={len(emp)}, lnM slope {coef[1]:.2f}, scatter {sc:.2f} dex_e)': summ(M_emp)}
for k, v in masses.items(): print(k, v)
# BD / WD semi-amplitudes for the low-precision RV design
MJ = 317.83
bdK = {f'{m} M_J': float(K_of(m*MJ, 0.370)) for m in (13, 20, 40, 60, 80)}
wdK = float(K_of(0.6*332946.0, 0.370))
print('BD K', bdK, 'WD K', wdK)
# ---------------- RV precision per 2400 s on this star, and epochs for a 5-sigma K ----------------
dV_756, dJ_756 = 15.61 - 14.6, 11.36 - 11.1
harps = 5.5*10**(0.2*dV_756)                    # photon-limited scaling from TOI-756
harps_674 = 7.3*10**(0.2*(15.61 - 14.21))*np.sqrt(1800/2400)  # TOI-674 V = 14.21 (assumed from the paper's photometry table)
nirps = 8.4*10**(0.2*dJ_756); nirps_bad = 15.4*10**(0.2*dJ_756)
esp = harps/np.sqrt((8.2/3.57)**2*1.25)
INSTR = [  # name, site lat, sigma per 2400 s (m/s), basis
 ('ESPRESSO (VLT, 1-UT HR)', -24.63, esp, 'HARPS-scaled by area x1.25 throughput (ASSUMED)'),
 ('NIRPS (3.6m, alone)', -29.26, nirps, 'TOI-756 photon noise scaled in J'),
 ('NIRPS (OH-affected epochs)', -29.26, nirps_bad, 'TOI-756 after line removal, scaled'),
 ('HARPS (3.6m)', -29.26, harps, 'TOI-756 photon noise scaled in V (TOI-674 scaling gives %.1f)' % harps_674),
 ('NIRPS + HARPS simultaneous', -29.26, 1/np.sqrt(1/nirps**2 + 1/harps**2), 'inverse-variance combination'),
]
rows = []
for name, lat, s, basis in INSTR:
    for jit in (5.0, 10.0, 15.0):
        seff = np.hypot(s, jit)
        for Kc, Klab in ((float(np.percentile(K_of(M_ck, Ms), 50)), 'K_CK_median'), (float(np.percentile(K_of(M_ck, Ms), 16)), 'K_CK_16th')):
            n5 = int(np.ceil(2*(5*seff/Kc)**2))
            rows.append(dict(instrument=name, sigma_2400s=round(s, 1), jitter=jit, sigma_eff=round(seff, 1), K_case=Klab, K=round(Kc, 1),
                             epochs_5sigma=n5, hours_on_target=round(n5*2400/3600*1.15, 1), basis=basis,
                             max_alt_deg=round(90 - abs(lat - TARGET['dec']), 1)))
RVF = pd.DataFrame(rows); RVF.to_csv(f'{OUT}/rv_feasibility.csv', index=False)
print(RVF[(RVF.K_case == 'K_CK_median')][['instrument', 'sigma_2400s', 'jitter', 'sigma_eff', 'K', 'epochs_5sigma', 'hours_on_target']].to_string())
# instruments that cannot reach Dec -51.3
NORTH = {'MAROON-X / KPF / SPIRou / IRD (Maunakea)': 19.82, 'CARMENES (Calar Alto)': 37.22, 'NEID (Kitt Peak)': 31.96,
         'HARPS-N / GIANO-B (La Palma)': 28.76, 'HPF (HET, fixed elevation 55 deg, Dec -10.3..+71.6)': 30.68, 'SOPHIE (OHP)': 43.93,
         'ESPRESSO/NIRPS/HARPS/CRIRES+ (Chile)': -24.63, 'PFS / MIKE (Magellan)': -29.01, 'Veloce (AAT)': -31.27, 'SALT HRS (fixed el. 53 deg; Dec +10..-75)': -32.38,
         'FEROS (MPG 2.2m)': -29.26, 'GHOST (Gemini-S)': -30.24, 'CHIRON (CTIO 1.5m)': -30.17}
reach = {k: dict(lat=v, max_alt_deg=round(90 - abs(v - TARGET['dec']), 1), min_airmass=round(1/np.sin(np.radians(max(90 - abs(v - TARGET['dec']), 1))), 2)) for k, v in NORTH.items()}
reach['HPF (HET, fixed elevation 55 deg, Dec -10.3..+71.6)']['note'] = 'cannot point to Dec -51'
# ---------------- BERV over the observable season (NIRPS OH/telluric overlap test) ----------------
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
loc = EarthLocation.from_geodetic(-70.738*u.deg, -29.2563*u.deg, 2400*u.m)
c = target_coord(2027.5)
from astropy.coordinates import SkyCoord
cc = SkyCoord(ra=c.ra, dec=c.dec)
days = Time('2027-01-01T04:00:00') + np.arange(0, 730, 5)*u.day
from astropy.coordinates import AltAz, get_sun
br = []
for t in days:
    hrs = t + np.arange(-6, 6.01, 0.5)*u.hour
    fr = AltAz(obstime=hrs, location=loc)
    alt = cc.transform_to(fr).alt.deg; sun = get_sun(hrs).transform_to(fr).alt.deg
    ok = (alt > 30) & (sun < -12)
    if ok.sum() == 0: continue
    berv = cc.radial_velocity_correction('barycentric', obstime=hrs[ok][len(hrs[ok])//2], location=loc).to(u.km/u.s).value
    br.append(dict(date=t.isot[:10], hours_observable=float(ok.sum()*0.5), berv_kms=float(berv), v_topo_kms=float(16.28 - berv)))
BR = pd.DataFrame(br); BR.to_csv(f'{OUT}/berv_season.csv', index=False)
frac_bad = float(np.mean(np.abs(BR.v_topo_kms) <= 10))
BR['month'] = BR.date.str[:7]
bad_months = sorted(set(BR[np.abs(BR.v_topo_kms) <= 10].month)); good_months = sorted(set(BR[np.abs(BR.v_topo_kms) > 10].month))
print('NIRPS OH-overlap months', bad_months, 'good months', good_months)
print('BERV range', BR.berv_kms.min(), BR.berv_kms.max(), 'fraction of observable dates with |Vsys-BERV|<=10 km/s:', frac_bad)
# ---------------- TSM / ESM ----------------
def B75(T):
    h, cc_, k = 6.626e-34, 2.998e8, 1.381e-23; l = 7.5e-6
    return 1/np.expm1(h*cc_/(l*k*T))
def tsm(Rp_, Mp_, Rs_, Teq_, J): return 1.15*Rp_**3*Teq_/(Mp_*Rs_**2)*10**(-J/5)
def esm(Rp_, Rs_, Teq_, Ts_, K): return 4.29e6*B75(1.10*Teq_)/B75(Ts_)*(Rp_*6.371e6/(Rs_*6.957e8))**2*10**(-K/5)
a = (G*Ms*MSUN_*(P*86400)**2/(4*np.pi**2))**(1/3)
Teq = Ts*np.sqrt(Rs*RSUN_/a)*0.25**0.25
TSMd = {k: [float(x) for x in np.percentile(tsm(Rp, M, Rs, Teq, 11.363), [16, 50, 84])] for k, M in (('CK17', M_ck), ('Otegi20', M_ot), ('Louie18', M_lo), ('empirical', M_emp))}
ESMd = [float(x) for x in np.percentile(esm(Rp, Rs, Teq, Ts, 10.482), [16, 50, 84])]
Teqd = [float(x) for x in np.percentile(Teq, [16, 50, 84])]
print('Teq', Teqd, 'TSM', TSMd, 'ESM', ESMd)
# known planets with uniform formulas
K2 = NEA.copy()
for col in ('pl_orbsmax', 'pl_eqt'):
    if col not in K2: K2[col] = np.nan
aK = (G*K2.st_mass*MSUN_*(K2.pl_orbper*86400)**2/(4*np.pi**2))**(1/3)
K2['Teq_calc'] = K2.st_teff*np.sqrt(K2.st_rad*RSUN_/aK)*0.25**0.25
limk = K2['pl_bmasselim'].fillna(0) if 'pl_bmasselim' in K2 else 0
K2['M_used'] = np.where(K2.pl_bmasse.notna() & (limk == 0), K2.pl_bmasse, 1.436*K2.pl_rade**1.70)
K2['TSM'] = tsm(K2.pl_rade, K2.M_used, K2.st_rad, K2.Teq_calc, K2.sy_jmag)
K2['ESM'] = esm(K2.pl_rade, K2.st_rad, K2.Teq_calc, K2.st_teff, K2.sy_kmag)
cand = dict(pl_name='TIC 214311908.01 (candidate)', hostname='TIC 214311908', pl_orbper=P, pl_rade=4.90, pl_bmasse=float(np.median(M_ck)),
            st_teff=3350, st_rad=0.383, st_mass=0.370, sy_jmag=11.363, sy_kmag=10.482, sy_dist=64.31, Teq_calc=float(np.median(Teq)),
            M_used=float(np.median(M_ck)), TSM=TSMd['CK17'][1], ESM=ESMd[1])
K2 = pd.concat([K2, pd.DataFrame([cand])], ignore_index=True).sort_values('TSM', ascending=False)
sel = K2[(K2.pl_rade >= 4) & (K2.pl_rade <= 7) & (K2.pl_orbper < 10)]
sel.to_csv(f'{OUT}/tsm_esm_rank.csv', index=False)
rank_tsm = int(np.where(sel.pl_name.values == cand['pl_name'])[0][0]) + 1
rank_esm = int(np.where(sel.sort_values('ESM', ascending=False).pl_name.values == cand['pl_name'])[0][0]) + 1
print(sel[['pl_name', 'pl_rade', 'M_used', 'st_teff', 'sy_jmag', 'Teq_calc', 'TSM', 'ESM']].round(1).to_string())
print('rank TSM', rank_tsm, 'of', len(sel), 'rank ESM', rank_esm)
jdump(dict(mass_K=masses, nea_source=nea_src, n_empirical=int(len(emp)), empirical_fit=dict(lnM_at_4p9=float(coef[0]), slope=float(coef[1]), scatter=sc,
                                                                                               planets=emp.pl_name.tolist()),
           BD_K_kms={k: v/1e3 for k, v in bdK.items()}, WD_K_kms=wdK/1e3,
           rv_precision_2400s=dict(ESPRESSO=esp, NIRPS=nirps, NIRPS_OH_affected=nirps_bad, HARPS=harps, HARPS_from_TOI674=harps_674),
           instrument_reach=reach, berv=dict(min=float(BR.berv_kms.min()), max=float(BR.berv_kms.max()), frac_dates_vtopo_le10=frac_bad, vsys=16.28, months_vtopo_le10=bad_months, months_ok=good_months),
           Teq_K=Teqd, TSM=TSMd, ESM=ESMd, rank_TSM=rank_tsm, rank_ESM=rank_esm, n_compared=int(len(sel)),
           thresholds=dict(TSM_1p5_10=90, ESM=7.5)), f'{OUT}/f08_rv_atmo.json')
