"""f01: transit predictions 2026-11-01 .. 2028-12-31 and observability from southern sites.
Usage: f01_transits.py [ephemeris.json]   (default: ../clean/out/ephemeris_final.json, fallback ../phot/out/timing.json)
* BJD_TDB -> UTC at each site: solve t_bjd = t_obs(TDB) + ltt_bary(t_obs, site) iteratively (astropy, target position
  propagated with Gaia DR3 proper motion), then TDB -> UTC.
* 1-sigma and 3-sigma timing errors from the full (T0, P) covariance.
* Observability window: T14 (1.70 h) + 2 x 3 sigma + 2 x 1.5 h baseline; 'full' = whole window at airmass < 2
  (alt > 30 deg) with the Sun below -18 deg (astronomical darkness); also the Sun < -12 deg variant and the
  fraction of the window that is usable; Moon illumination, separation and altitude.
* Comparison with ../ctx/out/followup_transits.csv.
Outputs: out/transits_all.csv, out/obs_<site>.csv (one per site), out/obs_summary.csv, out/f01_transits.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, warnings
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body
from astropy.utils import iers
iers.conf.auto_max_age = None
iers.conf.iers_degraded_accuracy = 'warn'
warnings.filterwarnings('ignore')

eph = load_ephemeris(sys.argv[1] if len(sys.argv) > 1 else None)
print('ephemeris', eph['label'], eph['T0'], eph['P'], np.sqrt(np.diag(eph['cov'])))
tstart, tend = Time('2026-11-01T00:00:00', scale='utc'), Time('2028-12-31T23:59:59', scale='utc')
E0 = int(np.ceil((tstart.tdb.jd - 0.01 - eph['T0'])/eph['P'])) + eph['Eref']
E1 = int(np.floor((tend.tdb.jd + 0.01 - eph['T0'])/eph['P'])) + eph['Eref']
E = np.arange(E0, E1 + 1)
tb, sig = predict(eph, E)
coord = target_coord(2027.5)
coord_noPM = target_coord(2016.0)

def bjd_to_utc(tbjd, loc):
    t = Time(tbjd, format='jd', scale='tdb', location=loc)
    for _ in range(3):
        ltt = t.light_travel_time(coord, kind='barycentric', location=loc)
        t = Time(tbjd, format='jd', scale='tdb', location=loc) - ltt
    return t.utc, ltt.to(u.s).value

geo = EarthLocation.from_geocentric(0*u.m, 0*u.m, 0*u.m)
tutc_geo, ltt_geo = bjd_to_utc(tb, EarthLocation.from_geodetic(0*u.deg, 0*u.deg, 0*u.m))
rows = pd.DataFrame(dict(epoch=E, bjd_tdb=tb, sigma_min=sig*1440, sigma3_min=3*sig*1440,
                         utc_mid_geocentre=[x.isot[:19] for x in tutc_geo], bary_corr_s=ltt_geo))
# light-time checks
ltt_all = Time(np.linspace(Time('2026-01-01').jd, Time('2029-01-01').jd, 2000), format='jd', scale='tdb').light_travel_time(coord, kind='barycentric', location=geo).to(u.s).value
from astropy.coordinates import BarycentricTrueEcliptic
ecl = coord.transform_to(BarycentricTrueEcliptic(equinox='J2000'))
tdb_utc = (Time('2027-06-01', scale='utc').tdb.jd - Time('2027-06-01', scale='utc').jd)*86400
check = dict(ltt_min_s=float(ltt_all.min()), ltt_max_s=float(ltt_all.max()), ecliptic_lat_deg=float(ecl.lat.deg),
             expected_amplitude_s=float(499.005*np.cos(np.radians(ecl.lat.deg))), TDB_minus_UTC_s_2027=float(tdb_utc),
             pm_effect_on_ltt_ms=float(np.max(np.abs(Time(tb[:50], format='jd', scale='tdb').light_travel_time(coord, 'barycentric', location=geo).to(u.s).value
                                                     - Time(tb[:50], format='jd', scale='tdb').light_travel_time(coord_noPM, 'barycentric', location=geo).to(u.s).value))*1e3),
             note='UTC = BJD_TDB - ltt_bary - (TDB-UTC); TDB-UTC = 32.184 + 37 leap s (+ <2 ms periodic) assuming no new leap second before 2029')
print(check)

# ---------------- observability ----------------
half_win_h = T14_H/2 + 1.5 + 3*sig*24
rec = []
for name, (lat, lon, h) in SITES.items():
    loc = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg, h*u.m)
    tutc, ltt = bjd_to_utc(tb, loc)
    for i in range(len(E)):
        hw = half_win_h[i]
        n = int(np.ceil(2*hw*60/4)) + 1                   # 4-min grid
        off = np.linspace(-hw, hw, n)*u.hour
        tt = tutc[i] + off
        frame = AltAz(obstime=tt, location=loc)
        alt = coord.transform_to(frame).alt.deg
        sun = get_sun(tt).transform_to(frame).alt.deg
        good18 = (alt > 30) & (sun < -18)
        good12 = (alt > 30) & (sun < -12)
        intr = np.abs(off.to(u.hour).value) <= (T14_H/2 + 3*sig[i]*24)
        # baseline: usable time before ingress and after egress (within the window), at Sun < -18
        pre = (off.to(u.hour).value < -(T14_H/2 + 3*sig[i]*24)) & good18
        post = (off.to(u.hour).value > (T14_H/2 + 3*sig[i]*24)) & good18
        dt = 2*hw/(n - 1)
        tm = tutc[i]
        fr = AltAz(obstime=tm, location=loc)
        moon = get_body('moon', tm, loc)
        sunm = get_sun(tm)
        elong = sunm.separation(moon).deg
        illum = (1 - np.cos(np.radians(elong)))/2            # illuminated fraction (1 + cos i)/2 with phase angle i ~ 180 - elongation
        rec.append(dict(site=name, epoch=int(E[i]), bjd_tdb=round(float(tb[i]), 5), utc_mid=tm.isot[:16].replace('T', ' '),
                        utc_win_start=(tm - hw*u.hour).isot[:16].replace('T', ' '), utc_win_end=(tm + hw*u.hour).isot[:16].replace('T', ' '),
                        sigma_min=round(float(sig[i]*1440), 2), sigma3_min=round(float(3*sig[i]*1440), 2), window_h=round(float(2*hw), 2),
                        bary_corr_s=round(float(ltt[i]), 2),
                        full_window_sun18=bool(good18.all()), full_window_sun12=bool(good12.all()),
                        transit_covered_sun18=bool(good18[intr].all()),
                        baseline_pre_h=round(float(pre.sum()*dt), 2), baseline_post_h=round(float(post.sum()*dt), 2),
                        frac_window_sun18=round(float(good18.mean()), 3),
                        alt_mid=round(float(coord.transform_to(fr).alt.deg), 1), alt_min_window=round(float(alt.min()), 1),
                        airmass_max_window=round(float(1/np.sin(np.radians(max(alt.min(), 5)))), 2),
                        sun_max_window=round(float(sun.max()), 1),
                        moon_illum=round(float(illum), 2), moon_sep_deg=round(float(moon.separation(coord).deg), 1),
                        moon_alt_mid=round(float(moon.transform_to(fr).alt.deg), 1)))
    print(name, 'done', flush=True)
R = pd.DataFrame(rec)
R['good_partial'] = R.transit_covered_sun18 & (R.baseline_pre_h >= 0.75) & (R.baseline_post_h >= 0.75) & ~R.full_window_sun18
R.to_csv(f'{OUT}/transits_obs_all_sites.csv', index=False)
rows.to_csv(f'{OUT}/transits_all.csv', index=False)
summ = []
for name in SITES:
    s = R[R.site == name]
    fn = f"{OUT}/obs_{SITE_SHORT[name]}.csv"
    s[s.frac_window_sun18 > 0].to_csv(fn, index=False)
    full = s[s.full_window_sun18]
    for yr in ('2026-11..2027-12', '2028'):
        sel = full[(full.utc_mid < '2028') if yr.startswith('2026') else (full.utc_mid >= '2028')]
        summ.append(dict(site=name, period=yr, n_full_sun18=len(sel), n_full_sun12=int(((s.utc_mid < '2028') if yr.startswith('2026') else (s.utc_mid >= '2028'))[s.full_window_sun12].sum()),
                         n_partial_good=int(s[((s.utc_mid < '2028') if yr.startswith('2026') else (s.utc_mid >= '2028')) & s.good_partial].shape[0]),
                         dark_moon_full=int((sel.moon_illum < 0.5).sum() + ((sel.moon_illum >= 0.5) & (sel.moon_sep_deg > 60)).sum()),
                         dates_utc_mid='; '.join(f"{r.utc_mid} (E{r.epoch}, moon {r.moon_illum:.2f}/{r.moon_sep_deg:.0f}d)" for r in sel.itertuples())))
S = pd.DataFrame(summ)
S.to_csv(f'{OUT}/obs_summary.csv', index=False)
print(S[['site', 'period', 'n_full_sun18', 'n_full_sun12', 'n_partial_good', 'dark_moon_full']].to_string())

# ---------------- comparison with ctx ----------------
C = pd.read_csv(f'{B}/ctx/out/followup_transits.csv')
cmp = []
for name_ctx, name in (('CTIO', 'CTIO (LCOGT lsc, SMARTS)'), ('SAAO Sutherland', 'SAAO Sutherland (LCOGT cpt, SALT)'),
                       ('Siding Spring', 'Siding Spring (LCOGT coj, FTS/MuSCAT4, AAT)'), ('La Silla', 'La Silla (TRAPPIST-S, ExTrA, 3.6m)')):
    c = C[C.site == name_ctx]
    m = R[R.site == name].merge(c, on='epoch', suffixes=('', '_ctx'))
    if len(m) == 0: continue
    dutc = (Time(list(m.utc_mid.str.replace(' ', 'T'))).jd - Time(list(m.utc_mid_ctx.str.replace(' ', 'T'))).jd)*1440
    cmp.append(dict(site=name_ctx, n=len(m), dBJD_min_range=[float(((m.bjd_tdb - m.bjd_tdb_ctx)*1440).min()), float(((m.bjd_tdb - m.bjd_tdb_ctx)*1440).max())],
                    dUTC_min_range=[float(dutc.min()), float(dutc.max())],
                    sigma_ratio_range=[float((m.sigma_min/m.sigma_min_ctx).min()), float((m.sigma_min/m.sigma_min_ctx).max())],
                    ctx_full_sun12=int(m.full_window_ok_sun12.sum()), new_full_sun12_same_epochs=int(m.full_window_sun12.sum()),
                    epochs_full_ctx_not_new=[int(x) for x in m.epoch[m.full_window_ok_sun12 & ~m.full_window_sun12]],
                    epochs_full_new_not_ctx=[int(x) for x in m.epoch[~m.full_window_ok_sun12 & m.full_window_sun12]]))
# precise check of ctx's barycentric UTC: convert ctx's own BJD with this code and compare to its minute-rounded utc_mid
c = C[C.site == 'CTIO']
loc_ctio = EarthLocation.from_geodetic(-70.8063*u.deg, -30.1690*u.deg, 2207*u.m)
u_ctx_conv, _ = bjd_to_utc(c.bjd_tdb.values, loc_ctio)
d_conv = (u_ctx_conv.jd - Time(list(c.utc_mid.str.replace(' ', 'T'))).jd)*1440
naive = (Time(c.bjd_tdb.values, format='jd', scale='tdb').utc.jd - Time(list(c.utc_mid.str.replace(' ', 'T'))).jd)*1440
cmp.append(dict(check='ctx BJD converted here minus ctx utc_mid (minute-rounded)', min=float(d_conv.min()), max=float(d_conv.max()),
                naive_no_bary_minus_ctx_range=[float(naive.min()), float(naive.max())]))
print(cmp)
# expected BJD difference from ephemeris change
ctx_eph = dict(T0=2460840.9153, P=5.6914968)
ex = [float((predict(eph, e)[0] - (ctx_eph['T0'] + ctx_eph['P']*e))*1440) for e in (E[0], E[-1])]
jdump(dict(ephemeris=dict(label=eph['label'], T0_BJD_TDB=eph['T0'], E_ref=eph['Eref'], P=eph['P'], cov=eph['cov']),
           n_transits=len(E), epoch_range=[int(E[0]), int(E[-1])],
           sigma_min_range=[float(sig.min()*1440), float(sig.max()*1440)],
           sigma_2027_min=[float(sig[(tb > 2461406.5) & (tb < 2461771.5)].min()*1440), float(sig[(tb > 2461406.5) & (tb < 2461771.5)].max()*1440)],
           sigma_2028_min=[float(sig[tb > 2461771.5].min()*1440), float(sig.max()*1440)],
           window_h_range=[float(2*half_win_h.min()), float(2*half_win_h.max())],
           light_time_check=check, ctx_comparison=cmp,
           ctx_ephemeris='ephemeris_v2: T0 2460840.9153 +- 0.0019, P 5.6914968 +- 0.0000145 (ctx/README)',
           expected_bjd_shift_from_ephemeris_change_min=ex,
           criteria='full = whole window (T14 1.70 h + 2x3 sigma + 2x1.5 h) at alt > 30 deg (airmass < 2) and Sun < -18 deg; '
                    'good_partial = transit +-3 sigma covered at Sun<-18/alt>30 with >=0.75 h baseline each side'),
      f'{OUT}/f01_transits.json')
