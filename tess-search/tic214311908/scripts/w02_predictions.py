"""w02: transit predictions 2026-11-01 .. 2028-12-31 from the ADOPTED ephemeris (writeup/ephemeris_adopted.json), with
observability per site. Adapted from followup/scripts/f01_transits.py (same window definition and site coordinates).

Conventions
* bjd_tdb_mid: target-barycentred BJD_TDB of mid-transit. sigma from the full (T0, P) covariance of the adopted ephemeris.
* utc_mid: GEOCENTRIC UTC, i.e. the UTC at which mid-transit light reaches the Earth's centre. Solved iteratively from
  t_BJD = t_geo(TDB) + light-travel time to the solar-system barycentre (astropy, Roemer delay, target position
  propagated with the Gaia DR3 proper motion), then TDB -> UTC (no new leap second assumed before 2029). Topocentric
  times differ by < 0.03 s at any site.
* utc_ingress / utc_egress = utc_mid -/+ T14/2 with T14 = 1.70 h (no timing margin).
* Full window at a site: the whole interval T14 + 2 x 3 sigma + 2 x 1.5 h baseline has the target at airmass < 2
  (altitude > 30 deg) and the Sun below -18 deg; 2-min grid. Also given: the same with Sun < -12 deg, whether the
  transit +- 3 sigma alone is covered, and 'good partial' (transit +- 3 sigma covered with >= 0.75 h baseline each side).
* Moon illumination from the exact Sun-Moon phase angle (geocentric) at mid-transit.
* 'CTIO/SOAR/Gemini-S' is evaluated at both Cerro Tololo and Cerro Pachon; a window counts as full only if it is full at
  both (they never differ in practice).
Outputs (writeup/): predicted_transits_2026_2028.csv, observability/obs_<site>.csv (+ observability/README.md), predicted_transits_2027.csv
(vetting-file format, epochs 108-149), scripts/w02_checks.json."""
import os, sys, json
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(k, '2')
import numpy as np, pandas as pd, warnings
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body, SkyCoord, Distance
from astropy.utils import iers
iers.conf.auto_max_age = None
iers.conf.iers_degraded_accuracy = 'warn'
warnings.filterwarnings('ignore')

W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the write-up folder (writeup/ or tic214311908/)
B = os.path.dirname(W) if os.path.isdir(os.path.join(os.path.dirname(W), 'clean')) else os.path.join(W, 'analysis')
REPO_2027 = '/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/predicted_transits_2027.csv'

ej = json.load(open(f'{W}/ephemeris_adopted.json'))
blk = ej[ej['adopted']]
T0, Eref, P, C = blk['T0_ref_BJD_TDB'], blk['reference_epoch'], blk['P_d'], np.array(blk['cov_T0ref_P'])
T14_H = blk['T14_h']

def predict(E):
    dE = np.asarray(E, float) - Eref
    return T0 + P*dE, np.sqrt(C[0, 0] + dE**2*C[1, 1] + 2*dE*C[0, 1])

# target: Gaia DR3 5925373272850784768 (as followup/scripts/flib.py)
TGT = dict(ra=260.969674, dec=-51.323821, plx=15.497306, pmra=-113.671152, pmdec=85.54355, rv=16.28)
c16 = SkyCoord(ra=TGT['ra']*u.deg, dec=TGT['dec']*u.deg, distance=Distance(parallax=TGT['plx']*u.mas),
               pm_ra_cosdec=TGT['pmra']*u.mas/u.yr, pm_dec=TGT['pmdec']*u.mas/u.yr, radial_velocity=TGT['rv']*u.km/u.s,
               obstime=Time(2016.0, format='jyear', scale='tcb'))
coord = c16.apply_space_motion(new_obstime=Time(2027.5, format='jyear', scale='tcb'))
GEO = EarthLocation.from_geocentric(0*u.m, 0*u.m, 0*u.m)

def bjd_to_utc(tbjd, loc=GEO):
    tb = Time(tbjd, format='jd', scale='tdb', location=loc)
    t = tb
    for _ in range(4):
        ltt = t.light_travel_time(coord, kind='barycentric', location=loc)
        t = tb - ltt
    return t.utc, ltt.to(u.s).value

def iso_s(t):          # ISO to the nearest second
    return [x.isot[:19] for x in (t + 0.5*u.s)] if t.isscalar is False else (t + 0.5*u.s).isot[:19]

def moon_illum_geo(t):
    s = get_sun(t); m = get_body('moon', t)
    elong = s.separation(m).rad
    rs, rm = s.distance.to(u.km).value, m.distance.to(u.km).value
    i = np.arctan2(rs*np.sin(elong), rm - rs*np.cos(elong))      # phase angle
    return (1 + np.cos(i))/2, m

# ---------------- epochs ----------------
tstart, tend = Time('2026-11-01T00:00:00', scale='utc'), Time('2028-12-31T23:59:59', scale='utc')
Es = np.arange(-50, 400)
tbs, _ = predict(Es)
E = Es[(tbs > tstart.tdb.jd + 0.01) & (tbs < tend.tdb.jd - 0.01)]
tb, sig = predict(E)
tutc, ltt = bjd_to_utc(tb)
illum, moon = moon_illum_geo(tutc)
moonsep = moon.separation(coord).deg

# ---------------- observability ----------------
SITES = {   # name: (lat, lon E, height m), as followup/scripts/flib.py
    'CTIO': (-30.1690, -70.8063, 2207.), 'Cerro_Pachon': (-30.2408, -70.7367, 2722.),
    'Las_Campanas': (-29.0146, -70.6926, 2380.), 'La_Silla': (-29.2563, -70.7380, 2400.),
    'Paranal': (-24.6272, -70.4042, 2635.), 'SAAO': (-32.3794, 20.8107, 1798.),
    'Siding_Spring': (-31.2733, 149.0614, 1165.), 'Dome_C': (-75.0997, 123.3325, 3233.)}
GROUPS = {   # output file -> (label, member sites)
    'CTIO_SOAR_GeminiS': ('CTIO / Cerro Pachon (SOAR, Gemini-S), Chile', ['CTIO', 'Cerro_Pachon']),
    'LCO': ('Las Campanas Observatory, Chile', ['Las_Campanas']),
    'La_Silla': ('La Silla (TRAPPIST-S, ExTrA, ESO 3.6 m), Chile', ['La_Silla']),
    'Paranal': ('Paranal (VLT, SPECULOOS-South), Chile', ['Paranal']),
    'SAAO': ('SAAO Sutherland (LCOGT, SALT), South Africa', ['SAAO']),
    'Siding_Spring': ('Siding Spring (FTS/MuSCAT4, LCOGT, AAT), Australia', ['Siding_Spring']),
    'Dome_C': ('Dome C (ASTEP+), Antarctica', ['Dome_C'])}
half_win_h = T14_H/2 + 1.5 + 3*sig*24
rec = {}
for name, (lat, lon, h) in SITES.items():
    loc = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg, h*u.m)
    rows = []
    for i in range(len(E)):
        hw = half_win_h[i]
        n = int(np.ceil(2*hw*60/2)) + 1                   # 2-min grid
        off = np.linspace(-hw, hw, n)
        tt = tutc[i] + off*u.hour
        fr = AltAz(obstime=tt, location=loc)
        alt = coord.transform_to(fr).alt.deg
        sun = get_sun(tt).transform_to(fr).alt.deg
        g18 = (alt > 30) & (sun < -18); g12 = (alt > 30) & (sun < -12)
        half_tr = T14_H/2 + 3*sig[i]*24
        intr = np.abs(off) <= half_tr
        dt = 2*hw/(n - 1)
        pre = ((off < -half_tr) & g18).sum()*dt; post = ((off > half_tr) & g18).sum()*dt
        k = n//2
        mloc = get_body('moon', tutc[i], loc).transform_to(AltAz(obstime=tutc[i], location=loc))
        rows.append(dict(full18=bool(g18.all()), full12=bool(g12.all()), transit18=bool(g18[intr].all()),
                         pre_h=round(float(pre), 2), post_h=round(float(post), 2), frac18=round(float(g18.mean()), 3),
                         alt_mid=round(float(alt[k]), 1), alt_min=round(float(alt.min()), 1),
                         airmass_max=round(float(1/np.sin(np.radians(max(alt.min(), 5)))), 2),
                         sun_max=round(float(sun.max()), 1), moon_alt=round(float(mloc.alt.deg), 1),
                         win_start=(tutc[i] - hw*u.hour), win_end=(tutc[i] + hw*u.hour), window_h=round(float(2*hw), 2)))
    rec[name] = rows
    print(name, 'done', flush=True)

def iso_min(t):     # nearest minute
    return (t + 30*u.s).isot[:16].replace('T', ' ')

# CTIO vs Cerro Pachon consistency
diff_cp = [int(E[i]) for i in range(len(E)) if rec['CTIO'][i]['full18'] != rec['Cerro_Pachon'][i]['full18']]
full_sites, full_sites12, nfull = [], [], []
for i in range(len(E)):
    fs, fs12 = [], []
    for g, (lab, mem) in GROUPS.items():
        if all(rec[m][i]['full18'] for m in mem): fs.append(g)
        if all(rec[m][i]['full12'] for m in mem): fs12.append(g)
    full_sites.append(';'.join(fs)); full_sites12.append(';'.join(fs12)); nfull.append(len(fs))

utc_mid = [x.isot[:19] for x in (tutc + 0.5*u.s)]
ing = tutc - T14_H/2*u.hour; egr = tutc + T14_H/2*u.hour
main = pd.DataFrame(dict(
    epoch=E, bjd_tdb_mid=np.round(tb, 6), sigma_min=np.round(sig*1440, 2), sigma3_min=np.round(3*sig*1440, 2),
    utc_mid=utc_mid, utc_ingress=[x.isot[:19] for x in (ing + 0.5*u.s)], utc_egress=[x.isot[:19] for x in (egr + 0.5*u.s)],
    bary_corr_s=np.round(ltt, 2), moon_illum=np.round(illum, 2), moon_sep_deg=np.round(moonsep, 1),
    n_full_window_sites=nfull, full_window_sites=full_sites, full_window_sites_sun12=full_sites12))
main.to_csv(f'{W}/predicted_transits_2026_2028.csv', index=False)

site_lines = []
for g, (lab, mem) in GROUPS.items():
    out = []
    for i in range(len(E)):
        r = [rec[m][i] for m in mem]
        a = r[0]
        if max(x['frac18'] for x in r) <= 0: continue
        out.append(dict(epoch=int(E[i]), bjd_tdb_mid=round(float(tb[i]), 6), utc_mid=utc_mid[i],
                        utc_window_start=iso_min(a['win_start']), utc_window_end=iso_min(a['win_end']),
                        utc_ingress=iso_min(ing[i]), utc_egress=iso_min(egr[i]),
                        sigma_min=round(float(sig[i]*1440), 2), sigma3_min=round(float(3*sig[i]*1440), 2), window_h=a['window_h'],
                        full_window=all(x['full18'] for x in r), full_window_sun12=all(x['full12'] for x in r),
                        transit_covered=all(x['transit18'] for x in r),
                        good_partial=(all(x['transit18'] for x in r) and min(x['pre_h'] for x in r) >= 0.75 and
                                      min(x['post_h'] for x in r) >= 0.75 and not all(x['full18'] for x in r)),
                        baseline_pre_h=min(x['pre_h'] for x in r), baseline_post_h=min(x['post_h'] for x in r),
                        frac_window_usable=min(x['frac18'] for x in r), alt_mid=a['alt_mid'],
                        alt_min_window=min(x['alt_min'] for x in r), airmass_max_window=max(x['airmass_max'] for x in r),
                        sun_max_window=max(x['sun_max'] for x in r), moon_illum=round(float(illum[i]), 2),
                        moon_sep_deg=round(float(moonsep[i]), 1), moon_alt_mid=a['moon_alt']))
    df = pd.DataFrame(out)
    df.to_csv(f'{W}/observability/obs_{g}.csv', index=False)      # plain CSV, header first; described in observability/README.md
    site_lines.append(f'| `obs_{g}.csv` | {lab} | ' + '; '.join(f"{m.replace('_', ' ')}: lat {SITES[m][0]}, lon {SITES[m][1]}, {SITES[m][2]:.0f} m" for m in mem) + ' |')

OBS_README = '''# Observability tables, TIC 214311908 transits 2026-11 to 2028-12

One plain CSV per site (header in the first line). Each lists the epochs with any usable time at the site (Sun below
-18 deg and the target at airmass < 2). All times come from the adopted ephemeris (`../ephemeris_adopted.json`) and were
computed by `../scripts/w02_predictions.py`. A window counts as full when T14 (1.70 h) + 2 x 3 sigma of timing margin +
2 x 1.5 h of baseline is entirely usable. AI-run analysis; not reviewed by a professional astronomer.

| File | Site | Coordinates (deg, deg, m) |
|---|---|---|
''' + '\n'.join(site_lines) + '''

Columns:
* `epoch`: transit number from the adopted T0; `bjd_tdb_mid`: mid-transit, target-barycentred BJD_TDB.
* `utc_mid`: geocentric UTC of mid-transit, to the second. The other UTC columns are rounded to the minute:
  `utc_window_start`, `utc_window_end` (full window), `utc_ingress`, `utc_egress` (T14 only, no margin).
* `sigma_min`, `sigma3_min`: 1 sigma and 3 sigma timing uncertainty (min); `window_h`: window length (h).
* `full_window`: the whole window is usable with the Sun below -18 deg (all member sites); `full_window_sun12`: the
  same with the Sun below -12 deg; `transit_covered`: transit +- 3 sigma usable; `good_partial`: transit +- 3 sigma
  usable with at least 0.75 h of baseline on each side, but not a full window.
* `baseline_pre_h`, `baseline_post_h`: usable baseline before and after (h); `frac_window_usable`: usable fraction
  of the window.
* `alt_mid`, `alt_min_window` (deg), `airmass_max_window`, `sun_max_window` (deg): target and Sun during the window.
* `moon_illum` (0-1), `moon_sep_deg`, `moon_alt_mid` (deg): the Moon at mid-transit.
'''
open(f'{W}/observability/README.md', 'w').write(OBS_README)

# ---------------- vetting-format 2027 file (epochs 108-149) ----------------
old = pd.read_csv(REPO_2027)
E27 = old.epoch.values
tb27, s27 = predict(E27)
u27, _ = bjd_to_utc(tb27)
new27 = pd.DataFrame(dict(epoch=E27, bjd_tdb=[round(float(x), 4) for x in tb27],
                          utc_mid=[(x + 30*u.s).isot[:16].replace('T', ' ') for x in u27],
                          sigma_min=[round(float(x*1440), 1) for x in s27]))
with open(f'{W}/predicted_transits_2027.csv', 'w') as f:
    f.write('epoch,bjd_tdb,utc_mid,sigma_min\n')
    for r in new27.itertuples():
        f.write(f"{r.epoch},{r.bjd_tdb!r},{r.utc_mid},{r.sigma_min}\n")

# ---------------- checks ----------------
chk = {}
# (a) the old repo file: its utc_mid against barycentric-correct UTC of its own bjd_tdb, and against naive TDB->UTC
uo, _ = bjd_to_utc(old.bjd_tdb.values)
old_utc = Time(list(old.utc_mid.str.replace(' ', 'T')), scale='utc')
naive = Time(old.bjd_tdb.values, format='jd', scale='tdb').utc
chk['old_repo_file'] = dict(
    file=REPO_2027,
    old_utc_minus_correct_min=[round(float(x), 2) for x in ((old_utc.jd - uo.jd)*1440)[[0, -1]]],
    old_utc_minus_correct_min_range=[round(float(((old_utc.jd - uo.jd)*1440).min()), 2), round(float(((old_utc.jd - uo.jd)*1440).max()), 2)],
    old_utc_minus_naive_tdb_to_utc_min_range=[round(float(((old_utc.jd - naive.jd)*1440).min()), 2), round(float(((old_utc.jd - naive.jd)*1440).max()), 2)],
    old_bjd_minus_adopted_min_range=[round(float(((old.bjd_tdb.values - tb27)*1440).min()), 2), round(float(((old.bjd_tdb.values - tb27)*1440).max()), 2)],
    old_utc_minus_new_utc_min_range=[round(float(((old_utc.jd - Time(list(new27.utc_mid.str.replace(' ', 'T')), scale='utc').jd)*1440).min()), 1),
                                     round(float(((old_utc.jd - Time(list(new27.utc_mid.str.replace(' ', 'T')), scale='utc').jd)*1440).max()), 1)],
    old_sigma_min=[float(old.sigma_min.min()), float(old.sigma_min.max())], new_sigma_min=[float(new27.sigma_min.min()), float(new27.sigma_min.max())])
# (b) barycentric conversion against followup/ and referee3_followup/ (same BJD in, compare UTC out)
fu = pd.read_csv(f'{B}/followup/out/transits_all.csv')
ufu, _ = bjd_to_utc(fu.bjd_tdb.values)
d_fu = (Time(list(fu.utc_mid_geocentre), scale='utc').jd - ufu.jd)*86400
rf = pd.read_csv(f'{B}/referee3_followup/out/r01_transits.csv')
urf, lrf = bjd_to_utc(rf.bjd_tdb.values)
d_rf = (Time(list(rf.utc_geo), scale='utc').jd - urf.jd)*86400
chk['bary_conversion_check'] = dict(
    followup_utc_minus_this_s=[round(float(d_fu.min()), 2), round(float(d_fu.max()), 2)],
    referee3_followup_utc_minus_this_s=[round(float(d_rf.min()), 2), round(float(d_rf.max()), 2)],
    referee3_delay_minus_this_s=[round(float((rf.delay_geo_s.values - lrf).min()), 3), round(float((rf.delay_geo_s.values - lrf).max()), 3)],
    note='both reference files give UTC truncated to the second, so differences of 0-1 s are rounding; referee3 includes Shapiro and parallax-curvature terms')
# (c) differences from the clean ephemeris used by followup/ for the same epochs
m = fu.merge(main[['epoch', 'bjd_tdb_mid', 'sigma_min']], on='epoch', suffixes=('_clean', ''))
chk['adopted_minus_clean_bjd_min'] = [round(float(((m.bjd_tdb_mid - m.bjd_tdb)*1440).min()), 2), round(float(((m.bjd_tdb_mid - m.bjd_tdb)*1440).max()), 2)]
chk['sigma_min_range'] = dict(all=[float(sig.min()*1440), float(sig.max()*1440)],
                              y2027=[float(sig[(tb > Time('2027-01-01').jd) & (tb < Time('2028-01-01').jd)].min()*1440), float(sig[(tb > Time('2027-01-01').jd) & (tb < Time('2028-01-01').jd)].max()*1440)],
                              y2028=[float(sig[tb > Time('2028-01-01').jd].min()*1440), float(sig.max()*1440)])
chk['n_transits'] = int(len(E)); chk['epoch_range'] = [int(E[0]), int(E[-1])]
chk['ctio_vs_cerro_pachon_full_differs_epochs'] = diff_cp
# network counts
def yr(i): return utc_mid[i][:4]
for y in ('2026', '2027', '2028'):
    idx = [i for i in range(len(E)) if yr(i) == y]
    chk[f'full_epochs_{y}_excl_DomeC'] = [int(E[i]) for i in idx if any(s != 'Dome_C' for s in full_sites[i].split(';') if s)]
    chk[f'full_epochs_{y}_incl_DomeC'] = [int(E[i]) for i in idx if full_sites[i]]
    chk[f'full_epochs_{y}_sun12_incl_DomeC'] = [int(E[i]) for i in idx if full_sites12[i]]
per_site = {}
for g in GROUPS:
    d = pd.read_csv(f'{W}/observability/obs_{g}.csv', comment='#')
    per_site[g] = {y: [f"E{r.epoch} {(Time(r.utc_mid, scale='utc') + 30*u.s).isot[:16]} moon {r.moon_illum:.2f}/{r.moon_sep_deg:.0f}deg" for r in d[d.full_window & d.utc_mid.str.startswith(y)].itertuples()] for y in ('2027', '2028')}
chk['full_windows_per_site'] = per_site
first = main[main.n_full_window_sites > 0].iloc[0]
chk['first_full_window'] = dict(epoch=int(first.epoch), utc_mid=first.utc_mid, sites=first.full_window_sites)
json.dump(chk, open(f'{W}/scripts/w02_checks.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in chk.items() if k != 'full_windows_per_site'}, indent=1))
print(json.dumps(per_site, indent=1))
