"""r01: independent transit predictions 2026-11-01..2028-12-31, BJD_TDB -> UTC at each site, and observability.
Own implementation (does not import the analyst's flib/f01):
 * ephemeris: clean/out/ephemeris_final.json 'combined' block (T0 at E_ref=-26, P, full covariance); also INV/REF blocks and
   phase-1 phot/out/timing.json for the systematic spread between ephemerides.
 * Roemer delay computed directly: d = (r_obs . n)/c with r_obs = barycentric Earth position (erfa epv00 via astropy
   get_body_barycentric) + GCRS observatory vector, n = unit vector to the star at the date (Gaia DR3 astrometry propagated
   with proper motion, parallax curvature term -(|r|^2-(r.n)^2)/(2 c d) included), Shapiro (Sun) included. Solve
   t_BJD = t_TDB(obs) + d(t) iteratively; TDB->UTC with astropy (leap seconds as currently tabulated).
 * Observability: 2-min grid; target altitude (with no refraction), Sun altitude, Moon from astropy get_body; 'full' window =
   T14 1.70 h + 2x3 sigma + 2x1.5 h baseline entirely at alt > 30 deg and Sun < -18 deg (also Sun < -12 deg).
Outputs: out/r01_transits.csv, out/r01_obs_all.csv, out/r01_summary.json"""
import os, sys, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, get_body_barycentric, get_body, get_sun, AltAz, SkyCoord, Distance
from astropy.utils import iers
iers.conf.auto_max_age = None
B='tic214311908/analysis'
O=f'{B}/referee3_followup/out'
eph_file = sys.argv[1] if len(sys.argv)>1 else f'{B}/clean/out/ephemeris_final.json'
E_all = json.load(open(eph_file))
def eph_block(name):
    b=E_all[name]; return dict(T0=b['T0_ref_BJD_TDB'], Eref=b['reference_epoch'], P=b['P_d'], C=np.array(b['cov_T0ref_P']))
eph = eph_block(E_all.get('adopted','combined'))
def pred(e, E):
    dE=np.asarray(E,float)-e['Eref']; t=e['T0']+e['P']*dE
    v=e['C'][0,0]+dE**2*e['C'][1,1]+2*dE*e['C'][0,1]; return t, np.sqrt(v)
# Gaia DR3 astrometry of the M3V (from data/gaia_cone25.csv, own query)
g=pd.read_csv(f'{B}/referee3_followup/data/gaia_cone25.csv'); g=g[g.source_id==5925373272850784768].iloc[0]
c_au_per_s = 299792.458/1.495978707e8
dist_au = 1000/g.parallax*206264.806
def unitvec(t_jyear):
    dt=t_jyear-2016.0
    dec=np.radians(g.dec+g.pmdec*dt/3.6e6); ra=np.radians(g.ra+g.pmra*dt/3.6e6/np.cos(np.radians(g.dec)))
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)*np.ones_like(ra)])
def delay_s(t_tdb, loc):
    """barycentric delay (s) for observer at loc at TDB time t (Time array)"""
    pe = get_body_barycentric('earth', t_tdb).xyz.to(u.au).value
    gc = loc.get_gcrs(t_tdb).cartesian.xyz.to(u.au).value
    r = pe + gc
    n = unitvec(t_tdb.jyear)
    rn = np.sum(r*n, axis=0); r2=np.sum(r*r,axis=0)
    roemer = rn/c_au_per_s - (r2-rn**2)/(2*dist_au)/c_au_per_s
    ps = get_body_barycentric('sun', t_tdb).xyz.to(u.au).value
    rs = r-ps; rsn=np.sum(rs*n,axis=0); rsm=np.sqrt(np.sum(rs*rs,axis=0))
    shapiro = 2*4.925490947e-6*np.log(1+rsn/rsm)
    return roemer + shapiro
def bjd2utc(tb, loc):
    t = Time(tb, format='jd', scale='tdb')
    for _ in range(4):
        d = delay_s(t, loc); t = Time(tb - d/86400., format='jd', scale='tdb')
    return t.utc, d
tstart, tend = Time('2026-11-01T00:00:00', scale='utc'), Time('2028-12-31T23:59:59', scale='utc')
Es = np.arange(-50, 400); tbs,_ = pred(eph, Es)
E = Es[(tbs > tstart.tdb.jd) & (tbs < tend.tdb.jd)]
tb, sig = pred(eph, E)
# alternative ephemerides
alts={}
for nm in ('inv','ref'):
    alts[nm]=pred(eph_block(nm),E)[0]
ph1=json.load(open(f'{B}/phot/out/timing.json'))
print('phase-1 keys', list(ph1.keys())[:10])
b=ph1['all_err_cal']; alts['phase1']=b['T0']+2457000+b['P']*E
geo=EarthLocation.from_geocentric(0*u.m,0*u.m,0*u.m)
ug, dg = bjd2utc(tb, geo)
T = pd.DataFrame(dict(epoch=E, bjd_tdb=tb, sig_min=sig*1440, utc_geo=[x.isot[:19] for x in ug], delay_geo_s=dg,
                      d_inv_min=(alts['inv']-tb)*1440, d_ref_min=(alts['ref']-tb)*1440, d_phase1_min=(alts['phase1']-tb)*1440))
# compare with analyst
A = pd.read_csv(f'{B}/followup/out/transits_all.csv')
M = T.merge(A, on='epoch', suffixes=('','_an'))
dutc = (Time(list(M.utc_geo)).jd - Time(list(M.utc_mid_geocentre)).jd)*86400
cmp = dict(n=len(M), dBJD_s_maxabs=float(np.max(np.abs(M.bjd_tdb-M.bjd_tdb_an))*86400), dUTC_geo_s_range=[float(dutc.min()), float(dutc.max())],
           ddelay_s_maxabs=float(np.max(np.abs(M.delay_geo_s-M.bary_corr_s))), dsig_maxabs_min=float(np.max(np.abs(M.sig_min-M.sigma_min))))
print('vs analyst', cmp)
print('ephemeris spread (min) inv/ref/phase1 range:', [ (round(T[c].min(),2), round(T[c].max(),2)) for c in ('d_inv_min','d_ref_min','d_phase1_min')])
print('sigma range 2027', T[(T.bjd_tdb<2461771.5)].sig_min.agg(['min','max']).values, '2028', T[T.bjd_tdb>2461771.5].sig_min.agg(['min','max']).values)
# ctx comparison
C = pd.read_csv(f'{B}/ctx/out/followup_transits.csv')
print(C.columns.tolist()); 
SITES = {'CTIO':(-30.1690,-70.8063,2207.),'Cerro_Pachon':(-30.2408,-70.7367,2722.),'Las_Campanas':(-29.0146,-70.6926,2380.),
         'La_Silla':(-29.2563,-70.7380,2400.),'Paranal':(-24.6272,-70.4042,2635.),'SAAO':(-32.3794,20.8107,1798.),
         'Siding_Spring':(-31.2733,149.0614,1165.),'Dome_C':(-75.0997,123.3325,3233.)}
T14=1.70
coord = SkyCoord(ra=g.ra*u.deg, dec=g.dec*u.deg, distance=Distance(parallax=g.parallax*u.mas), pm_ra_cosdec=g.pmra*u.mas/u.yr,
                 pm_dec=g.pmdec*u.mas/u.yr, radial_velocity=g.radial_velocity*u.km/u.s, obstime=Time(2016.0,format='jyear',scale='tcb')).apply_space_motion(new_obstime=Time(2027.5,format='jyear',scale='tcb'))
rows=[]
for s,(lat,lon,h) in SITES.items():
    loc=EarthLocation.from_geodetic(lon*u.deg,lat*u.deg,h*u.m)
    uu, dd = bjd2utc(tb, loc)
    for i in range(len(E)):
        hw = T14/2 + 3*sig[i]*24 + 1.5
        off = np.arange(-hw, hw+1e-9, 2/60.)
        tt = uu[i] + off*u.hour
        fr = AltAz(obstime=tt, location=loc, pressure=0*u.hPa)
        alt = coord.transform_to(fr).alt.deg
        sun = get_sun(tt).transform_to(fr).alt.deg
        tm = uu[i]; fm = AltAz(obstime=tm, location=loc, pressure=0*u.hPa)
        moon = get_body('moon', tm, loc); sunm = get_body('sun', tm, loc)
        # illuminated fraction from the Sun-Moon-observer phase angle (exact geometry)
        ms = (moon.cartesian.xyz - sunm.cartesian.xyz) if False else None
        mg = moon.transform_to(fm); 
        elong = sunm.separation(moon).rad
        dm = moon.distance.to(u.au).value; ds = sunm.distance.to(u.au).value
        pa = np.arctan2(ds*np.sin(elong), dm - ds*np.cos(elong))
        illum = (1+np.cos(pa))/2
        rows.append(dict(site=s, epoch=int(E[i]), utc_mid=tm.isot[:16], delay_s=dd[i], win_h=2*hw,
                         full18=bool(((alt>30)&(sun<-18)).all()), full12=bool(((alt>30)&(sun<-12)).all()),
                         frac18=float(((alt>30)&(sun<-18)).mean()), alt_min=float(alt.min()), sun_max=float(sun.max()),
                         moon_illum=float(illum), moon_sep=float(moon.separation(coord.transform_to(moon.frame) if False else coord).deg), moon_alt=float(mg.alt.deg)))
    print(s, 'done', flush=True)
R=pd.DataFrame(rows); R.to_csv(f'{O}/r01_obs_all.csv', index=False); T.to_csv(f'{O}/r01_transits.csv', index=False)
# compare with analyst per-site
AS = pd.read_csv(f'{B}/followup/out/transits_obs_all_sites.csv')
mp = {'CTIO':'CTIO (LCOGT lsc, SMARTS)','Cerro_Pachon':'Cerro Pachon (SOAR, Gemini-S)','Las_Campanas':'Las Campanas (LCO)','La_Silla':'La Silla (TRAPPIST-S, ExTrA, 3.6m)',
      'Paranal':'Paranal (VLT, SPECULOOS-S)','SAAO':'SAAO Sutherland (LCOGT cpt, SALT)','Siding_Spring':'Siding Spring (LCOGT coj, FTS/MuSCAT4, AAT)','Dome_C':'Dome C (ASTEP+)'}
summ={}
for s in SITES:
    r=R[R.site==s]; a=AS[AS.site==mp[s]].merge(r, on='epoch', suffixes=('_an',''))
    dd = (Time([x.replace(' ','T') for x in a.utc_mid_an]).jd - Time(list(a.utc_mid)).jd)*1440
    f = r[r.full18]
    summ[s]=dict(n_full18_2027=int((f.utc_mid<'2028').sum()), n_full18_2028=int((f.utc_mid>='2028').sum()),
                 n_full12_2027=int(((r.utc_mid<'2028')&r.full12).sum()), n_full12_2028=int(((r.utc_mid>='2028')&r.full12).sum()),
                 full18_list=[f"E{x.epoch} {x.utc_mid} moon {x.moon_illum:.2f}/{x.moon_sep:.0f}" for x in f.itertuples()],
                 disagree_full18_epochs=[int(x) for x in a.epoch[a.full_window_sun18!=a.full18]],
                 disagree_full12_epochs=[int(x) for x in a.epoch[a.full_window_sun12!=a.full12]],
                 dutc_min_range=[float(dd.min()), float(dd.max())], ddelay_s_max=float(np.max(np.abs(a.bary_corr_s-a.delay_s))),
                 dmoon_illum_max=float(np.max(np.abs(a.moon_illum_an-a.moon_illum))), dmoon_sep_max=float(np.max(np.abs(a.moon_sep_deg-a.moon_sep))),
                 first_any_dark_usable=r[r.frac18>0].utc_mid.min())
    print(s, {k:v for k,v in summ[s].items() if k!='full18_list'}); print('   ', summ[s]['full18_list'])
json.dump(dict(eph=dict(file=eph_file, T0=eph['T0'], Eref=eph['Eref'], P=eph['P']), vs_analyst=cmp, sites=summ,
               ephemeris_spread_min={c:[float(T[c].min()), float(T[c].max())] for c in ('d_inv_min','d_ref_min','d_phase1_min')},
               delay_geo_range_s=[float(T.delay_geo_s.min()), float(T.delay_geo_s.max())]), open(f'{O}/r01_summary.json','w'), indent=1)
