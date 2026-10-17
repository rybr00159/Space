"""r01b: (i) distinct epochs with a full window from at least one southern site (network view), per year;
(ii) earliest epochs with the full transit (T14 + 2x3 sigma) observable at alt > 30, Sun < -18 deg (partial-baseline events),
own 2-min grid; (iii) true UTC at Earth of the 2024-10-29 Gaia crossing (ctx 'utc' column is BJD_TDB on the UTC scale).
Output: out/r01b_network.json"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz, get_sun, SkyCoord
B='tic214311908/analysis'; R=f'{B}/referee3_followup'
O=pd.read_csv(f'{R}/out/r01_obs_all.csv'); T=pd.read_csv(f'{R}/out/r01_transits.csv')
res={}
for yr,sel in (('2027',O.utc_mid<'2028'),('2028',O.utc_mid>='2028')):
    s=O[sel&O.full18]
    res[f'distinct_full18_epochs_{yr}_excl_DomeC']=sorted(set(s[s.site!='Dome_C'].epoch.astype(int)))
    res[f'distinct_full18_epochs_{yr}_incl_DomeC']=sorted(set(s.epoch.astype(int)))
    s12=O[sel&O.full12]; res[f'distinct_full12_epochs_{yr}_incl_DomeC']=sorted(set(s12.epoch.astype(int)))
print({k:(len(v),v) for k,v in res.items()})
# earliest transit-only coverage
g=pd.read_csv(f'{R}/data/gaia_cone25.csv'); g=g[g.source_id==5925373272850784768].iloc[0]
c=SkyCoord(ra=g.ra*u.deg,dec=g.dec*u.deg)
SITES={'CTIO':(-30.1690,-70.8063,2207.),'SAAO':(-32.3794,20.8107,1798.),'Siding_Spring':(-31.2733,149.0614,1165.)}
early=[]
for s,(lat,lon,h) in SITES.items():
    loc=EarthLocation.from_geodetic(lon*u.deg,lat*u.deg,h*u.m)
    for r in O[(O.site==s)&(O.utc_mid<'2027-04-20')].itertuples():
        sig=T[T.epoch==r.epoch].sig_min.iloc[0]/1440
        hw=0.85/24+3*sig
        tt=Time(r.utc_mid)+np.arange(-hw,hw+1e-9,2/1440)*u.day
        fr=AltAz(obstime=tt,location=loc,pressure=0*u.hPa)
        ok=(c.transform_to(fr).alt.deg>30)&(get_sun(tt).transform_to(fr).alt.deg<-18)
        # baseline available before/after within 1.5 h
        tb=Time(r.utc_mid)+np.concatenate([np.arange(-hw-1.5/24,-hw,2/1440),np.arange(hw,hw+1.5/24,2/1440)])*u.day
        frb=AltAz(obstime=tb,location=loc,pressure=0*u.hPa)
        okb=(c.transform_to(frb).alt.deg>30)&(get_sun(tb).transform_to(frb).alt.deg<-18)
        if ok.all(): early.append(dict(site=s,epoch=int(r.epoch),utc=r.utc_mid,baseline_h=round(float(okb.sum()*2/60),2)))
print(early)
res['transit_fully_covered_before_2027-04-20']=early
# Gaia 2024-10-29 crossing: true UTC at the geocentre (Gaia at L2 differs by <~5 s)
bjd=2460613.265113589
t=Time(bjd,format='jd',scale='tdb')
for _ in range(3):
    d=t.light_travel_time(c,'barycentric',location=EarthLocation.from_geocentric(0*u.m,0*u.m,0*u.m)); t=Time(bjd,format='jd',scale='tdb')-d
res['gaia_2024-10-29_true_UTC_geocentric']=t.utc.isot; res['ctx_label']='2024-10-29 18:20:36.6 (BJD_TDB on UTC scale, no light-time)'
print(res['gaia_2024-10-29_true_UTC_geocentric'])
json.dump(res,open(f'{R}/out/r01b_network.json','w'),indent=1)
