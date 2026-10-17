import numpy as np, pandas as pd
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, get_body, get_sun
import astropy.units as u
W='tic214311908'
c=SkyCoord('17h23m52.92s','-51d19m27.1s')  # 2MASS position of the M3V
geo=EarthLocation.from_geocentric(0*u.m,0*u.m,0*u.m)
def bjd_to_utc(bjd):
    # iterate: guess UTC = bjd; compute ltt; refine
    t=Time(bjd,format='jd',scale='tdb',location=geo)
    for i in range(4):
        ltt=t.light_travel_time(c,kind='barycentric')
        t=Time(bjd,format='jd',scale='tdb',location=geo)-ltt
    return t.utc, ltt.to(u.s).value
df=pd.read_csv(W+'/predicted_transits_2026_2028.csv')
T0,P=2460840.915,5.691506
for E in [89,118,126,127,137,149,227]:
    r=df[df.epoch==E].iloc[0]
    bjd=T0+E*P
    t,ltt=bjd_to_utc(bjd)
    tt=Time(r.utc_mid,format='isot',scale='utc')
    sun=get_sun(t); moon=get_body('moon',t)
    elong=sun.separation(moon).deg
    illum=(1-np.cos(np.radians(elong)))/2
    sepm=moon.separation(c).deg
    print(E, 'bjd file %.6f calc %.6f'%(r.bjd_tdb_mid,bjd), 'utc file',r.utc_mid,'calc',t.isot[:19], 'diff_s %.1f'%((t-tt).to(u.s).value),
          'ltt_s %.1f file %.1f'%(ltt,r.bary_corr_s),'moon file %.2f/%.1f calc %.2f/%.1f'%(r.moon_illum,r.moon_sep_deg,illum,sepm))
old=pd.read_csv('/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/predicted_transits_2027.csv')
new=pd.read_csv(W+'/predicted_transits_2027.csv')
# old-file error: own bjd -> correct utc vs its utc_mid
errs=[]
for _,r in old.iterrows():
    t,_=bjd_to_utc(r.bjd_tdb)
    errs.append((t-Time(r.utc_mid.replace(' ','T'),format='isot',scale='utc')).to(u.min).value)
errs=np.array(errs); print('old file: correct UTC minus listed utc_mid (min): %.2f to %.2f'%(errs.min(),errs.max()), 'epochs', old.epoch.min(), old.epoch.max(), 'sigma', old.sigma_min.min(), old.sigma_min.max())
e2=[]
for _,r in new.iterrows():
    t,_=bjd_to_utc(r.bjd_tdb)
    e2.append((t-Time(r.utc_mid.replace(' ','T'),format='isot',scale='utc')).to(u.s).value)
e2=np.array(e2); print('new 2027 file: |calc-listed| max %.1f s'%np.abs(e2).max(), 'sigma', new.sigma_min.min(), new.sigma_min.max(), 'epochs',new.epoch.min(),new.epoch.max())
# new bjd vs adopted ephemeris
print('new bjd max diff from T0+EP (s):', np.abs(new.bjd_tdb-(T0+new.epoch*P)).max()*86400)
# old vs new utc difference
d=[(Time(b.replace(' ','T'))-Time(a.replace(' ','T'))).to(u.min).value for a,b in zip(old.utc_mid,new.utc_mid)]
print('new minus old utc_mid (min): %.1f to %.1f'%(min(d),max(d)))
