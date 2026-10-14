"""Target-specific barycentric correction for TESS FFI-based products (TESScut, QLP).
TESScut/QLP times carry the FFI's barycentric correction (computed for a reference direction, not
for our star).  We recompute r_TESS . n_target / c from the JPL Horizons TESS (-95) SSB ephemeris,
validate against the per-target SPOC TIMECORR in S1/S28/S68, and save corrections per sector."""
import numpy as np, glob, json, signal
from astropy.io import fits
from astroquery.jplhorizons import Horizons
from astropy.time import Time
OUT='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/data'
L='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut/'
R='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/ctoi_vet/recent/data/'
RA,DEC=np.deg2rad(330.887948),np.deg2rad(-42.835511)
n=np.array([np.cos(DEC)*np.cos(RA),np.cos(DEC)*np.sin(RA),np.sin(DEC)])
AU_S=149597870.700/299792.458
def handler(s,f): raise TimeoutError
signal.signal(signal.SIGALRM, handler)
def romer(tsc):  # tsc: spacecraft TDB in BTJD
    j0=tsc.min()+2457000-0.1; j1=tsc.max()+2457000+0.1
    signal.alarm(200)
    v=Horizons(id='-95',location='@0',epochs={'start':Time(j0,format='jd',scale='tdb').iso,'stop':Time(j1,format='jd',scale='tdb').iso,'step':'20m'}).vectors(refplane='earth')
    signal.alarm(0)
    jd=np.array(v['datetime_jd'])-2457000; X=np.vstack([np.array(v[c]) for c in 'xyz']).T
    proj=X@n
    return np.interp(tsc,jd,proj)*AU_S
res={}
for sec in [1,28,68,95,105,106]:
    tc=fits.open(glob.glob(L+f's{sec}/*.fits')[0])[1].data
    t=np.array(tc['TIME']); corr=np.array(tc['TIMECORR'])*86400; ok=np.isfinite(t)&np.isfinite(corr)
    tsc=t[ok]-corr[ok]/86400
    mine=romer(tsc)
    d=dict(sector=sec,tesscut_minus_true_s=float(np.median(corr[ok]-mine)),tesscut_minus_true_range=[float(np.min(corr[ok]-mine)),float(np.max(corr[ok]-mine))])
    np.save(f'{OUT}/bary_tc_s{sec}.npy',np.vstack([t[ok],mine-corr[ok]]))  # add (mine-corr)/86400 to tesscut TIME
    sp=glob.glob(R+f'mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')
    if sp:
        assert int(fits.getheader(sp[0])['TICID'])==71316387
        s=fits.open(sp[0])[1].data; ts=np.array(s['TIME']); cs=np.array(s['TIMECORR'])*86400; o=np.isfinite(ts)&np.isfinite(cs)
        m2=romer(ts[o]-cs[o]/86400)
        d['spoc_minus_mine_s']=[float(np.median(cs[o]-m2)),float(np.std(cs[o]-m2))]
    # QLP comparison: QLP TIME vs our corrected TESScut time at nearest spacecraft time
    qf=glob.glob(R+f'qlp_71316387_s{sec:04d}.fits')
    if qf:
        q=fits.open(qf[0])[1].data; tq=np.array(q['TIME']); tq=tq[np.isfinite(tq)]
        # QLP time minus FFI time: find TESScut cadence nearest in FFI time
        tt=t[ok]; i=np.clip(np.searchsorted(tt,tq),1,len(tt)-1); j=np.where(np.abs(tt[i]-tq)<np.abs(tt[i-1]-tq),i,i-1)
        dq=(tq-tt[j])*86400
        d['qlp_minus_tesscut_s']=float(np.median(dq))
        d['qlp_minus_true_s']=float(np.median(dq+corr[ok][j]-mine[j]))
    res[sec]=d; print(d, flush=True)
json.dump(res,open(f'{OUT}/bary_check.json','w'),indent=1)
