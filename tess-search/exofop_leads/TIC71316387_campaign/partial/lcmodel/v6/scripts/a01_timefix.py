"""Target-specific barycentric time correction for TESScut (FFI) data.
FFI/TESScut TIME = t_sc(TDB at spacecraft) + TIMECORR_ffi, where TIMECORR_ffi is the Roemer delay for the FFI
reference direction, not the target.  We compute the target Roemer delay from the JPL Horizons barycentric
state of TESS (id -95, ICRF equatorial frame) and validate the method against the target-specific TIMECORR
column of the TESS-SPOC light curves (S1, S28, S68).
Output: ext/horizons_s{sec}.npz (jd_tdb, xyz AU), results/timefix.json (validation + per-sector offsets)."""
import numpy as np, json, glob, os, signal
from astropy.io import fits
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
LOC = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut'
RD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/ctoi_vet/recent/data'
RA, DEC = np.deg2rad(330.887948), np.deg2rad(-42.835511)
N = np.array([np.cos(DEC) * np.cos(RA), np.cos(DEC) * np.sin(RA), np.sin(DEC)])
AU_S = 149597870.700e3 / 299792458.0   # light seconds per AU

def horizons(sec, t0, t1):
    fn = f'{V}/ext/horizons_s{sec}.npz'
    if os.path.exists(fn):
        z = np.load(fn); return z['jd'], z['xyz']
    from astroquery.jplhorizons import Horizons
    signal.alarm(120)
    o = Horizons(id='-95', location='@0', epochs={'start': f'JD{t0 + 2457000 - 0.5:.4f}', 'stop': f'JD{t1 + 2457000 + 0.5:.4f}', 'step': '30m'})
    v = o.vectors(refplane='earth')
    signal.alarm(0)
    jd = np.array(v['datetime_jd'], float); xyz = np.vstack([np.asarray(v[c], float) for c in ('x', 'y', 'z')]).T
    np.savez(fn, jd=jd, xyz=xyz)
    return jd, xyz

def roemer(sec, tsc):
    """Roemer delay (s) towards the target at spacecraft TDB times tsc (BTJD)."""
    jd, xyz = horizons(sec, np.nanmin(tsc), np.nanmax(tsc))
    x = jd - 2457000
    r = np.vstack([np.interp(tsc, x, xyz[:, i]) for i in range(3)]).T
    return r @ N * AU_S

if __name__ == '__main__':
    out = {}
    for sec in (1, 28, 68):
        f = glob.glob(f'{RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
        d = fits.open(f)[1].data
        t = np.array(d['TIME'], float); tc = np.array(d['TIMECORR'], float)
        m = np.isfinite(t) & np.isfinite(tc)
        tsc = t[m] - tc[m]
        mine = roemer(sec, tsc)
        diff = (tc[m] * 86400 - mine)
        out[f'spoc_validation_s{sec}'] = dict(mean_ms=float(diff.mean() * 1e3), rms_ms=float(diff.std() * 1e3), maxabs_ms=float(np.abs(diff).max() * 1e3))
        print(f'S{sec} SPOC TIMECORR - mine: mean {diff.mean()*1e3:.2f} ms, rms {diff.std()*1e3:.2f} ms, max {np.abs(diff).max()*1e3:.2f} ms')
    for fn in sorted(glob.glob(f'{LOC}/s*/*.fits')):
        h = fits.open(fn, memmap=True); sec = int(h[0].header['SECTOR'])
        d = h[1].data
        t = np.array(d['TIME'], float); tc = np.array(d['TIMECORR'], float)
        m = np.isfinite(t) & np.isfinite(tc) & (t > 0)
        tsc = t[m] - tc[m]
        mine = roemer(sec, tsc)
        corr = mine - tc[m] * 86400  # seconds to ADD to TESScut TIME
        np.save(f'{V}/data/timecorr_tc_s{sec}.npy', np.vstack([t[m], corr]))
        out[f'tesscut_s{sec}'] = dict(add_s_mean=float(corr.mean()), add_s_min=float(corr.min()), add_s_max=float(corr.max()))
        print(f'S{sec} TESScut: target-minus-FFI Roemer delay {corr.mean():+.2f} s (range {corr.min():+.2f} .. {corr.max():+.2f})')
    json.dump(out, open(f'{V}/results/timefix.json', 'w'), indent=1)
