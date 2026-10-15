"""Light-curve loading, outlier clipping, eclipse-window extraction and noise estimation (v5).
Products:
  SPOC : TESS-SPOC FFI light curves (PDCSAP, crowding-corrected by SPOC), S1 (1800 s), S28 (600 s), S68 (200 s)
  TC   : own TESScut aperture photometry (a02_phot.py; target barycentric times; PRF contamination removed)
  QLP  : QLP SAP flux (S1/S28/S68/S95/S105; times shifted by the measured QLP - target-BJD offset; S105 QLP times
         drift by ~-5 s/day relative to SPOC FFI times and are not used for timing)"""
import numpy as np, glob, json
from astropy.io import fits
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
RD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/ctoi_vet/recent/data'
P0 = 1.3036712; T00 = 3165.6697; T14 = 1.55 / 24
TEXP = {1: 1800 / 86400, 28: 600 / 86400, 68: 200 / 86400, 95: 200 / 86400, 105: 200 / 86400, 106: 200 / 86400}
QLP_OFF_S = {1: 1.85, 28: 1.98, 68: 1.98, 95: 1.96}   # QLP TIME - target BJD_TDB (a03 check)

def load_spoc(sec, col='PDCSAP_FLUX'):
    f = glob.glob(f'{RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
    h = fits.open(f); assert int(h[0].header['TICID']) == 71316387 and int(h[0].header['SECTOR']) == sec
    d = h[1].data
    t = np.array(d['TIME'], float); fl = np.array(d[col], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    t, fl = t[m], fl[m]
    return dict(t=t, f=fl / np.median(fl), sector=sec, prod='SPOC', texp=TEXP[sec], crowd=float(h[1].header['CROWDSAP']),
                flfr=float(h[1].header['FLFRCSAP']))

def load_tc(sec, raw=False):
    z = np.load(f'{V}/data/tc_s{sec}.npz')
    return dict(t=z['t'], f=(z['f_raw'] if raw else z['f']).copy(), sector=sec, prod='TC', texp=TEXP[sec], contam=float(z['contam']),
                sky=z['skyap'] / z['flux'], cx=z['cx'], cy=z['cy'])

def load_qlp(sec):
    f = f'{RD}/qlp_71316387_s{sec:04d}.fits'
    h = fits.open(f); assert int(h[0].header['TICID']) == 71316387 and int(h[0].header['SECTOR']) == sec
    d = h[1].data
    t = np.array(d['TIME'], float); fl = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    off = QLP_OFF_S.get(sec, np.nan)
    return dict(t=t[m] - (off if np.isfinite(off) else 0.0) / 86400, f=fl[m] / np.median(fl[m]), sector=sec, prod='QLP',
                texp=TEXP[sec], time_ok=bool(np.isfinite(off)))

def phase_dist(t, P=P0, T0=T00):
    ph = ((t - T0) / P) % 1.0
    return np.minimum(ph, 1 - ph) * P, np.abs(ph - 0.5) * P

def ecl_mask(t, w=0.8, P=P0, T0=T00):
    d1, d2 = phase_dist(t, P, T0)
    return (d1 < w * T14) | (d2 < w * T14)

def running_median(t, f, w):
    out = np.empty_like(f)
    j0 = 0; j1 = 0; n = t.size
    for i in range(n):
        while t[j0] < t[i] - w / 2: j0 += 1
        while j1 < n and t[j1] <= t[i] + w / 2: j1 += 1
        out[i] = np.median(f[j0:j1])
    return out

def clip(lc, nsig=5.0, w=0.3, P=P0, T0=T00):
    """5-sigma clip against a 0.3-d running median of the eclipse-masked flux; inside eclipses only upward outliers
    (> 5 sigma + 1%) are removed so that no eclipse point is clipped for being low."""
    t, f = lc['t'], lc['f']
    m = ecl_mask(t, P=P, T0=T0)
    rm = np.interp(t, t[~m], running_median(t[~m], f[~m], w))
    r = f - rm
    s = 1.4826 * np.median(np.abs(r[~m] - np.median(r[~m])))
    bad = (~m & (np.abs(r) > nsig * s)) | (m & (r > nsig * s + 0.01))
    out = {k: (v[~bad] if isinstance(v, np.ndarray) and v.shape == t.shape else v) for k, v in lc.items()}
    out['nclip'] = int(bad.sum()); out['sig_pt'] = float(s)
    return out

def windows(lc, P=P0, T0=T00, hw=0.20, phases=(0.0, 0.5), min_cov=0.7, shift=0.0):
    """Eclipse windows of +-hw d. Keep an event if the central |dt| < T14/2 holds >= min_cov of the expected
    cadences, each out-of-eclipse side (0.8 T14 < |dt| < hw) holds >= 40%, and no gap > 3.5 cadences lies within
    |dt| < T14/2 + 0.02 d."""
    t = lc['t']; cad = lc['texp']
    out = []
    n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
    for n in range(n0, n1 + 1):
        for ph in phases:
            tc = T0 + (n + ph + shift) * P
            idx = np.where(np.abs(t - tc) < hw)[0]
            if idx.size < 10: continue
            dt = t[idx] - tc
            nin = (np.abs(dt) < T14 / 2).sum()
            nl = (dt < -0.8 * T14).sum(); nr = (dt > 0.8 * T14).sum(); side = (hw - 0.8 * T14) / cad
            if nin < min_cov * T14 / cad or nl < 0.4 * side or nr < 0.4 * side: continue
            core = dt[np.abs(dt) < T14 / 2 + 0.02]
            if core.size > 1 and np.max(np.diff(core)) > 3.5 * cad: continue
            out.append(dict(n=n, ph=ph, tc=tc, idx=idx, sector=lc['sector'], prod=lc['prod']))
    return out

def beta_red(r, cad, scales_min=(30, 45, 60, 90)):
    """time-averaging (Pont et al. 2006) beta of a residual series r, median over bin scales"""
    b = []
    s1 = np.std(r)
    for sc in scales_min:
        n = int(round(sc / (cad * 1440)))
        if n < 2: continue
        m = r.size // n
        if m < 10: continue
        rb = r[:m * n].reshape(m, n).mean(1)
        b.append(np.std(rb) / (s1 / np.sqrt(n) * np.sqrt(m / (m - 1))))
    return float(np.median(b)) if b else 1.0

def build(lcs, hw=0.20, P=P0, T0=T00, verbose=True, inflate=True):
    """Concatenate eclipse windows of several light curves into one fitting data set."""
    cols = {k: [] for k in ('t', 'y', 'sig', 'ev', 'dt', 'sector', 'texp', 'ph', 'tc')}
    evtab = []; info = {}; eid = 0
    for lc in lcs:
        lc = clip(lc, P=P, T0=T0)
        ws = windows(lc, P=P, T0=T0, hw=hw)
        res = []
        for w in ws:
            i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.8 * T14
            c = np.polyfit(dt[oot], y[oot], 2); res.append(y[oot] / np.polyval(c, dt[oot]) - 1)
        r = np.concatenate(res)
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        # beta: residuals of each window side treated as short continuous series
        bet = beta_red(np.concatenate([x for x in res if x.size >= 5]), lc['texp']) if inflate else 1.0
        beta = max(1.0, bet)
        key = f"{lc['prod']}{lc['sector']}"
        info[key] = dict(prod=lc['prod'], sector=lc['sector'], sigma=float(sig), beta=float(bet), nev=len(ws), nclip=lc['nclip'],
                         n_prim=sum(1 for w in ws if w['ph'] == 0.0), n_sec=sum(1 for w in ws if w['ph'] == 0.5))
        if verbose:
            print(f"  {key}: {len(ws)} windows ({info[key]['n_prim']} prim, {info[key]['n_sec']} sec), sigma {sig*1e6:.0f} ppm/pt, beta {bet:.2f}, clipped {lc['nclip']}", flush=True)
        for w in ws:
            i = w['idx']
            cols['t'].append(lc['t'][i]); cols['y'].append(lc['f'][i]); cols['sig'].append(np.full(i.size, sig * beta))
            cols['ev'].append(np.full(i.size, eid)); cols['dt'].append(lc['t'][i] - w['tc']); cols['sector'].append(np.full(i.size, lc['sector']))
            cols['texp'].append(np.full(i.size, lc['texp'])); cols['ph'].append(np.full(i.size, w['ph'])); cols['tc'].append(np.full(i.size, w['tc']))
            evtab.append(dict(eid=eid, n=int(w['n']), ph=float(w['ph']), tc=float(w['tc']), sector=int(lc['sector']), prod=lc['prod'], key=key))
            eid += 1
    d = {k: np.concatenate(v) for k, v in cols.items()}
    d['evtab'] = evtab; d['info'] = info
    ev = d['ev']; d['starts'] = np.r_[0, np.where(np.diff(ev) != 0)[0] + 1]; d['ends'] = np.r_[d['starts'][1:], ev.size]
    d['key'] = np.array([evtab[e]['key'] for e in ev])
    return d
