"""Light-curve loading and eclipse-window extraction for TIC 71316387."""
import numpy as np, glob, json
from astropy.io import fits
HERE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel'
RD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/ctoi_vet/recent/data'
P0 = 1.3036712; T00 = 3165.6697; T14 = 1.55 / 24
EXPT = {1: 1800/86400, 28: 600/86400, 68: 200/86400, 95: 200/86400, 105: 200/86400, 106: 200/86400}
NSUB = {1: 15, 28: 7, 68: 3, 95: 3, 105: 3, 106: 3}

def load_spoc(sec):
    f = glob.glob(f'{RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
    h = fits.open(f); assert int(h[0].header['TICID']) == 71316387 and int(h[0].header['SECTOR']) == sec
    d = h[1].data
    t = np.array(d['TIME'], float); fl = np.array(d['PDCSAP_FLUX'], float); fe = np.array(d['PDCSAP_FLUX_ERR'], float)
    q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    t, fl, fe = t[m], fl[m], fe[m]
    med = np.median(fl)
    return dict(t=t, f=fl / med, e=fe / med, sector=sec, prod='SPOC-PDC', texp=EXPT[sec], nsub=NSUB[sec],
                crowd=float(h[1].header['CROWDSAP']))

def load_tc(sec, contam_correct=True):
    z = np.load(f'{HERE}/data/lc_tc_s{sec}.npz')
    t, f = z['t'], z['f']
    c = float(z['contam'])
    if contam_correct:
        f = (f - c) / (1 - c)
    return dict(t=t, f=f, e=np.full(t.size, np.nan), sector=sec, prod='TC-AP', texp=EXPT[sec], nsub=NSUB[sec], contam=c)

def load_qlp(sec):
    f = f'{RD}/qlp_71316387_s{sec:04d}.fits'
    h = fits.open(f); assert int(h[0].header['TICID']) == 71316387
    d = h[1].data
    t = np.array(d['TIME'], float); fl = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    bc = json.load(open(f'{HERE}/data/bary_check.json'))
    off = bc[str(sec)].get('qlp_minus_true_s', 0.0)
    return dict(t=t[m] - off / 86400, f=fl[m] / np.median(fl[m]), e=np.full(m.sum(), np.nan), sector=sec, prod='QLP-SAP',
                texp=EXPT[sec], nsub=NSUB[sec])

def default_set():
    """Primary analysis set: TESS-SPOC PDCSAP for S1/S28/S68, own TESScut aperture photometry for S95/S105/S106."""
    return [load_spoc(1), load_spoc(28), load_spoc(68), load_tc(95), load_tc(105), load_tc(106)]

def windows(lc, P=P0, T0=T00, hw=0.20, phases=(0.0, 0.5), ecl_hw=0.75 * T14, min_cov=0.7, shift=0.0):
    """Cut a +-hw window around each predicted eclipse (phase 0 and 0.5 + shift).
    Keep the event if (a) the in-eclipse part |dt|<T14/2 has >= min_cov of the expected cadences, and
    (b) each side of the out-of-eclipse baseline (ecl_hw<|dt|<hw) has >= 40% of expected cadences, and
    (c) no data gap > 3 cadences inside |dt|<T14/2+0.02 d."""
    t = lc['t']; cad = lc['texp']
    out = []
    n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
    for n in range(n0, n1 + 1):
        for ph in phases:
            tc = T0 + (n + ph + shift) * P
            m = np.abs(t - tc) < hw
            if m.sum() < 10:
                continue
            dt = t[m] - tc
            nin = (np.abs(dt) < T14 / 2).sum(); exp_in = T14 / cad
            nl = ((dt < -ecl_hw)).sum(); nr = ((dt > ecl_hw)).sum(); exp_side = (hw - ecl_hw) / cad
            if nin < min_cov * exp_in or nl < 0.4 * exp_side or nr < 0.4 * exp_side:
                continue
            core = dt[np.abs(dt) < T14 / 2 + 0.02]
            if core.size > 1 and np.max(np.diff(core)) > 3.5 * cad:
                continue
            idx = np.where(m)[0]
            out.append(dict(n=n, ph=ph, tc=tc, idx=idx, sector=lc['sector'], prod=lc['prod']))
    return out

def clip_outliers(lc, P=P0, T0=T00, nsig=5.0, w=0.3):
    """Remove isolated outliers: >nsig robust sigma from a running median (window w d) – applied to all points,
    but eclipse points are compared with the running median of eclipse-masked residual scatter only."""
    t, f = lc['t'], lc['f']
    keep = np.ones(t.size, bool)
    rm = np.empty_like(f)
    for i in range(t.size):
        s = np.abs(t - t[i]) < w / 2
        rm[i] = np.median(f[s])
    r = f - rm
    ph = ((t - T0) / P) % 1
    ecl = (np.minimum(ph, 1 - ph) * P < 0.8 * T14) | (np.abs(ph - 0.5) * P < 0.8 * T14)
    sig = 1.4826 * np.median(np.abs(r[~ecl] - np.median(r[~ecl])))
    # out of eclipse: symmetric clip; in eclipse: only clip positive (upward) outliers beyond 5 sig + depth margin
    bad = (~ecl & (np.abs(r) > nsig * sig)) | (ecl & (r > nsig * sig + 0.01))
    keep &= ~bad
    d = {k: (v[keep] if isinstance(v, np.ndarray) and v.size == t.size else v) for k, v in lc.items()}
    d['sig_pt'] = sig; d['nclip'] = int(bad.sum())
    return d
