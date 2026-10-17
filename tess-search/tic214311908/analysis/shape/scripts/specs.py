"""specs: synthetic photometry for composite-star checks.
Surface-flux spectra: BT-Settl CIFIST (logg 5.0, 2600-4100 K; copies in ../star/data/models made by the star agent)
below 4000 K, Castelli & Kurucz 2004 ATLAS9 ([M/H] = 0, logg 4.5 for 3500-50000 K, logg 5.0 for the 'compact'
option; STScI CDBS ck04models, data/ck04/) above; each spectrum renormalised to sigma T^4. Bands (SVO):
DECam g,r,i,z (AB), SkyMapper u,v,g (AB), TESS (photon-weighted relative flux). Interpolation in log Teff.
S_band(T) = band-averaged surface flux density (erg/s/cm2/Hz for AB bands); absolute AB magnitude of a sphere of
radius R (Rsun) at 10 pc: M = -2.5 log10(S (R/10pc)^2) - 48.6."""
import numpy as np, glob, os
from astropy.io import fits
SD = 'tic214311908/analysis/shape'
BT = 'tic214311908/analysis/star/data/models'
SIG = 5.6704e-5; C_A = 2.998e18; RSUN = 6.957e10; PC = 3.0857e18
BANDS = {'g': 'DECam.g', 'r': 'DECam.r', 'i': 'DECam.i', 'z': 'DECam.z', 'smu': 'SkyMapper.u', 'smv': 'SkyMapper.v', 'smg': 'SkyMapper.g', 'T': 'TESS.Red'}
FILT = {b: np.loadtxt(f'{SD}/data/filt_{f}.dat') for b, f in BANDS.items()}
EXT = {'g': 3.237, 'r': 2.176, 'i': 1.595, 'z': 1.217, 'smu': 4.294, 'smv': 4.026, 'smg': 2.986, 'T': 1.45}  # A_band/E(B-V) (Schlafly & Finkbeiner 2011, R_V = 3.1; TESS ~ between i and z)
def _norm(lam, fl, T):
    m = np.isfinite(fl) & (lam > 0)
    lam, fl = lam[m], fl[m]
    tot = np.trapezoid(fl, lam)
    # Rayleigh-Jeans tail beyond the last point
    tot += fl[-1]*lam[-1]/3
    return lam, fl*SIG*T**4/tot
def spec(T, compact=False):
    """surface-flux spectrum (lam [A], F_lambda) at Teff T"""
    if T < 4000 and not compact:
        grid = sorted(int(os.path.basename(f).split('_')[1]) for f in glob.glob(f'{BT}/btsettl_*_5.0.npy') if 'fbol' not in f)
        out = []
        lo = max([g for g in grid if g <= T] or [grid[0]]); hi = min([g for g in grid if g >= T] or [grid[-1]])
        specs_ = []
        for g in (lo, hi):
            a = np.load(f'{BT}/btsettl_{g}_5.0.npy').astype(float)
            specs_.append(_norm(a[0], a[1], g))
        lam = specs_[0][0]
        if lo == hi: return lam, specs_[0][1]*(T/lo)**4
        f2 = np.interp(lam, specs_[1][0], specs_[1][1])
        wgt = (np.log(T) - np.log(lo))/(np.log(hi) - np.log(lo))
        fl = np.exp((1-wgt)*np.log(np.clip(specs_[0][1], 1e-30, None)) + wgt*np.log(np.clip(f2, 1e-30, None)))
        return _norm(lam, fl, T)
    col = 'g50' if compact else 'g45'
    grid = sorted(int(os.path.basename(f).split('_')[1].split('.')[0]) for f in glob.glob(f'{SD}/data/ck04/ckp00_*.fits'))
    lo = max([g for g in grid if g <= T] or [grid[0]]); hi = min([g for g in grid if g >= T] or [grid[-1]])
    def rd(g):
        d = fits.open(f'{SD}/data/ck04/ckp00_{g}.fits')[1].data
        for c in [col, 'g50', 'g45', 'g40', 'g35']:
            fl = d[c].astype(float)
            if np.all(fl > 0) or np.sum(fl > 0) > 0.9*len(fl): break
        return _norm(d['WAVELENGTH'].astype(float), fl, g)
    l1, f1 = rd(lo)
    if lo == hi: return l1, f1*(T/lo)**4
    l2, f2 = rd(hi); f2 = np.interp(l1, l2, f2)
    wgt = (np.log(T) - np.log(lo))/(np.log(hi) - np.log(lo))
    fl = np.exp((1-wgt)*np.log(np.clip(f1, 1e-30, None)) + wgt*np.log(np.clip(f2, 1e-30, None)))
    return _norm(l1, fl, T)
def band_S(lam, fl, b):
    w, tr = FILT[b][:, 0], FILT[b][:, 1]
    Ti = np.interp(lam, w, tr, left=0, right=0); m = Ti > 0
    if b == 'T':  # photon-weighted mean F_lambda (relative use only)
        return np.trapezoid(fl[m]*Ti[m]*lam[m], lam[m])/np.trapezoid(Ti[m]*lam[m], lam[m])
    # AB: <f_nu> = int f_lam T lam dlam / int T c/lam dlam
    return np.trapezoid(fl[m]*Ti[m]*lam[m], lam[m])/np.trapezoid(Ti[m]*C_A/lam[m], lam[m])
_cache = {}
def S_table(compact=False):
    key = compact
    if key in _cache: return _cache[key]
    Ts = np.unique(np.concatenate([np.arange(2600, 4000, 100), np.array([4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6500, 7000, 7500, 8000, 9000, 10000, 11000, 12000, 14000, 16000, 18000, 20000, 25000, 30000, 35000, 40000, 45000, 50000])]))
    if compact: Ts = Ts[Ts >= 4000]
    tab = {b: np.array([band_S(*spec(T, compact), b) for T in Ts]) for b in BANDS}
    _cache[key] = (Ts, tab); return Ts, tab
def S(T, b, compact=False):
    Ts, tab = S_table(compact)
    return np.exp(np.interp(np.log(T), np.log(Ts), np.log(tab[b])))
def absmag(T, R, b, compact=False):
    return -2.5*np.log10(S(T, b, compact)*(R*RSUN/(10*PC))**2) - 48.6
