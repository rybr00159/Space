"""Synthetic photometry helpers: BT-Settl (CIFIST) grid, passbands, Mamajek dwarf table, LD.

Passbands (photon-counting response S(lambda); band-averaged F_lambda = int F S lam dlam / int S lam dlam):
  TESS, Johnson B/V (Bessell), Cousins R/I, SDSS g, Gaia G/BP/RP, 2MASS Ks, WISE W1/W2 from SVO;
  WASP = 400-700 nm filter (Pollacco et al. 2006) x e2v CCD42-40 back-illuminated QE (approx.), with
  alternatives for systematics.
"""
import re
import numpy as np
import lib3 as L

M = L.CH + '/data/models/'
_grid = None


def grid():
    global _grid
    if _grid is None:
        z = np.load(M + 'btsettl_grid.npz')
        _grid = dict(wave=z['wave'], edges=z['edges'], teff=z['teff'], logg=z['logg'], flux=z['flux'])
    return _grid


def spectrum(teff):
    """surface F_lambda (erg/s/cm2/A) at teff, log-linear interpolation in the grid (grid logg 5.0 for
    T<=4500, 4.5 above, i.e. dwarf gravities)."""
    g = grid()
    T = g['teff']
    teff = float(np.clip(teff, T.min(), T.max()))
    j = np.searchsorted(T, teff)
    if j == 0:
        return g['flux'][0]
    if T[j - 1] == teff:
        return g['flux'][j - 1]
    if j >= len(T):
        return g['flux'][-1]
    x = (teff - T[j - 1]) / (T[j] - T[j - 1])
    return np.exp((1 - x) * np.log(np.maximum(g['flux'][j - 1], 1e-30)) + x * np.log(np.maximum(g['flux'][j], 1e-30)))


def _svo(fn):
    s = open(M + fn).read()
    rows = re.findall(r'<TR>\s*<TD>([^<]+)</TD>\s*<TD>([^<]+)</TD>\s*</TR>', s)
    a = np.array(rows, float)
    o = np.argsort(a[:, 0])
    return a[o, 0], np.clip(a[o, 1], 0, None)


def wasp_band(kind='nominal'):
    lam = np.arange(3500, 7600, 10.0)
    qe_l = np.array([3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500])
    qe_v = np.array([0.30, 0.60, 0.75, 0.85, 0.90, 0.90, 0.88, 0.85, 0.78])
    qe = np.interp(lam, qe_l, qe_v)
    if kind == 'nominal':        # 400-700 nm, 10-nm linear edges
        filt = np.clip(np.minimum((lam - 3950) / 100, (7050 - lam) / 100), 0, 1)
        return lam, filt * qe
    if kind == 'flat':           # pure top-hat 400-700, no QE
        return lam, ((lam >= 4000) & (lam <= 7000)).astype(float)
    if kind == 'red':            # 420-700 nm (redder blue edge)
        filt = np.clip(np.minimum((lam - 4150) / 100, (7050 - lam) / 100), 0, 1)
        return lam, filt * qe
    if kind == 'blue':           # 400-680 nm
        filt = np.clip(np.minimum((lam - 3950) / 100, (6850 - lam) / 100), 0, 1)
        return lam, filt * qe
    raise ValueError(kind)


BANDS = {'TESS': 'filt_TESS_TESS.Red.xml', 'B': 'filt_Generic_Johnson.B.xml', 'V': 'filt_Generic_Johnson.V.xml',
         'R': 'filt_Generic_Cousins.R.xml', 'I': 'filt_Generic_Cousins.I.xml', 'g': 'filt_SLOAN_SDSS.g.xml',
         'G': 'filt_GAIA_GAIA3.G.xml', 'BP': 'filt_GAIA_GAIA3.Gbp.xml', 'RP': 'filt_GAIA_GAIA3.Grp.xml',
         'Ks': 'filt_2MASS_2MASS.Ks.xml', 'W1': 'filt_WISE_WISE.W1.xml', 'W2': 'filt_WISE_WISE.W2.xml'}
_bands = {}


def band(name):
    if name not in _bands:
        if name.startswith('WASP'):
            kind = name.split('_')[1] if '_' in name else 'nominal'
            _bands[name] = wasp_band(kind)
        else:
            _bands[name] = _svo(BANDS[name])
    return _bands[name]


def band_flux(spec, name):
    """photon-weighted mean F_lambda of a spectrum on the grid in band `name`"""
    g = grid()
    lam, S = band(name)
    w = g['wave']
    Si = np.interp(w, lam, S, left=0, right=0)
    dl = np.diff(g['edges'])
    return np.sum(spec * Si * w * dl) / np.sum(Si * w * dl)


def lam_eff(spec, name):
    g = grid(); lam, S = band(name); w = g['wave']; Si = np.interp(w, lam, S, left=0, right=0); dl = np.diff(g['edges'])
    return np.sum(spec * Si * w * w * dl) / np.sum(spec * Si * w * dl)


# Vega-system zero points (F_lambda, erg/s/cm2/A) for colour checks only (Bessell et al. 1998; Cohen 2003)
ZP = {'B': 6.32e-9, 'V': 3.631e-9, 'R': 2.177e-9, 'I': 1.126e-9, 'Ks': 4.283e-11}


def mamajek():
    rows = {}
    hdr = None
    for line in open(M + 'EEM_dwarf.txt'):
        if line.startswith('#SpT'):
            if hdr is None:
                hdr = line[1:].split()
            continue
        p = line.split()
        if hdr is None or len(p) < 10 or not re.match(r'^[OBAFGKMLTY]\d(\.\d)?V$', p[0]):
            continue
        d = {}
        for k, v in zip(hdr, p):
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = np.nan
        rows[p[0]] = d
    return rows


# --------------------------------------------------------------------------- limb darkening
_ld = None


def ld_quadratic(teff, logg, bandname):
    """quadratic LD (a, b): Claret (2011) ATLAS L-method Z=0 for B,V,R,I,K,g*; Claret (2017) TESS ATLAS.
    WASP: interpolate in effective wavelength between B (440 nm), V (550 nm) and R (640 nm) -> uses V and B/R mix;
    W1/W2: Claret K values (LD is weak in the IR)."""
    global _ld
    import pandas as pd
    if _ld is None:
        c11 = pd.read_csv(M + 'claret2011_quadratic.csv')
        c11 = c11[(c11.Mod == 'A') & (c11.Met == 'L') & (c11.Z == 0.0)]
        c17 = pd.read_csv(L.BASE + '/tic713/lcmodel/data/claret2017_tess_quadratic_atlas.csv')
        c17 = c17[(c17.Z == 0.0)]
        _ld = (c11, c17)
    c11, c17 = _ld

    def interp(tab, acol, bcol):
        tab = tab[np.isclose(tab.logg, min(max(logg, 4.5), 5.0)) | np.isclose(tab.logg, 4.5 if logg < 4.75 else 5.0)]
        tab = tab[np.isclose(tab.logg, (4.5 if logg < 4.75 else 5.0))]
        tab = tab.sort_values('Teff')
        return float(np.interp(teff, tab.Teff, tab[acol])), float(np.interp(teff, tab.Teff, tab[bcol]))
    if bandname == 'TESS':
        return interp(c17, 'aLSM', 'bLSM')
    m = {'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'Ks': 'K', 'W1': 'K', 'W2': 'K', 'g': 'g*'}
    if bandname in m:
        return interp(c11[c11.Filt == m[bandname]], 'a', 'b')
    if bandname.startswith('WASP'):
        aB, bB = interp(c11[c11.Filt == 'B'], 'a', 'b')
        aV, bV = interp(c11[c11.Filt == 'V'], 'a', 'b')
        aR, bR = interp(c11[c11.Filt == 'R'], 'a', 'b')
        # effective wavelength of WASP band for this Teff
        le = lam_eff(spectrum(teff), bandname)
        lB, lV, lR = 4400., 5500., 6400.
        if le < lV:
            x = (le - lB) / (lV - lB); return (1 - x) * aB + x * aV, (1 - x) * bB + x * bV
        x = (le - lV) / (lR - lV); return (1 - x) * aV + x * aR, (1 - x) * bV + x * bR
    raise ValueError(bandname)


# --------------------------------------------------------------------------- eclipse of a limb-darkened disk
def occulted_fraction(d, R1, R2, ld, nr=4000):
    """fraction of the flux of star 1 (radius R1, quadratic LD (a,b)) hidden by an opaque disk of radius R2
    at projected separation d (same units)."""
    a, b = ld
    r = (np.arange(nr) + 0.5) / nr * R1
    dr = R1 / nr
    mu = np.sqrt(1 - (r / R1) ** 2)
    I = 1 - a * (1 - mu) - b * (1 - mu) ** 2
    tot = np.sum(I * 2 * np.pi * r * dr)
    if d >= R1 + R2:
        return 0.0
    # angular fraction of the annulus at radius r inside the occultor
    with np.errstate(invalid='ignore', divide='ignore'):
        c = (r ** 2 + d ** 2 - R2 ** 2) / (2 * r * max(d, 1e-12))
    frac = np.where(c <= -1, 1.0, np.where(c >= 1, 0.0, np.arccos(np.clip(c, -1, 1)) / np.pi))
    if d < 1e-9:
        frac = (r < R2).astype(float)
    return float(np.sum(I * 2 * np.pi * r * dr * frac) / tot)


# --------------------------------------------------------------------------- empirically warped spectra
ANCHORS = ('B', 'V', 'R', 'I', 'Ks')
COLKEY = {'B': ('B-V', -1), 'R': ('V-Rc', +1), 'I': ('V-Ic', +1), 'Ks': ('V-Ks', +1)}


def synmag(spec, b):
    return -2.5 * np.log10(band_flux(spec, b) / ZP[b])


def warp(spec, target, niter=8):
    """multiply the spectrum by a smooth (piecewise-linear in lambda, in magnitudes) correction so that
    its synthetic B-V, V-Rc, V-Ic, V-Ks equal `target` (dict with those keys; missing keys -> no anchor).
    Returns warped spectrum and final residuals."""
    g = grid(); w = g['wave']
    sp = spec.copy()
    anchors = [b for b in ANCHORS if b == 'V' or (COLKEY[b][0] in target and np.isfinite(target[COLKEY[b][0]]))]
    for _ in range(niter):
        V = synmag(sp, 'V')
        lam = np.array([lam_eff(sp, b) for b in anchors])
        dm = []
        for b in anchors:
            if b == 'V':
                dm.append(0.0); continue
            key, sgn = COLKEY[b]
            want = V - sgn * target[key]          # B = V + (B-V); R = V - (V-R) ...
            dm.append(want - synmag(sp, b))
        dm = np.array(dm)
        o = np.argsort(lam)
        corr = np.interp(w, lam[o], dm[o])        # constant beyond the end anchors
        sp = sp * 10 ** (-0.4 * corr)
    V = synmag(sp, 'V')
    res = {COLKEY[b][0]: float((V - synmag(sp, b)) * COLKEY[b][1] - target[COLKEY[b][0]]) for b in anchors if b != 'V'}
    return sp, res


def stassun_T(G, c):
    """TIC-8 (Stassun et al. 2019) relation T = G - f(Bp-Rp) for dwarfs"""
    return G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473


def star_model(spt=None, teff=None, tab=None):
    """warped spectrum + absolute V for a Mamajek spectral type, or interpolated in Teff"""
    if tab is None:
        tab = mamajek()
    if spt is not None:
        r = tab[spt]
    else:
        keys = [k for k in tab if np.isfinite(tab[k]['Teff']) and 2900 < tab[k]['Teff'] < 7200]
        Ts = np.array([tab[k]['Teff'] for k in keys]); o = np.argsort(Ts); keys = [keys[i] for i in o]; Ts = Ts[o]
        r = {}
        for c in ('Teff', 'Mv', 'B-V', 'V-Rc', 'V-Ic', 'V-Ks', 'R_Rsun', 'Msun', 'Bp-Rp', 'G-V', 'M_G'):
            v = np.array([tab[k].get(c, np.nan) for k in keys])
            ok = np.isfinite(v)
            r[c] = float(np.interp(teff, Ts[ok], v[ok]))
        r['Teff'] = teff
    sp0 = spectrum(r['Teff'])
    sp, res = warp(sp0, r)
    return dict(row=r, spec=sp, spec_raw=sp0, warp_res=res)
