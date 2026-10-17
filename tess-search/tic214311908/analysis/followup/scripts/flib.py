"""flib: shared constants and helpers for the follow-up design (folder followup/).
Ephemeris loader (re-runnable with any ephemeris file), target astrometry, observing sites, Moffat PSF helpers."""
import os
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, numpy as np

B = 'tic214311908/analysis'
F = f'{B}/followup'
OUT = f'{F}/out'
DEFAULT_EPH = f'{B}/clean/out/ephemeris_final.json'
FALLBACK_EPH = f'{B}/phot/out/timing.json'

# ---------------- ephemeris ----------------
def load_ephemeris(path=None):
    """Return dict(T0 [BJD_TDB at reference epoch], Eref, P, cov 2x2 [(T0,P)], source, label).
    Accepts clean/out/ephemeris_final.json (uses the 'adopted' block), any json with the same block layout,
    or phase-1 phot/out/timing.json (uses 'all_err_cal', zero T0-P covariance, epoch 0)."""
    if path is None:
        path = DEFAULT_EPH if os.path.exists(DEFAULT_EPH) else FALLBACK_EPH
    d = json.load(open(path))
    if 'adopted' in d:
        blk = d[d['adopted']]
        T0 = blk['T0_ref_BJD_TDB']; Eref = blk['reference_epoch']; P = blk['P_d']
        cov = np.array(blk['cov_T0ref_P'])
        return dict(T0=T0, Eref=Eref, P=P, cov=cov, source=path, label=f"{os.path.basename(path)}:{d['adopted']}")
    if 'all_err_cal' in d:
        blk = d['all_err_cal']
        T0 = blk['T0'] + 2457000.0
        cov = np.diag([blk['T0_err']**2, blk['P_err']**2])
        return dict(T0=T0, Eref=0, P=blk['P'], cov=cov, source=path, label=f"{os.path.basename(path)}:all_err_cal")
    raise ValueError(f'unknown ephemeris layout in {path}')

def predict(eph, E):
    E = np.asarray(E, float)
    dE = E - eph['Eref']
    t = eph['T0'] + eph['P']*dE
    var = eph['cov'][0, 0] + dE**2*eph['cov'][1, 1] + 2*dE*eph['cov'][0, 1]
    return t, np.sqrt(var)

# ---------------- target (Gaia DR3 5925373272850784768) ----------------
TARGET = dict(ra=260.969674, dec=-51.323821, epoch=2016.0, plx=15.497306, pmra=-113.671152, pmdec=85.54355, rv=16.28)
T14_H = 1.70          # adopted planning duration (clean: 1.69-1.73 h)
DEPTH = 0.0152        # adopted TESS-band depth (clean/out/final_depth.json)
RSTAR, MSTAR, TEFF = 0.383, 0.370, 3350.
DIST_PC = 64.31

def target_coord(obstime_jyear=2027.5):
    from astropy.coordinates import SkyCoord, Distance
    from astropy.time import Time
    import astropy.units as u
    c = SkyCoord(ra=TARGET['ra']*u.deg, dec=TARGET['dec']*u.deg, distance=Distance(parallax=TARGET['plx']*u.mas),
                 pm_ra_cosdec=TARGET['pmra']*u.mas/u.yr, pm_dec=TARGET['pmdec']*u.mas/u.yr,
                 radial_velocity=TARGET['rv']*u.km/u.s, obstime=Time(TARGET['epoch'], format='jyear', scale='tcb'))
    return c.apply_space_motion(new_obstime=Time(obstime_jyear, format='jyear', scale='tcb'))

# ---------------- sites (geodetic lat, lon [deg E], height [m]) ----------------
SITES = {
    'CTIO (LCOGT lsc, SMARTS)': (-30.1690, -70.8063, 2207.),
    'Cerro Pachon (SOAR, Gemini-S)': (-30.2408, -70.7367, 2722.),
    'Las Campanas (LCO)': (-29.0146, -70.6926, 2380.),
    'La Silla (TRAPPIST-S, ExTrA, 3.6m)': (-29.2563, -70.7380, 2400.),
    'Paranal (VLT, SPECULOOS-S)': (-24.6272, -70.4042, 2635.),
    'SAAO Sutherland (LCOGT cpt, SALT)': (-32.3794, 20.8107, 1798.),
    'Siding Spring (LCOGT coj, FTS/MuSCAT4, AAT)': (-31.2733, 149.0614, 1165.),
    'Dome C (ASTEP+)': (-75.0997, 123.3325, 3233.),
}
SITE_SHORT = {k: k.split(' (')[0].replace(' ', '_').replace('/', '') for k in SITES}

def jdump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.bool_,)): return bool(o)
        return str(o)
    json.dump(obj, open(path, 'w'), indent=1, default=conv)

# ---------------- Moffat PSF ----------------
def moffat_alpha(fwhm, beta):
    return fwhm/(2*np.sqrt(2**(1/beta) - 1))

def moffat_frac_in_circle(d, r_ap, fwhm, beta, n=400):
    """Fraction of a unit-flux circular Moffat PSF centred at distance d (arcsec) that falls in a circular aperture of
    radius r_ap centred at the origin. Numerical integration in polar coordinates about the aperture centre."""
    a = moffat_alpha(fwhm, beta)
    r = (np.arange(n) + 0.5)/n*r_ap
    th = (np.arange(720) + 0.5)/720*2*np.pi
    R, TH = np.meshgrid(r, th, indexing='ij')
    x = R*np.cos(TH) - d
    y = R*np.sin(TH)
    I = (beta - 1)/(np.pi*a*a)*(1 + (x*x + y*y)/(a*a))**(-beta)
    return float(np.sum(I*R)*(r_ap/n)*(2*np.pi/720))

def moffat_enclosed(r_ap, fwhm, beta):
    a = moffat_alpha(fwhm, beta)
    return 1 - (1 + (r_ap/a)**2)**(1 - beta)

def moffat_peak_frac(fwhm, beta, pix):
    """fraction of flux in the central pixel (pix arcsec) for a PSF centred on the pixel (approximated by the peak
    surface brightness times the pixel area, integrated exactly over a square with 21x21 subsamples)"""
    a = moffat_alpha(fwhm, beta)
    s = (np.arange(21) + 0.5)/21*pix - pix/2
    X, Y = np.meshgrid(s, s)
    I = (beta - 1)/(np.pi*a*a)*(1 + (X*X + Y*Y)/(a*a))**(-beta)
    return float(I.mean()*pix*pix)
