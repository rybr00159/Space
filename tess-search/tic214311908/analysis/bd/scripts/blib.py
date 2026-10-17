"""blib: shared helpers for the BD / companion-hosted channel analysis of TIC 214311908 (folder bd/).
Own code. Reads phase-1/phase-2 data products read-only; writes only into bd/."""
import os
for _k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'): os.environ.setdefault(_k, '1')
import json, numpy as np
B = 'tic214311908/analysis'
BD = f'{B}/bd'
G = 6.674e-8; RSUN = 6.957e10; MSUN = 1.989e33; AU = 1.496e13; DAY = 86400.0; REARTH = 6.371e8
RJUP = 7.1492e9; MJUP = 1.898e30; YR = 3.156e7
T0, P = 3840.9141, 5.6915087          # phase-1 ephemeris (phot/out/timing.json)
FM3V = 979.0                          # M3V flux in PRF units (e-/s), referee_phot
RSTAR, MSTAR, RHO_STAR, RHO_ERR = 0.383, 0.370, 9.3, 0.9   # referee_star (single star; R 0.381-0.385)
DIST = 64.31
U1, U2, SU = 0.125, 0.533, 0.10       # quadratic LD of the M3V in the TESS band (as referee2_shape / referee2_fpp)

def rstd(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826*np.median(np.abs(x - np.median(x)))

def aR_of_rho(rho, per=P):
    return (np.asarray(rho)*G*(per*DAY)**2/(3*np.pi))**(1/3.)

def jdump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        return str(o)
    json.dump(obj, open(path, 'w'), indent=1, default=conv)

def jload(path):
    return json.load(open(path))

def g_ecc(e, w):
    """photo-eccentric factor g = (1 + e sin w)/sqrt(1 - e^2); rho_circ = rho_true * g^3 (Dawson & Johnson 2012);
    transit probability factor (1 + e sin w)/(1 - e^2)."""
    return (1 + e*np.sin(w))/np.sqrt(1 - e**2), (1 + e*np.sin(w))/(1 - e**2)
