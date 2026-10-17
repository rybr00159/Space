"""rlib3: shared constants/helpers for the referee-3 check of the bd/ analysis (own code)."""
import os
for _k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'): os.environ[_k] = '1'
import json, numpy as np
B = 'tic214311908/analysis'
R3 = f'{B}/referee3_bd'
G = 6.6743e-8; RSUN = 6.957e10; MSUN = 1.98892e33; AU = 1.495978707e13; DAY = 86400.0; REARTH = 6.3781e8
RJUP = 7.1492e9; MJUP = 1.89813e30; YR = 3.15576e7
# final cleaned ephemeris (clean/out/ephemeris_final.json, combined): T0 3840.91488 BTJD_TDB, P 5.6915020
T0C, PC = 3840.91488, 5.6915020
MSTAR, RSTAR, RHO_STAR, RHO_ERR = 0.370, 0.383, 9.42, 0.88     # star/ (Mann+15/19), referee_star R 0.381-0.385
U1, U2, SU = 0.125, 0.533, 0.10      # quadratic LD (TESS) used by referee2_shape / referee2_fpp / bd

def jdump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        return str(o)
    json.dump(obj, open(path, 'w'), indent=1, default=conv)

def rstd(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826*np.median(np.abs(x - np.median(x)))

def aR_of_rho(rho, per=PC):
    return (np.asarray(rho)*G*(per*DAY)**2/(3*np.pi))**(1/3.)
