"""panel_fpp library: imports referee2_fpp/scripts/r05_fpp.py READ-ONLY (nothing is written to referee2_fpp/) and applies
the photo-eccentric sign fix (rho_c = rho_host * g^+3) in the same way as referee3_bd/s02, with planet and BD e laws.
Own wrapper; used to (i) reproduce the adopted sign-fixed numbers, (ii) re-run the prior-dimension variants and joint
prior draws that were only run before the sign fix."""
import os, sys
for _k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'): os.environ[_k] = '2'
import numpy as np
B = 'tic214311908/analysis'
sys.path.insert(0, f'{B}/referee2_fpp/scripts')
_cwd = os.getcwd(); os.chdir(f'{B}/referee2_fpp')
import r05_fpp as R
os.chdir(_cwd)
N_E = int(os.environ.get('PF_NE', 5000))
def edraw(kind, n, rng):
    if kind == 'kipping': return rng.beta(0.867, 3.03, n)
    if kind == 'vaneylen': return np.abs(rng.normal(0, 0.32, n))
    if kind == 'circ': return np.zeros(n)
    if kind == 'rayleigh19': return rng.rayleigh(0.19, n)
    if kind == 'bdsp':
        a = np.clip(rng.normal(0.426, 0.092, n), 0.05, None); b = np.clip(rng.normal(2.354, 0.73, n), 0.1, None); return rng.beta(a, b)
    if kind == 'bdlp':
        a = np.clip(rng.normal(1.879, 0.92, n), 0.05, None); b = np.clip(rng.normal(2.470, 1.25, n), 0.1, None); return rng.beta(a, b)
    raise ValueError(kind)
def ecc_fixed(kind, n=2500):
    rng = np.random.default_rng(12345)
    e = np.clip(edraw(kind, N_E, rng), 0, 0.9); w = rng.uniform(0, 2*np.pi, N_E)
    g = (1 + e*np.sin(w))/np.sqrt(1 - e**2); geo = (1 + e*np.sin(w))/(1 - e**2)
    return 1/g, geo          # r05 evaluates rho_h * ge^-3 -> rho_h * g^+3 (correct sign)
_orig_build = R.build
def build_fixed(ds, kind, alpha):
    if '+' in kind:
        kp, kb = kind.split('+')
        tp = _orig_build(ds, kp, alpha); tb = _orig_build(ds, kb, alpha)
        out = dict(tp); out['BD'] = tb['BD']; return out
    return _orig_build(ds, kind, alpha)
R.ecc = ecc_fixed; R.build = build_fixed
def run(p_over, rvf=1.0, seed=11):
    R.rng = np.random.default_rng(seed)
    p = dict(R.BASE); p.update(p_over)
    r = R.evaluate(p)
    po = dict(r['post'])
    if rvf != 1.0:
        po['BD_on_M3V_single'] *= rvf; po['BD_on_primary_of_binary'] *= rvf
        s = sum(po.values()); po = {k: v/s for k, v in po.items()}
    planet = po['TP_single'] + po['TP_diluted'] + po['TP_on_twin']
    bd = po['BD_on_M3V_single'] + po['BD_on_primary_of_binary'] + po['BDHEB_on_companion']
    return dict(FPP=1 - planet, FPP_strict=1 - po['TP_single'] - po['TP_diluted'], NFPP=po['NTP_pairA'],
                FPP_system=1 - planet - po['PTP_companion'], PTP=po['PTP_companion'], BD=bd, NTP=po['NTP_pairA'],
                WD=po['EB_WD_on_M3V'], HEB=po['HEB_MSMS'], BG=po['BG_all'], TP_twin=po['TP_on_twin'],
                P_planet_bound_incl_pairA=planet + po['PTP_companion'] + po['NTP_pairA'])
ADOPT = dict(ecc='kipping+bdsp', f_BD=1.3e-4)
