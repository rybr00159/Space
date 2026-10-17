"""Referee-2 helpers for the 'search' audit of TIC 214311908 (independent code; only data products are read).

Conventions: BTJD_TDB; fractional fluxes rescaled to the PRF-unit M3V flux 979 e-/s:
  f = 1 + (935.97/979)(flux_norm - 1).
"""
import os, json, warnings
for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(k, "2")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

BASE = "tic214311908/analysis"
ME = f"{BASE}/referee2_search"
SCALE = 935.97 / 979.0
FREF = 979.0
T0, P = 3840.914139360305, 5.691508695210822      # phot global_B (adopted in phase 1)
SECTORS = (12, 39, 66, 93, 103, 104)
P_RR_GAIA = 0.6485275345290323                      # Gaia DR3 vari_rrlyrae pf of 5925372615693467904
RSTAR, MSTAR, RHO = 0.381, 0.370, 9.4
RE_RSUN = 0.0091577


def jdump(obj, fn):
    def conv(o):
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        raise TypeError(type(o))
    with open(fn, "w") as f:
        json.dump(obj, f, indent=1, default=conv)


def rstd(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def phase_dist(t, P_=P, T0_=T0):
    """signed time from nearest mid-transit (d)."""
    return (((t - T0_) / P_ + 0.5) % 1 - 0.5) * P_


def load_m3v(model="B"):
    d = pd.read_csv(f"{BASE}/phot/lc_M3V_all.csv")
    d = d[d.quality_ok == 1].copy()
    if model == "B":
        det, raw = d.flux_norm.values, d.flux_raw_norm.values
    else:
        det, raw = d.flux_norm_modelA.values, d.flux_raw_norm_modelA.values
    out = pd.DataFrame(dict(t=d.time_btjd_tdb.values, s=d.sector.values.astype(int), cad=d.cadence_s.values.astype(float),
                            det=1 + SCALE * (det - 1), raw=1 + SCALE * (raw - 1), dx=d.jitter_dx_px.values,
                            dy=d.jitter_dy_px.values, bkg=d.bkg_e_s_px.values))
    out = out[np.isfinite(out.det)].reset_index(drop=True)
    out["intr"] = np.abs(phase_dist(out.t.values)) < 0.1
    return out


def binlc(t, y, s, cad, wmin=10.0):
    """bin within sector to max(wmin, cadence); returns t, y, n, s."""
    T, Y, N, S = [], [], [], []
    for sec in np.unique(s):
        m = s == sec
        w = max(wmin * 60, np.median(cad[m])) / 86400
        k = np.floor((t[m] - t[m].min()) / w).astype(int)
        u, inv = np.unique(k, return_inverse=True)
        n = np.bincount(inv)
        T.append(np.bincount(inv, t[m]) / n); Y.append(np.bincount(inv, y[m]) / n); N.append(n); S.append(np.full(len(u), sec))
    return np.concatenate(T), np.concatenate(Y), np.concatenate(N), np.concatenate(S)


def harmX(t, Pp, nh=3):
    c = [np.ones(len(t))]
    for h in range(1, nh + 1):
        c += [np.sin(2 * np.pi * h * t / Pp), np.cos(2 * np.pi * h * t / Pp)]
    return np.column_stack(c)


def trapz_model(t, tc, T14=1.6 / 24, flat=0.7, cad_s=None):
    """unit-depth trapezoid (1 in the flat part), exposure-averaged with 9 sub-samples."""
    if cad_s is None: cad_s = np.full(len(t), 120.0)
    out = np.zeros(len(t))
    for c in np.unique(cad_s):
        m = cad_s == c
        acc = np.zeros(m.sum())
        for o in np.linspace(-0.5, 0.5, 9) * c / 86400:
            x = np.abs(t[m] + o - tc)
            ing = 0.5 * T14 * (1 - flat)
            acc += np.clip((0.5 * T14 - x) / ing, 0, 1)
        out[m] = acc / 9
    return out


def trapz_periodic(t, Pp, T0_, cad_s, T14=1.6 / 24, flat=0.7):
    x = phase_dist(t, Pp, T0_)
    return trapz_model(x, 0.0, T14, flat, cad_s)


def a_R(Pd, rho=RHO):
    return (6.674e-8 * rho * (Pd * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def t14c(Pd, rho=RHO):
    return Pd / np.pi * np.arcsin(np.clip(1 / a_R(Pd, rho), 0, 1))
