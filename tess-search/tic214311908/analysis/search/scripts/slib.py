"""Shared helpers for the additional-planet / rotation / flare search on TIC 214311908 (M3V).

Inputs are read-only from ../phot and ../referee_phot. Nothing outside ../search is written.
Conventions: times in BTJD_TDB (BJD_TDB - 2457000), target-barycentred (phot pipeline).
Fractional fluxes are rescaled to the PRF-unit M3V flux of 979 e-/s (referee correction): f = 1 + 0.955 (flux_norm - 1).
"""
import os, json, warnings
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2"); os.environ.setdefault("NUMBA_NUM_THREADS", "2")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import numba as nb

BASE = "tic214311908/analysis"
ROOT = f"{BASE}/search"
PHOT = f"{BASE}/phot"
SECTORS = (12, 39, 66, 93, 103, 104)
SCALE = 935.97 / 979.0            # PRF-unit rescaling of fractional fluxes (flux_norm uses tmag_flux(13.0135) = 935.97 e-/s)
FREF = 979.0
# referee-adopted ephemeris (phot/out/transits.json global_B)
T0_B, P_B = 3840.914139360305, 5.691508695210822
T14_B = 0.0764                       # longest of the quoted T14 values (1.83 h)
MASK_HALF = 0.1                      # +-0.1 d around each 5.69-d transit (= phot in_transit_window)
RSTAR, MSTAR, RHO = 0.381, 0.370, 9.4          # star/star_params.json (single-star case)
REARTH_RSUN = 0.0091577


def save_json(obj, fn):
    def conv(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.bool_,)): return bool(o)
        raise TypeError(type(o))
    with open(fn, "w") as f: json.dump(obj, f, indent=1, default=conv)


def robust_std(x):
    x = np.asarray(x); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def transit_mask(t, P=P_B, T0=T0_B, half=MASK_HALF):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    return np.abs(ph * P) < half


def load_star(star="M3V", model="B"):
    """Full-resolution light curve of one deblended star from ../phot (quality_ok only), rescaled.
    Returns DataFrame with t, sector, cad (s), raw (undetrended, normalised), det (detrended), err, intr."""
    d = pd.read_csv(f"{PHOT}/lc_{star}_all.csv")
    d = d[d.quality_ok == 1].copy()
    if star == "M3V":
        if model == "B":
            raw, det, err = d.flux_raw_norm, d.flux_norm, d.flux_err
        else:
            raw, det, err = d.flux_raw_norm_modelA, d.flux_norm_modelA, d.flux_err_modelA
        raw = 1 + SCALE * (raw - 1); det = 1 + SCALE * (det - 1); err = SCALE * err
    else:
        if model in ("A",):
            raw, det = d.flux_raw_norm_modelA, d.flux_norm_modelA
        else:
            raw, det = d.flux_raw_norm, d.flux_norm
        err = pd.Series(np.full(len(d), np.nan), index=d.index)
    out = pd.DataFrame(dict(t=d.time_btjd_tdb.values, sector=d.sector.values, cad=d.cadence_s.values,
                            raw=np.asarray(raw, float), det=np.asarray(det, float), err=np.asarray(err, float)))
    out["intr"] = transit_mask(out.t.values)
    out = out[np.isfinite(out.det)].reset_index(drop=True)
    return out


def bin_lc(t, f, sector, cad, width_min=10.0):
    """Bin within each sector to max(width, cadence). Returns t_b, f_b, n_b, sector_b, cad_b (s of bin)."""
    tb, fb, nb_, sb, cb = [], [], [], [], []
    for s in np.unique(sector):
        m = sector == s
        w = max(width_min * 60.0, float(np.median(cad[m]))) / 86400.0
        tt, ff = t[m], f[m]
        k = np.floor((tt - tt.min()) / w + 1e-6).astype(int)
        u, inv = np.unique(k, return_inverse=True)
        n = np.bincount(inv)
        tb.append(np.bincount(inv, tt) / n); fb.append(np.bincount(inv, ff) / n); nb_.append(n)
        sb.append(np.full(len(u), s)); cb.append(np.full(len(u), w * 86400.0))
    return (np.concatenate(tb), np.concatenate(fb), np.concatenate(nb_), np.concatenate(sb), np.concatenate(cb))


def detrend(t, f, mask, window=0.75):
    """wotan biweight as in phot/p06 (break tolerance 0.5 d, edge_cutoff 0), masked points interpolated."""
    from wotan import flatten
    tr = np.full(len(t), np.nan)
    m = np.isfinite(f) & ~mask
    trm = flatten(t[m], f[m], method="biweight", window_length=window, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
    k = np.isfinite(trm)
    tr[:] = np.interp(t, t[m][k], trm[k])
    return tr


# ----------------------------------------------------------------------------------- transit duration model
def a_over_r(P, rho=RHO):
    G = 6.674e-8
    return (G * rho * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1 / 3)


def t14_central(P, rho=RHO, k=0.0):
    aR = a_over_r(P, rho)
    return P / np.pi * np.arcsin(np.clip((1 + k) / aR, 0, 1))


DUR_FACS = np.array([0.3, 0.45, 0.65, 0.9, 1.25])      # x central-transit duration


def freq_grid(pmin, pmax, baseline, os_=2.0, dfac=0.45, rho=RHO):
    """Frequencies such that the phase drift over the baseline is < (dfac * T14c(P)) / os_ (Ofir-like, duration-aware)."""
    f = 1.0 / pmax; out = [f]
    fmax = 1.0 / pmin
    while f < fmax:
        P = 1.0 / f
        D = max(dfac * t14_central(P, rho), 0.4 / 24)
        f = f + D * f / (baseline * os_)
        out.append(f)
    return np.array(out)


# ----------------------------------------------------------------------------------- numba BLS
@nb.njit(cache=True, fastmath=True)
def _bls_one(t, y, w, freq, durs, nbin_per_dmin, tref):
    """Box search at one frequency. y = flux - 1 (mean-subtracted with weights), w = 1/sigma^2.
    Returns best dchi2 (dips only), depth, phase of mid-transit, duration index."""
    P = 1.0 / freq
    dmin = durs[0]
    nb_ = int(np.ceil(P / (dmin / nbin_per_dmin)))
    if nb_ < 8: nb_ = 8
    sw = np.zeros(nb_); swy = np.zeros(nb_)
    for i in range(t.shape[0]):
        ph = (t[i] - tref) * freq
        ph = ph - np.floor(ph)
        b = int(ph * nb_)
        if b >= nb_: b = nb_ - 1
        sw[b] += w[i]; swy[b] += w[i] * y[i]
    W = sw.sum()
    best = 0.0; bd = 0.0; bph = 0.0; bk = -1
    for k in range(durs.shape[0]):
        m = int(np.round(durs[k] / P * nb_))
        if m < 1: m = 1
        if m >= nb_ // 2: continue
        # running window sums over circular bins
        a = 0.0; c = 0.0
        for j in range(m):
            a += sw[j]; c += swy[j]
        for j in range(nb_):
            if a > 0 and a < W and c < 0:
                d2 = c * c / (a * (1.0 - a / W))
                if d2 > best:
                    best = d2; bd = -c / a / (1.0 - a / W); bph = (j + 0.5 * m) / nb_; bk = k
            # slide
            jo = j; jn = (j + m) % nb_
            a += sw[jn] - sw[jo]; c += swy[jn] - swy[jo]
    return best, bd, bph, bk


@nb.njit(parallel=True, cache=True)
def bls_grid(t, y, w, freqs, durfac, rho, nbin_per_dmin, tref):
    n = freqs.shape[0]
    D2 = np.zeros(n); DEP = np.zeros(n); PH = np.zeros(n); K = np.zeros(n, dtype=np.int64)
    G = 6.674e-8
    for i in nb.prange(n):
        P = 1.0 / freqs[i]
        aR = (G * rho * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)
        x = 1.0 / aR
        if x > 1: x = 1.0
        tc = P / np.pi * np.arcsin(x)
        durs = durfac * tc
        for k in range(durs.shape[0]):
            if durs[k] < 0.35 / 24: durs[k] = 0.35 / 24
        d2, dep, ph, kk = _bls_one(t, y, w, freqs[i], durs, nbin_per_dmin, tref)
        D2[i] = d2; DEP[i] = dep; PH[i] = ph; K[i] = kk
    return D2, DEP, PH, K


def prep_weights(y, sector):
    """Per-sector white-noise weights from the robust std of the (binned) data; y mean-subtracted per sector."""
    y = y.copy(); w = np.zeros(len(y))
    for s in np.unique(sector):
        m = sector == s
        y[m] -= np.median(y[m])
        w[m] = 1.0 / robust_std(y[m]) ** 2
    y -= np.sum(w * y) / np.sum(w)
    return y, w


def run_bls(t, y, w, freqs, durfac=DUR_FACS, rho=RHO, nbin=3, tref=None):
    tref = float(np.min(t)) if tref is None else tref
    D2, DEP, PH, K = bls_grid(t.astype(np.float64), y.astype(np.float64), w.astype(np.float64), freqs.astype(np.float64),
                              durfac.astype(np.float64), float(rho), int(nbin), float(tref))
    T0 = tref + PH / freqs
    return dict(freq=freqs, dchi2=D2, depth=DEP, t0=T0, kdur=K)


def sde(power):
    return (power - np.mean(power)) / np.std(power)


# ----------------------------------------------------------------------------------- transit model for injections
def batman_model(t, t0, P, rp_rs, b, cad_s, u=(0.25, 0.35), rho=RHO):
    import batman
    pm = batman.TransitParams()
    pm.t0 = t0; pm.per = P; pm.rp = rp_rs; pm.a = a_over_r(P, rho)
    inc = np.degrees(np.arccos(np.clip(b / pm.a, 0, 1)))
    pm.inc = inc; pm.ecc = 0.; pm.w = 90.; pm.u = list(u); pm.limb_dark = "quadratic"
    out = np.ones(len(t))
    for c in np.unique(cad_s):
        m = cad_s == c
        sf = max(1, int(round(c / 60.0)))
        mdl = batman.TransitModel(pm, t[m], supersample_factor=sf, exp_time=c / 86400.0)
        out[m] = mdl.light_curve(pm)
    return out
