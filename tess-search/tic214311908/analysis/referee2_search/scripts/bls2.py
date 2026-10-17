"""Referee-2 box search (independent implementation; not the analyst's slib).

Statistic: for a box of duration D at period P and epoch T, with per-point weights w = 1/sigma^2 and y mean-subtracted
per sector, dchi2 = S^2 / (W_in (1 - W_in / W)), S = -sum_in w y (dips only, S > 0).
Phase handling: every point is assigned to fine phase bins of width D_min/5 at each trial frequency; boxes are
contiguous runs of fine bins (cumulative sums), all start bins tested.
Durations: 6 values log-spaced in 0.35-1.4 x T14c(P) for rho = 9.4 g/cm3 (covers the twin case rho = 15.5 and
b up to ~0.9), floored at 0.5 h.
Frequency grid: df = D_min(P) / (3 B P), i.e. the box drifts by < D_min/3 over the baseline B.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "2")
import numpy as np
import numba as nb

G = 6.674e-8


def t14c(P, rho=9.4):
    aR = (G * rho * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1 / 3)
    return P / np.pi * np.arcsin(np.minimum(1.0, 1 / aR))


DFAC = np.exp(np.linspace(np.log(0.35), np.log(1.4), 6))
DMIN_D = 0.5 / 24


def durations(P, rho=9.4):
    d = DFAC * t14c(P, rho)
    return np.maximum(d, DMIN_D)


def freq_grid(pmin, pmax, B, rho=9.4, over=3.0):
    f = 1.0 / pmax; out = [f]
    while f < 1.0 / pmin:
        Pc = 1.0 / f
        dmin = max(DFAC[0] * t14c(Pc, rho), DMIN_D)
        f += dmin / (over * B * Pc)
        out.append(f)
    return np.array(out)


@nb.njit(cache=True)
def _one(t, y, w, f, dfac, dmin_floor, rho, tref):
    P = 1.0 / f
    aR = (6.674e-8 * rho * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)
    x = 1.0 / aR
    if x > 1.0: x = 1.0
    tc = P / np.pi * np.arcsin(x)
    nd = dfac.shape[0]
    durs = np.empty(nd)
    for k in range(nd):
        durs[k] = max(dfac[k] * tc, dmin_floor)
    binw = durs[0] / 5.0
    nb_ = int(P / binw) + 1
    sw = np.zeros(nb_); swy = np.zeros(nb_)
    for i in range(t.shape[0]):
        ph = (t[i] - tref) * f
        ph -= np.floor(ph)
        b = int(ph * nb_)
        if b >= nb_: b = nb_ - 1
        sw[b] += w[i]; swy[b] += w[i] * y[i]
    W = 0.0
    for b in range(nb_): W += sw[b]
    # cumulative sums over a doubled array (circular)
    cs = np.zeros(2 * nb_ + 1); cy = np.zeros(2 * nb_ + 1)
    for b in range(2 * nb_):
        cs[b + 1] = cs[b] + sw[b % nb_]; cy[b + 1] = cy[b] + swy[b % nb_]
    best = 0.0; bdep = 0.0; bph = 0.0; bk = -1
    for k in range(nd):
        m = int(round(durs[k] / P * nb_))
        if m < 1: m = 1
        if m > nb_ // 2: continue
        for j in range(nb_):
            a = cs[j + m] - cs[j]
            c = cy[j + m] - cy[j]
            if a <= 0.0 or c >= 0.0 or a >= W: continue
            d2 = c * c / (a * (1.0 - a / W))
            if d2 > best:
                best = d2; bdep = -c / a / (1.0 - a / W); bph = (j + 0.5 * m) / nb_; bk = k
    return best, bdep, bph, bk


@nb.njit(parallel=True, cache=True)
def bls(t, y, w, freqs, dfac, dmin_floor, rho, tref):
    n = freqs.shape[0]
    D2 = np.zeros(n); DEP = np.zeros(n); PH = np.zeros(n); K = np.zeros(n, dtype=np.int64)
    for i in nb.prange(n):
        a, b, c, d = _one(t, y, w, freqs[i], dfac, dmin_floor, rho, tref)
        D2[i] = a; DEP[i] = b; PH[i] = c; K[i] = d
    return D2, DEP, PH, K


def run(t, y, w, freqs, rho=9.4):
    tref = float(t.min())
    D2, DEP, PH, K = bls(t.astype(np.float64), y.astype(np.float64), w.astype(np.float64), freqs.astype(np.float64),
                         DFAC.astype(np.float64), float(DMIN_D), float(rho), tref)
    return dict(freq=freqs, dchi2=D2, depth=DEP, t0=tref + PH / freqs, k=K)


def weights(y, s):
    """per-sector median subtraction and robust white weights."""
    y = y.copy(); w = np.zeros(len(y))
    for q in np.unique(s):
        m = s == q
        y[m] -= np.median(y[m])
        r = 1.4826 * np.median(np.abs(y[m] - np.median(y[m])))
        w[m] = 1 / r ** 2
    y -= np.sum(w * y) / np.sum(w)
    return y, w
