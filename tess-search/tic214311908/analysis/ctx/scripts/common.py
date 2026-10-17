"""Shared helpers for the ctx investigation of TIC 214311908.

Ephemeris, transit template (trapezoid integrated over the exposure), detrending with the transits masked,
linear depth fits, and an empirical noise calibration from depths measured at random (signal-free) phases.
"""
import json
import numpy as np
from wotan import flatten

EPH = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))
P = EPH["P"]; T0 = EPH["T0"]; P_ERR = EPH["P_err"]; T0_ERR = EPH["T0_err"]
QLP_SHIFT = {int(k): v for k, v in EPH["qlp_shift_days"].items()}   # add to QLP S103/S104 TIME
T14 = 1.63 / 24.0          # d, total duration (earlier fits 1.61-1.65 h)
FLAT = 0.70                # T23/T14 (earlier: flat part 60-80% of the transit)


def trap(t, tc, depth=1.0, t14=T14, flat=FLAT):
    """Trapezoid of unit depth (returns fractional dimming, positive in transit)."""
    t23 = flat * t14
    x = np.abs(t - tc)
    y = np.clip((t14 / 2 - x) / ((t14 - t23) / 2), 0, 1)
    return depth * y


def template(t, exptime_d, per=P, t0=T0, phase0=0.0, t14=T14, flat=FLAT, nsub=15):
    """Fractional dimming template (unit depth) for an array of mid-exposure times, integrated over exptime."""
    off = (np.arange(nsub) + 0.5) / nsub - 0.5
    tt = t[:, None] + off[None, :] * exptime_d
    ph = ((tt - t0) / per - phase0 + 0.5) % 1 - 0.5
    return trap(ph * per, 0.0, 1.0, t14, flat).mean(axis=1)


def epoch(t, per=P, t0=T0, phase0=0.0):
    return np.round((t - t0) / per - phase0).astype(int)


def detrend(t, f, mask, window=1.0, break_tol=0.4):
    """Biweight trend fitted on points with mask==False, interpolated across masked points (same as m3v_checks.py)."""
    ok = ~mask & np.isfinite(f)
    tr = flatten(t[ok], f[ok], method="biweight", window_length=window, return_trend=True, break_tolerance=break_tol)[1]
    k = np.isfinite(tr)
    trend = np.interp(t, t[ok][k], tr[k])
    return f / trend, trend


def in_window(t, per=P, t0=T0, phase0=0.0, half=1.0 * T14):
    ph = ((t - t0) / per - phase0 + 0.5) % 1 - 0.5
    return np.abs(ph * per) < half


def robust_std(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def fit_depth(f, s, tmpl):
    """Weighted LSQ depth for f = 1 - d*tmpl (f already normalised). Returns depth, error."""
    w = 1.0 / s ** 2
    a = np.sum(w * tmpl * tmpl)
    if a <= 0:
        return np.nan, np.nan
    d = np.sum(w * tmpl * (1 - f)) / a
    return d, 1 / np.sqrt(a)


def random_phase_depths(t, f, s, exptime_d, n=400, per=P, t0=T0, avoid=(0.0, 0.5), rng=None, min_pts=3):
    """Depths measured at random phases (with the same period) away from phase 0 and 0.5: empirical noise of
    the depth estimator including red noise and detrending effects."""
    rng = np.random.default_rng(1) if rng is None else rng
    out = []
    dph = 2.5 * T14 / per
    tries = 0
    while len(out) < n and tries < 20 * n:
        tries += 1
        ph0 = rng.uniform(0, 1)
        if any(min(abs(ph0 - a), 1 - abs(ph0 - a)) < dph for a in avoid):
            continue
        tm = template(t, exptime_d, per, t0, ph0)
        if np.sum(tm > 0.5) < min_pts:
            continue
        out.append(fit_depth(f, s, tm)[0])
    return np.array(out)


def detrend_additive(t, f, mask, window=1.0, break_tol=0.4):
    """Additive version of detrend(): returns the trend in the units of f. wotan's biweight is not offset-invariant
    when a light curve crosses zero (found with TGLC pairB: 35.7 vs 13.5 e/s), so a large constant is added first."""
    off = 10 * np.nanmax(np.abs(f)) + 1.0
    _, tr = detrend(t, (f + off) / off, mask, window, break_tol)
    return tr * off - off
