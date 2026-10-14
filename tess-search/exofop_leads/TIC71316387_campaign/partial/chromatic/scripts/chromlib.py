"""Shared helpers for the chromatic eclipse-depth study of TIC 71316387.

Conventions
-----------
* Times are BJD_TDB (full JD) unless stated otherwise.
* Ephemeris (deep eclipse): T0 = 2460165.6697 BJD_TDB, P = 1.3036712 +- 0.0000005 d.
* Eclipse shape: unit-depth "generalised V" template s(x) = (1-|x|^p)^q, |x|<1,
  x = (t - tc)/(T/2), integrated over the exposure time.  Shape parameters are fitted to
  short-cadence TESS data (tess_ref.py) and frozen for the ground surveys.
* Linear model for any light curve (relative flux, normalised to ~1):
      f = sum_seg [c0_seg + c1_seg*(t - tmid_seg) (+ c2_seg*(t-tmid)^2)] - A_d*s_d(t) - A_s*s_s(t)
  i.e. per-segment (night) polynomial baselines fitted SIMULTANEOUSLY with the eclipse
  amplitudes (no masking bias).  A_d, A_s are the fractional depths at eclipse centre.
"""
import numpy as np

P = 1.3036712
P_ERR = 0.0000005
T0 = 2460165.6697
RA, DEC = 330.887948, -42.835511


def phase(t, t0=T0, p=P):
    return ((t - t0) / p) % 1.0


def vshape(x, p, q):
    ax = np.abs(x)
    out = np.zeros_like(ax)
    m = ax < 1
    out[m] = (1.0 - ax[m] ** p) ** q
    return out


def template(t, tc_phase, T, p, q, exptime=0.0, nsub=1, t0=T0, per=P):
    """Unit-depth eclipse template centred at phase tc_phase (cycles), full width T (days).
    Exposure-time integration via nsub sub-samples."""
    t = np.asarray(t, float)
    if exptime > 0 and nsub > 1:
        off = (np.arange(nsub) + 0.5) / nsub - 0.5
        tt = t[:, None] + off[None, :] * exptime
    else:
        tt = t[:, None]
    ph = ((tt - t0) / per - tc_phase) % 1.0
    ph[ph > 0.5] -= 1.0
    x = ph * per / (T / 2.0)
    return vshape(x, p, q).mean(axis=1)


def segments(t, gap=0.3, extra=None):
    """Split into contiguous segments (nights); `extra` (e.g. camera id) also splits."""
    idx = np.argsort(t)
    seg = np.zeros(len(t), int)
    ts = t[idx]
    brk = np.r_[0, np.cumsum(np.diff(ts) > gap)]
    seg[idx] = brk
    if extra is not None:
        _, ex = np.unique(extra, return_inverse=True)
        seg = seg * (ex.max() + 1) + ex
        _, seg = np.unique(seg, return_inverse=True)
    return seg


def design(t, seg, templates, order=1):
    """Design matrix: per-segment polynomial columns + negative template columns."""
    nseg = seg.max() + 1
    cols = []
    for k in range(nseg):
        m = seg == k
        tm = np.median(t[m])
        for o in range(order + 1):
            c = np.zeros(len(t))
            c[m] = (t[m] - tm) ** o
            cols.append(c)
    nb = len(cols)
    for s in templates:
        cols.append(-s)
    return np.vstack(cols).T, nb


def lsq(t, f, w, seg, templates, order=1):
    """Weighted linear least squares with per-segment polynomial baselines profiled out
    analytically (block elimination / Schur complement); exact same solution as the dense
    design-matrix fit (lsq_dense) but O(N) instead of O(N*nseg^2).
    Returns amplitudes (template coefficients, model = baseline - sum a_i s_i),
    formal errors (assuming weights = 1/sigma^2), residuals, model."""
    t = np.asarray(t, float); f = np.asarray(f, float); w = np.asarray(w, float)
    seg = np.asarray(seg)
    _, sg = np.unique(seg, return_inverse=True)
    nseg = sg.max() + 1
    tm = np.zeros(nseg)
    # segment reference time (median is slow per segment; mean is equivalent for the fit)
    tm = np.bincount(sg, weights=t, minlength=nseg) / np.maximum(np.bincount(sg, minlength=nseg), 1)
    x = t - tm[sg]
    no = order + 1
    B = np.vstack([x ** o for o in range(no)]).T            # N x no
    M = np.zeros((nseg, no, no))
    for i in range(no):
        for j in range(i, no):
            v = np.bincount(sg, weights=w * B[:, i] * B[:, j], minlength=nseg)
            M[:, i, j] = v; M[:, j, i] = v
    Minv = np.linalg.pinv(M)

    def proj(v):
        V = np.vstack([np.bincount(sg, weights=w * B[:, i] * v, minlength=nseg) for i in range(no)]).T  # nseg x no
        c = np.einsum('kij,kj->ki', Minv, V)
        return v - np.sum(B * c[sg], axis=1), c
    fp, cf = proj(f)
    S = [-np.asarray(s, float) for s in templates]   # model uses -A*s
    Sp = [proj(s)[0] for s in S]
    m = len(S)
    G = np.array([[np.sum(w * Sp[i] * Sp[j]) for j in range(m)] for i in range(m)])
    b = np.array([np.sum(w * Sp[i] * fp) for i in range(m)])
    try:
        Gi = np.linalg.inv(G)
    except np.linalg.LinAlgError:
        Gi = np.linalg.pinv(G)
    a = Gi @ b
    res = fp - sum(a[i] * Sp[i] for i in range(m)) if m else fp
    model = f - res
    err = np.sqrt(np.abs(np.diag(Gi)))
    return a, err, res, model


def lsq_dense(t, f, w, seg, templates, order=1):
    """Reference dense implementation (slow); used only to validate lsq()."""
    X, nb = design(t, seg, templates, order)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    fw = f * sw
    # drop empty columns (segments too short for the order)
    good = np.any(Xw != 0, axis=0)
    Xg = Xw[:, good]
    beta, *_ = np.linalg.lstsq(Xg, fw, rcond=None)
    full = np.zeros(X.shape[1])
    full[good] = beta
    model = X @ full
    res = f - model
    try:
        cov = np.linalg.pinv(Xg.T @ Xg)
        err_g = np.sqrt(np.diag(cov))
        err = np.full(X.shape[1], np.nan)
        err[good] = err_g
    except np.linalg.LinAlgError:
        err = np.full(X.shape[1], np.nan)
    return full[nb:], err[nb:], res, model


def mad_std(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def hjd_utc_to_bjd_tdb(hjd, ra=RA, dec=DEC):
    """Convert HJD(UTC) to BJD_TDB for a geocentric observer (error < 0.1 s for ground sites)."""
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    c = SkyCoord(ra * u.deg, dec * u.deg)
    loc = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    t = Time(np.asarray(hjd), format='jd', scale='utc', location=loc)
    # iterate: geocentric JD = HJD - helio correction(JD)
    jd = t
    for _ in range(3):
        lt_h = jd.light_travel_time(c, kind='heliocentric')
        jd = Time(np.asarray(hjd) - lt_h.to_value(u.day), format='jd', scale='utc', location=loc)
    lt_b = jd.light_travel_time(c, kind='barycentric')
    return (jd.tdb + lt_b).jd


def jd_utc_geo_to_bjd_tdb(jd_utc, ra=RA, dec=DEC):
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    c = SkyCoord(ra * u.deg, dec * u.deg)
    loc = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    t = Time(np.asarray(jd_utc), format='jd', scale='utc', location=loc)
    return (t.tdb + t.light_travel_time(c, kind='barycentric')).jd


def fake_phases(n, rng, avoid=(0.0, 0.5), halfwidth=0.08):
    """Random phase offsets dphi such that fake deep (dphi) and fake shallow (0.5+dphi)
    windows avoid the real eclipses by > halfwidth in phase."""
    out = []
    while len(out) < n:
        d = rng.uniform(0, 1)
        ok = True
        for a in avoid:
            for c in (d, (d + 0.5) % 1):
                dd = abs(((c - a) + 0.5) % 1 - 0.5)
                if dd < halfwidth:
                    ok = False
        if ok:
            out.append(d)
    return np.array(out)
