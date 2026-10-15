"""lib3: core library for the v3 chromatic eclipse-depth study of TIC 71316387.

Independent re-implementation (v3).  Conventions
  * times BJD_TDB (full JD)
  * deep-eclipse ephemeris T0 = 2460165.6697, P = 1.3036712 d (re-checked in ephem3.py)
  * eclipse template: unit-depth generalised V  s(x) = (1-|x|^p)^q for |x|<1, x=(t-tc)/(T/2),
    integrated over the exposure.  "Depth" A = amplitude of the template = central fractional depth.
  * light-curve model (relative flux):  f = sum_k poly_k(t) - A_d s_d(t) - A_s s_s(t)
    poly_k = polynomial baseline per segment (night x camera), fitted SIMULTANEOUSLY with A.
"""
import numpy as np

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad'
CH = BASE + '/tic713/chromatic'
V3 = CH + '/v3'
P = 1.3036712
P_ERR = 0.0000005
T0 = 2460165.6697
RA, DEC = 330.887948, -42.835511


# ----------------------------------------------------------------------------- time systems
def _sky():
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    return SkyCoord(RA * u.deg, DEC * u.deg)


def hjd_utc_to_bjd_tdb(hjd):
    """HJD(UTC) -> BJD_TDB for a geocentric observer (iterative; < 0.1 s for ground sites)."""
    from astropy.time import Time
    from astropy.coordinates import EarthLocation
    import astropy.units as u
    c = _sky()
    loc = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    hjd = np.asarray(hjd, float)
    jd = Time(hjd, format='jd', scale='utc', location=loc)
    for _ in range(3):
        lt = jd.light_travel_time(c, kind='heliocentric').to_value(u.day)
        jd = Time(hjd - lt, format='jd', scale='utc', location=loc)
    return (jd.tdb + jd.light_travel_time(c, kind='barycentric')).jd


def jd_utc_to_bjd_tdb(jd_utc):
    from astropy.time import Time
    from astropy.coordinates import EarthLocation
    import astropy.units as u
    loc = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    t = Time(np.asarray(jd_utc, float), format='jd', scale='utc', location=loc)
    return (t.tdb + t.light_travel_time(_sky(), kind='barycentric')).jd


# ----------------------------------------------------------------------------- template
def vshape(x, p, q):
    ax = np.abs(x)
    out = np.zeros_like(ax, dtype=float)
    m = ax < 1
    out[m] = (1.0 - ax[m] ** p) ** q
    return out


def template(t, phc, T, p, q, exptime=0.0, t0=T0, per=P):
    """unit-depth template centred at phase phc; exposure-integrated with ~60 s sub-samples."""
    t = np.asarray(t, float)
    nsub = max(1, int(np.ceil(exptime / (60.0 / 86400.0)))) if exptime > 0 else 1
    off = ((np.arange(nsub) + 0.5) / nsub - 0.5) * exptime if nsub > 1 else np.zeros(1)
    out = np.zeros(len(t))
    for o in off:
        ph = ((t + o - t0) / per - phc) % 1.0
        ph[ph > 0.5] -= 1.0
        out += vshape(ph * per / (T / 2.0), p, q)
    return out / len(off)


def tmpl(t, shp, which, dphi=0.0, exptime=0.0, per=None, t0=None):
    """template for eclipse `which` ('d'/'s'); a shape dict may carry its own ephemeris ('per', 't0')"""
    s = shp[which]
    per = s.get('per', P) if per is None else per
    t0 = s.get('t0', T0) if t0 is None else t0
    return template(t, s['phc'] + dphi, s['T'], s['p'], s['q'], exptime, t0=t0, per=per)


def survey_shapes(shp, ephem):
    """copy of the TESS template shapes placed on the refined linear ephemeris (ephem3.json 'joint'):
    deep at T0 + E*P' + oc_d, shallow at T0 + (E+0.5)*P' + oc_s"""
    j = ephem['joint']
    out = {}
    for k, oc in (('d', j['oc0_deep']), ('s', j['oc0_shallow'])):
        out[k] = dict(shp[k], per=j['P'], t0=T0, phc=(0.0 if k == 'd' else 0.5) + oc / 1440.0 / j['P'])
    return out


def dphase(t, phc, per=P, t0=T0):
    """signed time (d) from nearest occurrence of phase phc"""
    ph = ((t - t0) / per - phc + 0.5) % 1.0 - 0.5
    return ph * per


# ----------------------------------------------------------------------------- segments
def segments(t, gap=0.25, extra=None):
    o = np.argsort(t)
    brk = np.r_[0, np.cumsum(np.diff(t[o]) > gap)]
    seg = np.empty(len(t), int); seg[o] = brk
    if extra is not None:
        _, ex = np.unique(extra, return_inverse=True)
        seg = seg * (ex.max() + 1) + ex
    _, seg = np.unique(seg, return_inverse=True)
    return seg


def mad_std(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


# ----------------------------------------------------------------------------- linear fit
class Projector:
    """Weighted projection that removes a per-segment polynomial (order `order`) from vectors.
    v_perp = v - B (B^T W B)^-1 B^T W v, block-diagonal over segments.  Exact."""

    def __init__(self, t, w, seg, order=1):
        _, self.sg = np.unique(seg, return_inverse=True)
        self.n = self.sg.max() + 1
        self.w = np.asarray(w, float)
        cnt = np.bincount(self.sg, minlength=self.n)
        tm = np.bincount(self.sg, weights=t, minlength=self.n) / np.maximum(cnt, 1)
        x = t - tm[self.sg]
        sc = np.sqrt(np.bincount(self.sg, weights=x * x, minlength=self.n) / np.maximum(cnt, 1)) + 1e-12
        x = x / sc[self.sg]                     # scaled for conditioning
        self.no = order + 1
        self.B = np.vstack([x ** o for o in range(self.no)]).T
        M = np.zeros((self.n, self.no, self.no))
        for i in range(self.no):
            for j in range(i, self.no):
                v = np.bincount(self.sg, weights=self.w * self.B[:, i] * self.B[:, j], minlength=self.n)
                M[:, i, j] = v; M[:, j, i] = v
        self.Minv = np.linalg.pinv(M)

    def __call__(self, v):
        V = np.vstack([np.bincount(self.sg, weights=self.w * self.B[:, i] * v, minlength=self.n)
                       for i in range(self.no)]).T
        c = np.einsum('kij,kj->ki', self.Minv, V)
        return v - np.sum(self.B * c[self.sg], axis=1)


def fit_amps(t, f, w, seg, cols, order=1, proj=None):
    """Least-squares amplitudes a for model f = baseline - sum a_i cols_i.
    Returns a, formal errors (if w = 1/sigma^2), residuals, projector."""
    if proj is None:
        proj = Projector(t, w, seg, order)
    fp = proj(np.asarray(f, float))
    S = [proj(-np.asarray(c, float)) for c in cols]
    m = len(S)
    G = np.array([[np.sum(proj.w * S[i] * S[j]) for j in range(m)] for i in range(m)])
    b = np.array([np.sum(proj.w * S[i] * fp) for i in range(m)])
    Gi = np.linalg.pinv(G)
    a = Gi @ b
    res = fp - sum(a[i] * S[i] for i in range(m))
    return a, np.sqrt(np.abs(np.diag(Gi))), res, proj


def fit_dense(t, f, w, seg, cols, order=1):
    """reference dense solver (validation only)"""
    _, sg = np.unique(seg, return_inverse=True)
    X = []
    for k in range(sg.max() + 1):
        m = sg == k
        tm = t[m].mean()
        for o in range(order + 1):
            c = np.zeros(len(t)); c[m] = (t[m] - tm) ** o; X.append(c)
    nb = len(X)
    X += [-np.asarray(c) for c in cols]
    X = np.array(X).T
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], f * sw, rcond=None)
    return beta[nb:]


# ----------------------------------------------------------------------------- fake phases
def fake_offsets(n, rng, guard=0.08):
    """phase offsets d such that fake eclipses at d and 0.5+d lie > guard (phase) from the real
    eclipses at 0 and 0.5 (guard 0.08 = 2.5 h vs eclipse half-width ~0.8 h)."""
    out = []
    while len(out) < n:
        d = rng.uniform(0, 1)
        dd = abs(((d % 0.5) + 0.25) % 0.5 - 0.25)    # distance of d from {0, 0.5} mod 0.5
        if dd > guard:
            out.append(d)
    return np.array(out)


def coverage(t, shp, which, dphi=0.0, exptime=0.0, thr=0.5):
    """number of points and distinct events with template > thr"""
    s = tmpl(t, shp, which, dphi, exptime)
    m = s > thr
    phc = shp[which]['phc'] + dphi
    ev = np.unique(np.round((t[m] - shp[which].get('t0', T0)) / shp[which].get('per', P) - phc))
    return int(m.sum()), int(len(ev))


# ----------------------------------------------------------------------------- full survey measurement
def measure_at(t, f, w, proj, shp, exptime, dphi):
    cols = [tmpl(t, shp, 'd', dphi, exptime), tmpl(t, shp, 's', dphi, exptime)]
    if dphi != 0.0:
        cols += [tmpl(t, shp, 'd', 0.0, exptime), tmpl(t, shp, 's', 0.0, exptime)]
    if np.sum(cols[0] > 0.5) < 3 or np.sum(cols[1] > 0.5) < 3:
        return None
    a, e, r, _ = fit_amps(t, f, w, None, cols, proj=proj)
    return a[:2], e[:2]


def clean_and_weight(t, f, seg, shp, exptime, order=1, nclip=4.0, niter=8, minpts=10, floor=2e-4, extra_dphi=()):
    """iterative sigma clipping on residuals of the full (baseline + real eclipse) model;
    per-segment robust sigma -> weights (with a floor).  Returns mask, segment ids (re-indexed), weights.
    extra_dphi: additional eclipse-pair phases included in the cleaning model (used by end-to-end
    injections so injected eclipses are treated exactly like the real ones)."""
    keep = np.ones(len(t), bool)
    cols_all = [tmpl(t, shp, 'd', 0, exptime), tmpl(t, shp, 's', 0, exptime)]
    for d in extra_dphi:
        cols_all += [tmpl(t, shp, 'd', d, exptime), tmpl(t, shp, 's', d, exptime)]
    for it in range(niter):
        cnt = np.bincount(seg[keep], minlength=seg.max() + 1)
        keep &= cnt[seg] >= minpts
        tt, ff, sg = t[keep], f[keep], seg[keep]
        a, e, r, pr = fit_amps(tt, ff, np.ones(len(tt)), sg, [c[keep] for c in cols_all], order)
        sig = np.zeros(len(tt))
        for k in np.unique(sg):
            m = sg == k
            sig[m] = max(mad_std(r[m]), floor)
        bad = np.abs(r) > nclip * sig
        if not bad.any():
            break
        idx = np.where(keep)[0]
        keep[idx[bad]] = False
    cnt = np.bincount(seg[keep], minlength=seg.max() + 1)
    keep &= cnt[seg] >= minpts
    tt, ff, sg = t[keep], f[keep], seg[keep]
    a, e, r, pr = fit_amps(tt, ff, np.ones(len(tt)), sg, [c[keep] for c in cols_all], order)
    _, sg2 = np.unique(sg, return_inverse=True)
    w = np.zeros(len(tt))
    for k in np.unique(sg2):
        m = sg2 == k
        w[m] = 1.0 / max(mad_std(r[m]), floor) ** 2
    return keep, sg2, w


def run_survey(t, f, seg, shp, exptime, order=1, nfake=500, ninj=150, inj_depths=(0.003, 0.006),
               s_over_d=0.464, seed=11, nclip=4.0, minpts=10, floor=2e-4, verbose=False):
    """Full measurement: clean, fit real eclipses, fake-phase null (real eclipses modelled),
    injection-recovery at fake phases.  Returns dict and cleaned arrays."""
    rng = np.random.default_rng(seed)
    keep, sg, w = clean_and_weight(t, f, seg, shp, exptime, order, nclip=nclip, minpts=minpts, floor=floor)
    t, f = t[keep], f[keep]
    proj = Projector(t, w, sg, order)
    r0 = measure_at(t, f, w, proj, shp, exptime, 0.0)
    if r0 is None:
        return None, None
    a, e = r0
    nd, ed = coverage(t, shp, 'd', 0, exptime)
    ns, es = coverage(t, shp, 's', 0, exptime)
    res = dict(n=int(len(t)), nseg=int(sg.max() + 1), order=order, nclip=nclip,
               n_in_deep=nd, ev_deep=ed, n_in_shallow=ns, ev_shallow=es,
               rms_med=float(np.median(1 / np.sqrt(w))),
               deep=float(a[0]), shallow=float(a[1]), deep_formal=float(e[0]), shallow_formal=float(e[1]))
    fk, fe, offs = [], [], []
    for d in fake_offsets(nfake, rng):
        r_ = measure_at(t, f, w, proj, shp, exptime, d)
        if r_ is None:
            continue
        fk.append(r_[0]); fe.append(r_[1]); offs.append(d)
    fk, fe = np.array(fk), np.array(fe)
    for j, k in enumerate(('deep', 'shallow')):
        v = fk[:, j]
        res[k + '_null_mean'] = float(v.mean())
        res[k + '_null_sd'] = float(v.std(ddof=1))
        res[k + '_null_mad'] = float(mad_std(v))
        res[k + '_beta'] = float(np.std(v / fe[:, j], ddof=1))     # null sd in units of formal error
        # adopted error: null sd, scaled to this phase's formal error via beta (coverage-aware); take max
        res[k + '_err'] = float(max(res[k + '_null_sd'], res[k + '_beta'] * e[j]))
        res[k + '_p_ge'] = float((np.sum(v >= a[j]) + 1) / (len(v) + 1))
        res[k + '_p_le'] = float((np.sum(v <= a[j]) + 1) / (len(v) + 1))
        res[k + '_null_pct'] = np.percentile(v, [0.5, 2.5, 16, 50, 84, 97.5, 99.5]).tolist()
        res[k + '_null_kurt'] = float(np.mean(((v - v.mean()) / v.std()) ** 4) - 3)
    res['null_corr'] = float(np.corrcoef(fk[:, 0], fk[:, 1])[0, 1])
    res['n_fake'] = int(len(fk))
    inj = {}
    for D in inj_depths:
        rec = []
        for d in fake_offsets(ninj, rng):
            fi = f * (1 - D * tmpl(t, shp, 'd', d, exptime) - s_over_d * D * tmpl(t, shp, 's', d, exptime))
            r_ = measure_at(t, fi, w, proj, shp, exptime, d)
            if r_ is None:
                continue
            r0 = measure_at(t, f, w, proj, shp, exptime, d)
            rec.append(np.r_[r_[0], r0[0]])
        rec = np.array(rec)
        dd = rec[:, 0] - rec[:, 2]; ds = rec[:, 1] - rec[:, 3]       # paired difference: gain
        inj[str(D)] = dict(deep_in=D, shallow_in=s_over_d * D, n=int(len(rec)),
                           deep_rec_mean=float(rec[:, 0].mean()), deep_rec_sd=float(rec[:, 0].std(ddof=1)),
                           shallow_rec_mean=float(rec[:, 1].mean()), shallow_rec_sd=float(rec[:, 1].std(ddof=1)),
                           deep_gain=float(dd.mean() / D), deep_gain_sd=float(dd.std(ddof=1) / D),
                           shallow_gain=float(ds.mean() / (s_over_d * D)))
    res['injection'] = inj
    res['null_samples_deep'] = fk[:, 0].tolist()
    res['null_samples_shallow'] = fk[:, 1].tolist()
    return res, dict(t=t, f=f, w=w, seg=sg, proj=proj)


def phase_fold_binned(t, f, model_base, phc, per=P, half_h=3.0, bin_h=0.25):
    """binned residual flux around phase phc after dividing by baseline"""
    x = dphase(t, phc, per) * 24
    m = np.abs(x) < half_h
    bins = np.arange(-half_h, half_h + 1e-9, bin_h)
    idx = np.digitize(x[m], bins)
    y = (f / model_base - 1)[m]
    out = []
    for i in range(1, len(bins)):
        s = idx == i
        if s.sum() >= 3:
            out.append((0.5 * (bins[i - 1] + bins[i]), y[s].mean(), y[s].std(ddof=1) / np.sqrt(s.sum()), s.sum()))
    return np.array(out)


def e2e_injection(t, f, seg, shp, exptime, order, offsets, depth, s_over_d, minpts=10, floor=2e-4, nclip=4.0,
                  profile=None):
    """End-to-end injection-recovery: inject an eclipse pair at each fake offset into the RAW (normalised)
    series, re-run cleaning (with that phase modelled, as the real eclipses are), weighting and the fit,
    and compare with the uninjected series processed identically (paired).  `profile` = optional
    (tgrid_days, deep_unit_profile, shallow_unit_profile) physical shapes; default = the template.
    Returns arrays (recovered amplitude, paired gain) for deep and shallow."""
    out = []
    for d in offsets:
        if profile is None:
            sd_ = tmpl(t, shp, 'd', d, exptime); ss_ = tmpl(t, shp, 's', d, exptime)
            norm_d = norm_s = 1.0
        else:
            tg, pd_, ps_ = profile
            sd_ = np.interp(dphase(t, shp['d']['phc'] + d, shp['d'].get('per', P), shp['d'].get('t0', T0)), tg, pd_, left=0, right=0)
            ss_ = np.interp(dphase(t, shp['s']['phc'] + d, shp['s'].get('per', P), shp['s'].get('t0', T0)), tg, ps_, left=0, right=0)
            # expected template amplitude of the unit-peak profile
            st = template(tg, 0.0, shp['d']['T'], shp['d']['p'], shp['d']['q'], t0=0.0, per=P)
            norm_d = np.sum(pd_ * st) / np.sum(st * st)
            st = template(tg, 0.0, shp['s']['T'], shp['s']['p'], shp['s']['q'], t0=0.0, per=P)
            norm_s = np.sum(ps_ * st) / np.sum(st * st)
        fi = f * (1 - depth * sd_ - s_over_d * depth * ss_)
        res = []
        for ff in (fi, f):
            keep, sg, w = clean_and_weight(t, ff, seg, shp, exptime, order, nclip=nclip, minpts=minpts, floor=floor, extra_dphi=(d,))
            tt, fk = t[keep], ff[keep]
            pr = Projector(tt, w, sg, order)
            r = measure_at(tt, fk, w, pr, shp, exptime, d)
            res.append(r[0] if r is not None else np.array([np.nan, np.nan]))
        a1, a0 = res
        out.append((a1[0], a1[1], (a1[0] - a0[0]) / (depth * norm_d), (a1[1] - a0[1]) / (s_over_d * depth * norm_s)))
    return np.array(out)
