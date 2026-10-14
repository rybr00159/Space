"""v2 library for the chromatic eclipse-depth study of TIC 71316387 (re-derivation).

Re-uses the validated pieces of ../../scripts/chromlib.py (block least squares, V template;
lsq() checked against the dense solver to 1e-13) and adds:
  * eclipse-window extraction with real-eclipse masking for nulls,
  * generic measurement with fake-phase null in which the REAL eclipses are always modelled
    (so fake windows overlapping real eclipses are not biased),
  * error-scaling factor beta = std(a_fake / sigma_formal_fake), so phases with different
    coverage are normalised properly.
Conventions: times BJD_TDB; deep-eclipse ephemeris T0 = 2460165.6697, P = 1.3036712 d.
Model: f = per-segment polynomial baseline - A_d*s_d(t) - A_s*s_s(t), s = unit-depth V template
(1-|x|^p)^q integrated over the exposure.  A is the central (maximum) fractional depth.
"""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/scripts')
import chromlib as C   # noqa: E402

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad'
CH = BASE + '/tic713/chromatic'
V2 = CH + '/v2'
P, T0, P_ERR = C.P, C.T0, C.P_ERR
lsq, mad_std, segments = C.lsq, C.mad_std, C.segments


def load_shapes():
    return json.load(open(V2 + '/results/tess_ref2.json'))['shapes']


def tmpl(t, shp, which, dphi=0.0, exptime=0.0):
    s = shp[which]
    nsub = max(1, int(np.ceil(exptime / (120 / 86400.)))) if exptime > 0 else 1
    return C.template(t, s['phc'] + dphi, s['T'], s['p'], s['q'], exptime, nsub)


def dist_to_phase(t, phc):
    """|time| (d) from the nearest occurrence of phase phc."""
    ph = ((t - T0) / P - phc) % 1.0
    ph = np.where(ph > 0.5, ph - 1, ph)
    return np.abs(ph) * P


def fake_phases(n, rng, avoid=(0.0, 0.5), halfwidth=0.06):
    """phase offsets d such that the fake deep (0+d) and fake shallow (0.5+d) windows are
    > halfwidth (phase) from the real eclipses. halfwidth 0.06 = 1.88 h vs eclipse half-width 0.77 h."""
    out = []
    while len(out) < n:
        d = rng.uniform(0, 1)
        ok = True
        for a in avoid:
            for c in (d, (d + 0.5) % 1):
                if abs(((c - a) + 0.5) % 1 - 0.5) < halfwidth:
                    ok = False
        if ok:
            out.append(d)
    return np.array(out)


def measure(t, f, w, seg, shp, exptime, order=1, dphi=0.0):
    """Fit amplitudes at phase offset dphi; for dphi != 0 the real eclipses are fitted too."""
    tem = [tmpl(t, shp, 'd', dphi, exptime), tmpl(t, shp, 's', dphi, exptime)]
    if dphi != 0.0:
        tem += [tmpl(t, shp, 'd', 0.0, exptime), tmpl(t, shp, 's', 0.0, exptime)]
    # guard: a template column identically zero -> singular; drop and return nan
    good = [np.any(s > 0) for s in tem]
    if not all(good[:2]):
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), None
    tem_g = [s for s, g in zip(tem, good) if g]
    a, ae, r, m = lsq(t, f, w, seg, tem_g, order)
    return a[:2], ae[:2], r


def clean(t, f, seg, shp, exptime, order=1, nclip=4.0, niter=6, minpts=8, wfloor=1e-4):
    """Iterative clipping on residuals of full model; per-segment robust sigma -> weights."""
    keep = np.ones(len(t), bool)
    for _ in range(niter):
        tt, ff = t[keep], f[keep]
        _, sg = np.unique(seg[keep], return_inverse=True)
        a, ae, r = measure(tt, ff, np.ones_like(tt), sg, shp, exptime, order)
        sig = np.empty(len(tt))
        for k in np.unique(sg):
            m = sg == k
            sig[m] = mad_std(r[m]) if m.sum() >= 5 else np.inf
        bad = np.abs(r) > nclip * sig
        idx = np.where(keep)[0]
        nk = keep.copy(); nk[idx[bad]] = False
        cnt = np.bincount(seg[nk], minlength=seg.max() + 1)
        nk &= cnt[seg] >= minpts
        if nk.sum() == keep.sum():
            break
        keep = nk
    tt, ff = t[keep], f[keep]
    _, sg = np.unique(seg[keep], return_inverse=True)
    a, ae, r = measure(tt, ff, np.ones_like(tt), sg, shp, exptime, order)
    w = np.empty(len(tt))
    for k in np.unique(sg):
        m = sg == k
        w[m] = 1.0 / max(mad_std(r[m]), wfloor) ** 2
    return keep, sg, w


def run_survey(name, t, f, seg, shp, exptime, order=1, nfake=400, ninj=200, inj_depths=(0.003, 0.006),
               seed=1, nclip=4.0, minpts=8, ratio_s=0.47, wfloor=1e-4):
    """Clean, measure, fake-phase null (with beta scaling), injection-recovery."""
    rng = np.random.default_rng(seed)
    keep, sg, w = clean(t, f, seg, shp, exptime, order, nclip=nclip, minpts=minpts, wfloor=wfloor)
    t, f = t[keep], f[keep]
    a, ae, r = measure(t, f, w, sg, shp, exptime, order)
    chi2r = np.sum(w * r ** 2) / (len(t) - (order + 1) * (sg.max() + 1) - 2)
    sd = tmpl(t, shp, 'd', 0, exptime); ss = tmpl(t, shp, 's', 0, exptime)
    res = dict(name=name, n=int(len(t)), nseg=int(sg.max() + 1), order=order, nclip=nclip,
               n_in_deep=int((sd > 0.1).sum()), n_in_shallow=int((ss > 0.1).sum()),
               n_deep_events=int(len(np.unique(np.round((t[sd > 0.1] - T0) / P)))),
               n_shallow_events=int(len(np.unique(np.round((t[ss > 0.1] - T0) / P - 0.5)))),
               rms_median=float(np.median(1 / np.sqrt(w))), chi2r=float(chi2r),
               deep=float(a[0]), shallow=float(a[1]), deep_err_formal=float(ae[0]), shallow_err_formal=float(ae[1]))
    fk, fke = [], []
    for d in fake_phases(nfake, rng):
        aa, ee, _ = measure(t, f, w, sg, shp, exptime, order, d)
        if np.all(np.isfinite(aa)):
            fk.append(aa); fke.append(ee)
    fk, fke = np.array(fk), np.array(fke)
    z = fk / fke
    for j, k in enumerate(('deep', 'shallow')):
        res[f'{k}_null_mean'] = float(np.mean(fk[:, j]))
        res[f'{k}_null_sd'] = float(np.std(fk[:, j], ddof=1))
        res[f'{k}_null_mad'] = float(mad_std(fk[:, j]))
        res[f'{k}_beta'] = float(np.std(z[:, j], ddof=1))
        res[f'{k}_zmean'] = float(np.mean(z[:, j]))
        # adopted error: larger of null sd and beta-scaled formal error
        res[f'{k}_err'] = float(max(res[f'{k}_null_sd'], res[f'{k}_beta'] * ae[j]))
        res[f'{k}_p_null_upper'] = float((np.sum(fk[:, j] >= a[j]) + 1) / (len(fk) + 1))
        res[f'{k}_p_null_lower'] = float((np.sum(fk[:, j] <= a[j]) + 1) / (len(fk) + 1))
        q = np.percentile(fk[:, j], [0.5, 2.5, 16, 50, 84, 97.5, 99.5])
        res[f'{k}_null_pct'] = q.tolist()
    res['null_corr'] = float(np.corrcoef(fk[:, 0], fk[:, 1])[0, 1])
    res['n_fake'] = int(len(fk))
    res['null_samples'] = fk.tolist()
    inj = {}
    for D in inj_depths:
        rec = []
        for d in fake_phases(ninj, rng):
            fi = f * (1 - D * tmpl(t, shp, 'd', d, exptime) - ratio_s * D * tmpl(t, shp, 's', d, exptime))
            aa, _, _ = measure(t, fi, w, sg, shp, exptime, order, d)
            if np.all(np.isfinite(aa)):
                rec.append(aa)
        rec = np.array(rec)
        inj[str(D)] = dict(deep_in=D, shallow_in=ratio_s * D,
                           deep_rec_mean=float(np.mean(rec[:, 0])), deep_rec_sd=float(np.std(rec[:, 0], ddof=1)),
                           shallow_rec_mean=float(np.mean(rec[:, 1])), shallow_rec_sd=float(np.std(rec[:, 1], ddof=1)),
                           deep_gain=float((np.mean(rec[:, 0]) - res['deep_null_mean']) / D),
                           shallow_gain=float((np.mean(rec[:, 1]) - res['shallow_null_mean']) / (ratio_s * D)),
                           deep_gain_err=float(np.std(rec[:, 0], ddof=1) / np.sqrt(len(rec)) / D))
    res['injection'] = inj
    return res, (t, f, w, sg)
