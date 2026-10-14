"""Generic eclipse-depth measurement for a ground-based / sparse light curve with
fake-phase null and injection-recovery.

Usage (as a module):  run(name, t_bjd, f, seg, exptime, ...)
  t_bjd  : BJD_TDB
  f      : relative flux (normalised to ~1 per segment is NOT required: baselines are free)
  seg    : integer segment id (night x camera); each gets its own polynomial baseline
Outputs a dict with measured depths, empirical-null errors, injection-recovery results.
"""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/scripts')
import chromlib as C

OUT = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic'
REF = json.load(open(OUT + '/results/tess_ref.json'))
SH = REF['shapes']
PHC = REF['phase_centres']


def tmpl(t, which, dphi=0.0, exptime=0.0, shp=None, t0=C.T0, per=C.P):
    s = (shp or SH)[which]
    nsub = max(1, int(round(exptime / (60 / 86400))))
    return C.template(t, PHC[which] + dphi, s['T_d'], s['p'], s['q'], exptime, nsub, t0=t0, per=per)


def clean_and_weight(t, f, seg, exptime, order=1, nclip=4.0, niter=5, minpts=8):
    """Iterative sigma clipping on residuals of the full model, per-segment robust weights."""
    keep = np.ones(len(t), bool)
    for it in range(niter):
        tt, ff, sg = t[keep], f[keep], seg[keep]
        _, sg = np.unique(sg, return_inverse=True)
        tem = [tmpl(tt, 'd', 0, exptime), tmpl(tt, 's', 0, exptime)]
        a, ae, r, m = C.lsq(tt, ff, np.ones_like(tt), sg, tem, order)
        sig = np.zeros(len(tt))
        for k in np.unique(sg):
            mm = sg == k
            sig[mm] = C.mad_std(r[mm]) if mm.sum() >= 5 else np.inf
        bad = np.abs(r) > nclip * sig
        idx = np.where(keep)[0]
        newkeep = keep.copy()
        newkeep[idx[bad]] = False
        if newkeep.sum() == keep.sum():
            break
        keep = newkeep
    # drop segments with too few points
    sgk = seg[keep]
    cnt = np.bincount(sgk, minlength=seg.max() + 1)
    keep &= (cnt[seg] >= minpts)
    # final per-segment sigma
    tt, ff, sg = t[keep], f[keep], seg[keep]
    _, sg = np.unique(sg, return_inverse=True)
    tem = [tmpl(tt, 'd', 0, exptime), tmpl(tt, 's', 0, exptime)]
    a, ae, r, m = C.lsq(tt, ff, np.ones_like(tt), sg, tem, order)
    w = np.zeros(len(tt))
    for k in np.unique(sg):
        mm = sg == k
        s = C.mad_std(r[mm])
        w[mm] = 1.0 / max(s, 1e-4) ** 2
    return keep, sg, w


def measure(t, f, w, seg, exptime, order=1, dphi=0.0, real_too=True, shp=None):
    tem = [tmpl(t, 'd', dphi, exptime, shp), tmpl(t, 's', dphi, exptime, shp)]
    if real_too and dphi != 0.0:
        tem += [tmpl(t, 'd', 0.0, exptime, shp), tmpl(t, 's', 0.0, exptime, shp)]
    a, ae, r, m = C.lsq(t, f, w, seg, tem, order)
    return a, ae, r


def run(name, t, f, seg, exptime, order=1, nfake=300, inj_depths=(0.003, 0.006), ninj=150, seed=42,
        extra_info=None, minpts=8, nclip=4.0):
    rng = np.random.default_rng(seed)
    assert abs(np.median(f) - 1) < 0.1, 'flux must be normalised to ~1 (fractional depths)'
    keep, sg, w = clean_and_weight(t, f, seg, exptime, order, nclip=nclip, minpts=minpts)
    t, f = t[keep], f[keep]
    a, ae, r = measure(t, f, w, sg, exptime, order)
    # in-eclipse sample counts
    sd = tmpl(t, 'd', 0, exptime); ss = tmpl(t, 's', 0, exptime)
    res = dict(name=name, n=int(len(t)), nseg=int(sg.max() + 1), order=order,
               n_in_deep=int((sd > 0.05).sum()), n_in_shallow=int((ss > 0.05).sum()),
               n_deep_events=int(len(np.unique(np.round((t[sd > 0.05] - C.T0) / C.P)))),
               n_shallow_events=int(len(np.unique(np.round((t[ss > 0.05] - C.T0) / C.P - 0.5)))),
               sum_w_s2_deep=float(np.sum(w * sd ** 2)), rms_median=float(np.median(1 / np.sqrt(w))),
               deep=float(a[0]), shallow=float(a[1]), deep_err_formal=float(ae[0]), shallow_err_formal=float(ae[1]))
    # fake-phase null
    fk = []
    for dphi in C.fake_phases(nfake, rng):
        aa, _, _ = measure(t, f, w, sg, exptime, order, dphi)
        fk.append(aa[:2])
    fk = np.array(fk)
    res['null_deep_mean'] = float(np.mean(fk[:, 0])); res['null_shallow_mean'] = float(np.mean(fk[:, 1]))
    res['deep_err_null'] = float(np.std(fk[:, 0], ddof=1)); res['shallow_err_null'] = float(np.std(fk[:, 1], ddof=1))
    res['deep_err_null_mad'] = float(C.mad_std(fk[:, 0])); res['shallow_err_null_mad'] = float(C.mad_std(fk[:, 1]))
    res['deep_null_p'] = float((np.sum(fk[:, 0] >= a[0]) + 1) / (len(fk) + 1))
    res['shallow_null_p'] = float((np.sum(fk[:, 1] >= a[1]) + 1) / (len(fk) + 1))
    res['null_samples_deep'] = fk[:, 0].tolist(); res['null_samples_shallow'] = fk[:, 1].tolist()
    # injection-recovery at fake phases: multiplicative injection of TESS-shaped eclipses
    inj = {}
    for D in inj_depths:
        rec = []
        for dphi in C.fake_phases(ninj, rng):
            fi = f * (1 - D * tmpl(t, 'd', dphi, exptime) - 0.47 * D * tmpl(t, 's', dphi, exptime))
            aa, _, _ = measure(t, fi, w, sg, exptime, order, dphi)
            rec.append(aa[:2])
        rec = np.array(rec)
        inj[str(D)] = dict(deep_in=D, shallow_in=0.47 * D, deep_rec_mean=float(np.mean(rec[:, 0])),
                           deep_rec_sd=float(np.std(rec[:, 0], ddof=1)), shallow_rec_mean=float(np.mean(rec[:, 1])),
                           shallow_rec_sd=float(np.std(rec[:, 1], ddof=1)),
                           deep_ratio=float(np.mean(rec[:, 0]) / D), shallow_ratio=float(np.mean(rec[:, 1]) / (0.47 * D)))
    res['injection'] = inj
    # null-corrected depth: subtract null mean (baseline bias), error = null sd
    res['deep_corr'] = res['deep'] - res['null_deep_mean']
    res['shallow_corr'] = res['shallow'] - res['null_shallow_mean']
    if extra_info:
        res.update(extra_info)
    return res, (t, f, w, sg)


def phase_scan(t, f, w, sg, exptime, order=1, grid=None):
    """Deep-eclipse amplitude vs trial phase offset (ephemeris check)."""
    if grid is None:
        grid = np.linspace(-0.06, 0.06, 121)
    out = []
    for g in grid:
        tem = [tmpl(t, 'd', g, exptime), tmpl(t, 's', g, exptime)]
        a, ae, r = C.lsq(t, f, w, sg, tem, order)[:3]
        out.append((g, a[0], a[1], np.sum(w * r ** 2)))
    return np.array(out)
