"""Check 4: time system, in-transit selection, baselines and quality masking.

(a) TESScut and QLP time keywords; TESScut TIME vs QLP TIME at matched FFI cadences (same CADENCENO and nearest time).
(b) Empirical mid-transit offset in the pixel data: pair-aperture (r <= 2 px) background-subtracted light curve,
    per-transit linear baseline from 0.75-2.0 T14, fixed-shape trapezoid (T14 0.1875 d, T12 0.045 d), grid over
    dt0 in [-0.06, +0.06] d with free depth; error from delta chi2 = 1 after scaling chi2/dof to 1.
(c) Depth vs an artificial ephemeris shift (-0.4..+0.4 d): the in-transit selection must peak at zero.
(d) Quality flags: fraction of flagged cadences in-transit vs out-of-transit; bits involved.
(e) centroid.py windows (dur = 0.1249 d, the BLS duration): fraction of its 'OOT' cadences that fall inside the true
    transit (|dt| < T14/2 with the CONTEXT ephemeris).
"""
import numpy as np
from rcommon import Cube, SECTORS, OUT, P, T0, T14, dump

T12 = 0.045
P_V, T0_V, DUR_V = 4.136190053729093, 1658.411589268281, 0.12491872185913565


def trap(dt):
    return np.clip((T14 / 2 - np.abs(dt)) / T12, 0, 1)


def pair_lc(c):
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    ap = np.hypot(xx - c.xt, yy - c.yt) <= 2.0
    f = c.f.reshape(len(c.t), -1)
    bkg = np.nanpercentile(f, 10, axis=1)
    lc = (c.f[:, ap] - bkg[:, None]).sum(1)
    return lc


def detrended(c, lc, P_, T0_, shift=0.0):
    """Per-transit linear baseline (both sides, 0.75-2.0 T14); returns t-rel, relative flux, transit id."""
    g = c.good
    ep = np.round((c.t - T0_ - shift) / P_); dt = c.t - T0_ - shift - ep * P_
    rel, dts, ids = [], [], []
    for e in np.unique(ep[g & (np.abs(dt) < 2.0 * T14)]):
        m = g & (ep == e) & (np.abs(dt) < 2.0 * T14)
        o = m & (np.abs(dt) >= 0.75 * T14)
        if o.sum() < 6 or (m & (np.abs(dt) < 0.4 * T14)).sum() < 2 or (o & (dt < 0)).sum() < 2 or (o & (dt > 0)).sum() < 2:
            continue
        A = np.vstack([np.ones(o.sum()), dt[o]]).T
        coef, *_ = np.linalg.lstsq(A, lc[o], rcond=None)
        base = coef[0] + coef[1] * dt[m]
        rel.append(lc[m] / base - 1); dts.append(dt[m]); ids.append(np.full(m.sum(), e))
    return np.concatenate(dts), np.concatenate(rel), np.concatenate(ids)


def fit_dt0(dts, rel):
    grid = np.linspace(-0.06, 0.06, 241)
    chi = []
    for g0 in grid:
        m = trap(dts - g0)
        d = (m @ -rel) / (m @ m)
        chi.append(((rel + d * m) ** 2).sum())
    chi = np.array(chi)
    k = int(np.argmin(chi))
    m = trap(dts - grid[k]); d = (m @ -rel) / (m @ m)
    s2 = chi[k] / (len(rel) - 2)
    dchi = (chi - chi[k]) / s2
    inside = grid[dchi <= 1]
    return float(grid[k]), float((inside.max() - inside.min()) / 2), float(d)


def main():
    qlp = np.load('/home/user/work/tic87/data/qlp_raw_columns.npy', allow_pickle=True).item()
    res = {}
    for s in SECTORS:
        c = Cube(s)
        r = {}
        r['tesscut_time_keys'] = {k: c.h1.get(k) for k in ('TIMEREF', 'TASSIGN', 'TIMESYS', 'BJDREFI', 'BJDREFF', 'TIMEUNIT')}
        # (a) TESScut vs QLP times
        q = qlp[s]
        tq = np.asarray(q['TIME'], float); cq = np.asarray(q['CADENCENO'])
        with __import__('astropy.io.fits', fromlist=['open']).open(
                __import__('rcommon').cutout_path(s)) as h:
            cad = np.array(h[1].data['CADENCENO'])
        common, i1, i2 = np.intersect1d(cad, cq, return_indices=True)
        if len(common) > 10:
            d = (c.t[i1] - tq[i2]) * 86400.0
            d = d[np.isfinite(d)]
            r['tesscut_minus_qlp_same_cadenceno_s'] = dict(n=int(len(d)), median=float(np.median(d)),
                                                         p1=float(np.percentile(d, 1)), p99=float(np.percentile(d, 99)))
        idx = np.searchsorted(tq, c.t).clip(1, len(tq) - 1)
        near = np.where(np.abs(tq[idx] - c.t) < np.abs(tq[idx - 1] - c.t), idx, idx - 1)
        d2 = (c.t - tq[near]) * 86400.0
        d2 = d2[np.isfinite(d2)]
        r['tesscut_minus_qlp_nearest_s'] = dict(median=float(np.median(d2)), mad=float(np.median(np.abs(d2 - np.median(d2)))))
        # (b) transit timing in the pixels
        lc = pair_lc(c)
        dts, rel, ids = detrended(c, lc, P, T0)
        dt0, edt0, depth = fit_dt0(dts, rel)
        r['pair_aperture_dt0_d'] = dt0; r['pair_aperture_dt0_err_d'] = edt0; r['pair_aperture_depth_ppt'] = 1e3 * depth
        r['n_transits_used'] = int(len(np.unique(ids)))
        # predicted ephemeris uncertainty at this sector
        n_ep = np.round((np.median(c.t) - T0) / P)
        r['ephemeris_sigma_d'] = float(np.hypot(0.0027, n_ep * 7e-6))
        # (c) depth vs artificial shift
        shifts = np.round(np.arange(-0.4, 0.401, 0.05), 3)
        dep = []
        for sh in shifts:
            dd, rr, _ = detrended(c, lc, P, T0, shift=sh)
            m = trap(dd)
            dep.append(float(1e3 * (m @ -rr) / (m @ m)))
        r['depth_vs_shift'] = dict(shift_d=shifts.tolist(), depth_ppt=dep)
        # (d) quality flags
        ep = np.round((c.t - T0) / P); dt = c.t - T0 - ep * P
        fin = np.isfinite(c.t)
        intr = fin & (np.abs(dt) < T14 / 2); near_oot = fin & (np.abs(dt) >= 0.75 * T14) & (np.abs(dt) < 2.5 * T14)
        flagged = c.q != 0
        r['flag_frac_in_transit'] = float(flagged[intr].mean()); r['flag_frac_near_oot'] = float(flagged[near_oot].mean())
        r['flag_frac_all'] = float(flagged[fin].mean())
        bits = {}
        for b in range(0, 17):
            nb = int(((c.q[intr] >> b) & 1).sum())
            if nb:
                bits[str(2 ** b)] = nb
        r['flag_bits_in_transit_counts'] = bits
        r['n_in_transit_cadences'] = int(intr.sum())
        # (e) centroid.py windows
        epv = np.round((c.t - T0_V) / P_V); dtv = c.t - T0_V - epv * P_V
        ootv = (np.abs(dtv) > 0.75 * DUR_V) & (np.abs(dtv) < 1.75 * DUR_V) & (c.q == 0)
        r['centroidpy_oot_frac_inside_true_transit'] = float((np.abs(dt[ootv]) < T14 / 2).mean())
        r['centroidpy_oot_frac_inside_transit_T14_0.21'] = float((np.abs(dt[ootv]) < 0.21 / 2).mean())
        r['centroidpy_vs_context_ephemeris_offset_min'] = float(
            ((T0_V + np.round((np.median(c.t) - T0_V) / P_V) * P_V) - (T0 + n_ep * P)) * 1440)
        res[str(s)] = r
        print(f"S{s}: keys {r['tesscut_time_keys']}; TESScut-QLP (same CADENCENO) "
              f"{r.get('tesscut_minus_qlp_same_cadenceno_s')} nearest {r['tesscut_minus_qlp_nearest_s']}; "
              f"pixel dt0 = {dt0*1440:.1f} +- {edt0*1440:.1f} min (ephemeris sigma {r['ephemeris_sigma_d']*1440:.1f} min), "
              f"depth {1e3*depth:.2f} ppt, n_tr {r['n_transits_used']}; peak shift {shifts[int(np.argmax(dep))]:+.2f} d; "
              f"flags in/oot/all {r['flag_frac_in_transit']:.3f}/{r['flag_frac_near_oot']:.3f}/{r['flag_frac_all']:.3f} bits {bits}; "
              f"centroid.py OOT inside transit {r['centroidpy_oot_frac_inside_true_transit']:.3f}, "
              f"ephem offset {r['centroidpy_vs_context_ephemeris_offset_min']:.1f} min", flush=True)
    dump(res, f'{OUT}/r4_timing.json')


if __name__ == '__main__':
    main()
