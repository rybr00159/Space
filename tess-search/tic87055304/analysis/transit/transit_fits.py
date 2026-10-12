"""Per-transit fits of TIC 87055304 with a fixed best-fit template: mid-times (TTV analysis) and depths.

For every predicted transit:
  window |t - Tc_pred| < W (0.25 d primary; 0.40 d check, env TF_WINDOW), restricted to the momentum-dump/gap segment containing mid-transit;
  coverage = (in-transit cadences x exptime) / T14; pre/post = out-of-transit time span before ingress/after egress.
  Selection for timing: coverage >= 0.70 and >= 0.5 h of out-of-transit data on BOTH sides (needed so the local
  linear baseline does not trade off against the mid-time).
  Timing fit: model = (c0 + c1 (t - Tc_pred)) * template(t - tc); c0, c1 solved linearly for each trial tc; tc on a
  0.5-min grid over +-0.1 d, refined with a bounded scalar minimiser; formal 1-sigma from the Delta chi^2 = 1 points
  with the white noise set to the rms of the best-fit residuals.
  Depth fit: model = (c0 + c1 (t - Tc)) * (1 - A (1 - template(t - Tc))), Tc fixed at the linear ephemeris, A free.
Injection calibration (per sector): the template is injected (multiplicatively into the same flux) at random times
whose whole window avoids real transits; the same selection and fits are applied; the robust scatter (1.4826 MAD) and
the standard deviation of (recovered - injected) times and depths, relative to the formal errors, give per-sector
error scale factors (applied as max(1, scale)).
Usage: python transit_fits.py [template_json] [flux_key f_raw|f] [tag]
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize_scalar
import tmodel as tm

HERE = os.path.dirname(os.path.abspath(__file__))
W = float(os.environ.get("TF_WINDOW", "0.25"))  # half-window [d]; 0.25 primary, 0.40 as a check
COV_MIN, OOT_MIN = 0.70, 0.5 / 24
GRID_HALF, GRID_STEP = 0.10, 0.5 / 1440


def template_from(json_path):
    js = json.load(open(json_path))
    if "best" in js:   # mcmc summary: use posterior medians of the shape parameters
        p = js["params"]
        th = [p[n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]]
    else:
        th = js["x"][:7]
    return np.array(th, float)


class Template:
    def __init__(self, th7):
        self.th = th7
        g = tm.gm()
        T0, P, k, b, rho, u1, u2 = th7
        self.prof = g.profile(P, k, tm.a_from_rho(rho, P), b, u1, u2)
        self.g = g.g; self.dt = g.dt
        C = np.concatenate([[0.0], np.cumsum(0.5 * (self.prof[1:] + self.prof[:-1]) * self.dt)])
        self.F = {}
        for E in (200.0, 600.0, 1800.0):
            m = int(round(E / 86400 / 2 / self.dt)); Fe = np.ones_like(self.prof)
            Fe[m:-m] = (C[2 * m:] - C[:-2 * m]) / (2 * m * self.dt); self.F[E] = Fe
        der = tm.derived(np.array([th7]))
        self.T14 = float(der["T14_h"][0]) / 24
        self.depth = 1 - self.prof.min()

    def __call__(self, tau, E):
        return np.interp(tau, self.g, self.F[E], left=1.0, right=1.0)


def lin_solve(t, f, mdl, tref):
    X = np.vstack([mdl, mdl * (t - tref)]).T
    c, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ c
    return c, r


def window_data(t, f, seg, Tc, T14, E):
    m = np.abs(t - Tc) < W
    if m.sum() < 5:
        return None
    tt, ff, ss = t[m], f[m], seg[m]
    inn = np.abs(tt - Tc) < T14 / 2
    if inn.sum() == 0:
        return None
    # segment containing mid-transit (majority of in-transit cadences)
    vals, cnt = np.unique(ss[inn], return_counts=True)
    k = vals[np.argmax(cnt)]
    keep = ss == k
    tt, ff = tt[keep], ff[keep]
    inn = np.abs(tt - Tc) < T14 / 2
    cov = min(1.0, inn.sum() * E / 86400 / T14)
    pre = tt[tt < Tc - T14 / 2]; post = tt[tt > Tc + T14 / 2]
    pre_span = (Tc - T14 / 2 - pre.min()) if len(pre) else 0.0
    post_span = (post.max() - (Tc + T14 / 2)) if len(post) else 0.0
    # require actual cadences on each side: span counts only if data exist
    pre_n, post_n = len(pre), len(post)
    return dict(t=tt, f=ff, cov=cov, pre=pre_span if pre_n > 2 else 0.0, post=post_span if post_n > 2 else 0.0)


def fit_time(w, Tc, tmpl, E):
    t, f = w["t"], w["f"]
    def chi(tc):
        mdl = tmpl(t - tc, E)
        c, r = lin_solve(t, f, mdl, Tc)
        return np.sum(r ** 2)
    grid = Tc + np.arange(-GRID_HALF, GRID_HALF + GRID_STEP / 2, GRID_STEP)
    ch = np.array([chi(x) for x in grid])
    i = int(np.argmin(ch))
    lo, hi = grid[max(i - 2, 0)], grid[min(i + 2, len(grid) - 1)]
    res = minimize_scalar(chi, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-6))
    tc = res.x; c2 = res.fun
    sig2 = c2 / max(1, len(t) - 3)
    dchi = (ch - c2) / sig2
    # Delta chi2 = 1 crossings (using the grid + interpolation)
    def cross(side):
        idx = np.arange(i, -1, -1) if side < 0 else np.arange(i, len(grid))
        for j0, j1 in zip(idx[:-1], idx[1:]):
            if dchi[j1] >= 1 > dchi[j0]:
                return grid[j0] + (1 - dchi[j0]) / (dchi[j1] - dchi[j0]) * (grid[j1] - grid[j0])
        return np.nan
    l, h = cross(-1), cross(+1)
    err = np.nanmean([tc - l, h - tc]) if np.isfinite(l) or np.isfinite(h) else np.nan
    at_edge = (i == 0) or (i == len(grid) - 1)
    # depth S/N of this event at best tc: compare with flat model
    c0, r0 = lin_solve(t, f, np.ones_like(t), Tc)
    dchi_null = (np.sum(r0 ** 2) - c2) / sig2
    return dict(tc=float(tc), err=float(err), err_lo=float(tc - l), err_hi=float(h - tc), rms=float(np.sqrt(sig2)),
                at_edge=bool(at_edge), dchi_null=float(dchi_null), n=int(len(t)))


def fit_depth(w, Tc, tmpl, E):
    t, f = w["t"], w["f"]
    s = 1 - tmpl(t - Tc, E)  # positive in transit
    # model = (c0 + c1 dt) * (1 - A s) = c0 + c1 dt - c0 A s - c1 A s dt ; linearise (c1 A s dt negligible)
    X = np.vstack([np.ones_like(t), t - Tc, -s]).T
    c, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ c
    sig2 = np.sum(r ** 2) / max(1, len(t) - 3)
    cov = np.linalg.inv(X.T @ X) * sig2
    A = c[2] / c[0]
    return dict(A=float(A), A_err=float(np.sqrt(cov[2, 2]) / c[0]), rms=float(np.sqrt(sig2)))


def run(template_json, key="f_raw", tag="raw", n_inj=400, seed=7):
    th = template_from(template_json)
    tmpl = Template(th)
    T0, P = th[0], th[1]
    lc = np.load(os.path.join(HERE, "lc_prepared.npy"), allow_pickle=True).item()
    rng = np.random.default_rng(seed)
    rows, inj = [], {}
    for s in tm.SECTORS:
        v = lc[s]; E = tm.EXPT[s]
        t, f, seg = v["t"], v[key], v["seg"]
        n_all = np.arange(np.ceil((t.min() - W - T0) / P), np.floor((t.max() + W - T0) / P) + 1).astype(int)
        for n in n_all:
            Tc = T0 + n * P
            w = window_data(t, f, seg, Tc, tmpl.T14, E)
            if w is None:
                rows.append(dict(sector=s, epoch=int(n), Tc_pred=float(Tc), cov=0.0, used=False, reason="no data")); continue
            row = dict(sector=s, epoch=int(n), Tc_pred=float(Tc), cov=float(w["cov"]), pre_h=float(w["pre"] * 24), post_h=float(w["post"] * 24))
            ok = w["cov"] >= COV_MIN and w["pre"] >= OOT_MIN and w["post"] >= OOT_MIN
            row["used"] = bool(ok)
            row["reason"] = "" if ok else ("coverage" if w["cov"] < COV_MIN else "baseline")
            if w["cov"] >= 0.3 and (w["pre"] > 0 or w["post"] > 0):
                row.update({f"time_{k}": val for k, val in fit_time(w, Tc, tmpl, E).items()})
                row.update({f"depth_{k}": val for k, val in fit_depth(w, Tc, tmpl, E).items()})
            rows.append(row)
        # injections
        real = T0 + n_all * P
        res = []
        tries = 0
        while len(res) < n_inj and tries < 20 * n_inj:
            tries += 1
            Ti = rng.uniform(t.min() + 0.1, t.max() - 0.1)
            if np.min(np.abs(real - Ti)) < W + tmpl.T14 / 2 + 0.05:
                continue
            fi = f * tmpl(t - Ti, E)
            w = window_data(t, fi, seg, Ti, tmpl.T14, E)
            if w is None or not (w["cov"] >= COV_MIN and w["pre"] >= OOT_MIN and w["post"] >= OOT_MIN):
                continue
            ft = fit_time(w, Ti, tmpl, E); fd = fit_depth(w, Ti, tmpl, E)
            res.append((Ti, ft["tc"] - Ti, ft["err"], fd["A"] - 1, fd["A_err"], ft["at_edge"]))
        res = np.array(res)
        dt, et, dA, eA, edge = res[:, 1], res[:, 2], res[:, 3], res[:, 4], res[:, 5].astype(bool)
        mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
        inj[s] = dict(n=int(len(res)), n_edge=int(edge.sum()), time_bias_min=float(np.median(dt) * 1440),
                      time_mad_min=float(mad(dt) * 1440), time_sd_min=float(np.std(dt) * 1440),
                      time_formal_med_min=float(np.nanmedian(et) * 1440),
                      time_scale_mad=float(np.nanmedian(np.abs(dt) / et) / 0.6745),
                      time_scale_sd=float(np.std(dt / et)),
                      depth_bias=float(np.median(dA)), depth_mad=float(mad(dA)), depth_sd=float(np.std(dA)),
                      depth_formal_med=float(np.median(eA)), depth_scale_mad=float(np.median(np.abs(dA) / eA) / 0.6745),
                      depth_scale_sd=float(np.std(dA / eA)))
        np.save(os.path.join(HERE, f"inj_{tag}_S{s}.npy"), res)
        print(s, inj[s], flush=True)
    out = dict(template=th.tolist(), template_T14_h=tmpl.T14 * 24, template_depth_ppm=tmpl.depth * 1e6, key=key,
               W=W, cov_min=COV_MIN, oot_min_h=OOT_MIN * 24, transits=rows, injection=inj)
    json.dump(out, open(os.path.join(HERE, f"transit_fits_{tag}.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    tj = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "ml_start.json")
    key = sys.argv[2] if len(sys.argv) > 2 else "f_raw"
    tag = sys.argv[3] if len(sys.argv) > 3 else "raw"
    nin = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    run(tj, key, tag, n_inj=nin)
