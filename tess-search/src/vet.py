"""Detailed vetting of a candidate: transit fit, per-sector depths, held-out prediction, odd/even,
secondary eclipse, Gaia neighbours able to mimic the depth, and plots."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import batman
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
import stage2
from bls import diagnostics, robust_std

G = 6.674e-11; RSUN = 6.957e8; MSUN = 1.989e30; REARTH = 6.371e6


def model(t, P, t0, rp, a, b, exptime):
    p = batman.TransitParams()
    inc = np.degrees(np.arccos(min(b / a, 1.0)))
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = t0, P, rp, a, inc, 0.0, 90.0
    p.u, p.limb_dark = [0.45, 0.2], "quadratic"
    ss = max(1, int(exptime / 120))
    return batman.TransitModel(p, t, supersample_factor=ss, exp_time=exptime / 86400).light_curve(p)


def window(secs, P, t0, dur, k=3.0):
    """Data within k durations of each transit (keeps the fit fast)."""
    out = []
    for s in secs:
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5
        m = np.abs(ph * P) < k * dur
        if m.sum() > 10: out.append(dict(s, t=s["t"][m], f=s["f"][m]))
    return out


def fit(secs, P, t0, dur, depth):
    w = window(secs, P, t0, dur)
    exp = [600.0 if s["exptime"] <= 600 else s["exptime"] for s in w]
    a0 = max(2.0, P / (np.pi * dur))

    def resid(x):
        P_, t0_, rp, la, b = x
        r = []
        for s, e in zip(w, exp):
            r.append((s["f"] - model(s["t"], P_, t0_, rp, 10 ** la, b, e)) / s["sigma"])
        return np.concatenate(r)

    x0 = [P, t0, np.sqrt(max(depth, 1e-5)), np.log10(a0), 0.3]
    lo = [P * 0.999, t0 - dur, 1e-3, 0.0, 0.0]; hi = [P * 1.001, t0 + dur, 0.5, 2.5, 1.2]
    best = None
    for b0 in (0.1, 0.5, 0.85):
        x0[4] = b0
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
        if best is None or r.cost < best.cost: best = r
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * max(1.0, 2 * best.cost / max(1, len(best.fun) - 5))
    err = np.sqrt(np.diag(cov))
    P_, t0_, rp, la, b = best.x
    chi_null = np.sum(np.concatenate([(s["f"] - 1) / s["sigma"] for s in w]) ** 2)
    return dict(P=P_, P_err=err[0], t0=t0_, t0_err=err[1], rp=rp, rp_err=err[2], a=10 ** la, a_err=10 ** la * np.log(10) * err[3],
                b=b, b_err=err[4], dchi2=chi_null - 2 * best.cost, npts=len(best.fun))


def rho_star(a, P):
    """Stellar density (solar units) implied by a/R* and P for a circular orbit."""
    rho = 3 * np.pi * a ** 3 / (G * (P * 86400) ** 2)
    return rho / (MSUN / (4 / 3 * np.pi * RSUN ** 3))


def held_out(secs, P, dur, rel=0.002):
    """For each sector, get the ephemeris from the other sectors and measure the depth at the predicted times."""
    out = []
    for i, s in enumerate(secs):
        rest = [x for j, x in enumerate(secs) if j != i]
        if not rest: continue
        j = stage2.joint_fit(rest, P, dur, rel=rel)
        d = diagnostics(s["t"], s["f"], j["P"], j["t0"], j["dur"], s["sigma"])
        # timing uncertainty propagated is not available from BLS; report the depth S/N at the prediction
        out.append(dict(sector=s["sector"], prov=s["prov"], depth=d["depth"], snr=d["depth"] / d["depth_err"] if d["depth_err"] > 0 else np.nan,
                        n_tr=d["n_tr"], P_pred=j["P"]))
    return out


def gaia_neighbours(ra, dec, tmag, depth, radius=63.0):
    """Gaia DR3 stars within `radius` arcsec bright enough to cause `depth` if fully eclipsed."""
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = 500
    q = f"""SELECT source_id, ra, dec, phot_g_mean_mag, phot_rp_mean_mag, phot_bp_mean_mag, ruwe, parallax,
           DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep
           FROM gaiadr3.gaia_source WHERE CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{radius / 3600}))=1
           ORDER BY sep"""
    g = Gaia.launch_job(q).get_results().to_pandas()
    # approximate TESS magnitude from Gaia RP (good to ~0.1-0.2 mag for most stars)
    g["tmag_est"] = g.phot_rp_mean_mag - 0.1
    dm_max = -2.5 * np.log10(depth)  # a neighbour fainter than tmag + dm_max cannot produce the depth
    g["can_mimic"] = g.tmag_est <= tmag + dm_max
    return g


def plot(tic, secs, P, t0, dur, fitres, outpng, title=""):
    ns = len(secs)
    fig, axes = plt.subplots(2, max(ns, 2), figsize=(3.2 * max(ns, 2), 6.0), squeeze=False)
    exp = 600.0
    for i, s in enumerate(secs):
        ax = axes[0, i]
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5
        h = ph * P * 24
        m = np.abs(h) < 4 * dur * 24
        ax.plot(h[m], (s["f"][m] - 1) * 1e6, ".", ms=2, color="0.6")
        bins = np.linspace(-4 * dur * 24, 4 * dur * 24, 25)
        idx = np.digitize(h[m], bins)
        bc = [(bins[k - 1] + bins[k]) / 2 for k in range(1, len(bins)) if np.sum(idx == k) > 2]
        bv = [np.mean((s["f"][m][idx == k] - 1) * 1e6) for k in range(1, len(bins)) if np.sum(idx == k) > 2]
        ax.plot(bc, bv, "o", color="C0", ms=4)
        tt = np.linspace(-4 * dur, 4 * dur, 400)
        ax.plot(tt * 24, (model(tt + fitres["t0"], fitres["P"], fitres["t0"], fitres["rp"], fitres["a"], fitres["b"], exp) - 1) * 1e6, "C3-", lw=1)
        ax.set_title(f"S{s['sector']} {s['prov']}", fontsize=9); ax.set_xlabel("hours from mid-transit", fontsize=8)
        lim = max(4 * fitres["rp"] ** 2 * 1e6, 5 * s["sigma"] * 1e6 / 3)
        ax.set_ylim(-lim, lim / 1.5)
    # combined fold, odd/even, secondary
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    ph = ((t - t0) / P + 0.5) % 1 - 0.5; ep = np.round((t - t0) / P)
    ax = axes[1, 0]
    for k, c, lab in ((0, "C0", "even"), (1, "C1", "odd")):
        m = (ep % 2 == k) & (np.abs(ph * P) < 3 * dur)
        hb = np.linspace(-3 * dur * 24, 3 * dur * 24, 31); idx = np.digitize(ph[m] * P * 24, hb)
        bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if np.sum(idx == q) > 3]
        bv = [np.mean((f[m][idx == q] - 1) * 1e6) for q in range(1, len(hb)) if np.sum(idx == q) > 3]
        ax.plot(bc, bv, "o-", color=c, ms=3, lw=0.8, label=lab)
    ax.legend(fontsize=7); ax.set_title("odd vs even", fontsize=9)
    ax = axes[1, 1]
    ph2 = (ph + 1.0) % 1.0 - 0.5  # secondary at 0
    m = np.abs(ph2 * P) < 3 * dur
    hb = np.linspace(-3 * dur * 24, 3 * dur * 24, 31); idx = np.digitize(ph2[m] * P * 24, hb)
    bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if np.sum(idx == q) > 3]
    bv = [np.mean((f[m][idx == q] - 1) * 1e6) for q in range(1, len(hb)) if np.sum(idx == q) > 3]
    ax.plot(bc, bv, "o-", color="C2", ms=3, lw=0.8); ax.axhline(0, color="k", lw=0.5)
    ax.set_title("phase 0.5 (secondary)", fontsize=9)
    for k in range(2, axes.shape[1]):
        ax = axes[1, k]
        if k == 2:
            phw = ph
            m = np.abs(phw) < 0.5
            hb = np.linspace(-0.5, 0.5, 201); idx = np.digitize(phw[m], hb)
            bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if np.sum(idx == q) > 3]
            bv = [np.mean((f[m][idx == q] - 1) * 1e6) for q in range(1, len(hb)) if np.sum(idx == q) > 3]
            ax.plot(bc, bv, ".", color="k", ms=3); ax.set_title("full phase", fontsize=9)
        else:
            ax.axis("off")
    fig.suptitle(title, fontsize=10); fig.tight_layout(); fig.savefig(outpng, dpi=110); plt.close(fig)
