"""Characterise eclipsing-binary candidates: true period, primary/secondary eclipses, novelty, localisation."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
import stage2, known
from bls import _bls
from centroid import localise

W = "/home/user/work"


def fold_profile(t, f, P, t0, nb=400):
    ph = ((t - t0) / P) % 1.0
    idx = np.minimum((ph * nb).astype(int), nb - 1)
    c = np.bincount(idx, minlength=nb); s = np.bincount(idx, f, minlength=nb)
    prof = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    return (np.arange(nb) + 0.5) / nb, prof, c


def refine_period(secs, P, dur, rel=0.002):
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    w = np.concatenate([np.full(len(s["t"]), 1 / s["sigma"] ** 2) for s in secs])
    base = t.max() - t.min(); dP = P * (dur / 8) / base
    n = int(min(30000, 2 * rel * P / dP) + 1)
    fr = 1 / np.linspace(P * (1 - rel), P * (1 + rel), n)
    y = f - 1; y = y - np.sum(w * y) / np.sum(w)
    pw, dp, t0, du = _bls(t - t.min(), y, w, fr, np.array([dur / 1.5, dur, dur * 1.5]), min(dur / 8, 5 / 1440))
    i = int(np.argmax(pw))
    return 1 / fr[i], t0[i] + t.min(), du[i]


def eclipse_times(secs, thresh_frac=0.5):
    """Times and depths of individual eclipses (runs of low flux), refined with a parabola."""
    ev = []
    for s in secs:
        t, f = s["t"], s["f"]
        dmax = 1 - np.percentile(f, 0.2)
        if dmax < 8 * s["sigma"]: continue
        idx = np.where(f < 1 - thresh_frac * dmax * 0.3)[0]
        if not len(idx): continue
        for g in np.split(idx, np.where(np.diff(t[idx]) > 0.15)[0] + 1):
            if len(g) < 3: continue
            k = g[np.argmin(f[g])]
            sel = g[np.argsort(f[g])[:min(7, len(g))]]
            c = np.polyfit(t[sel] - t[k], f[sel], 2)
            tm = t[k] - c[1] / (2 * c[0]) if c[0] > 0 and abs(c[1] / (2 * c[0])) < 0.1 else t[k]
            ev.append(dict(sector=s["sector"], t=float(tm), depth=float(1 - f[g].min())))
    return ev


def timing_ephemeris(ev, P_guess):
    """Linear ephemeris of the deeper eclipse family by cycle counting, refined sector by sector."""
    if len(ev) < 3: return None
    d = np.array([e["depth"] for e in ev]); t = np.array([e["t"] for e in ev])
    deep = d > 0.5 * (d.max() + np.median(d)) if (d.max() - d.min()) > 0.15 * d.max() else np.ones(len(d), bool)
    tp = np.sort(t[deep])
    # local period from the closest pairs, then bootstrap the cycle count outward in time
    gaps = np.diff(tp); P = np.median(gaps[gaps < 1.5 * np.median(gaps)]) if len(gaps) else P_guess
    if not np.isfinite(P) or P <= 0: P = P_guess
    T0 = tp[0]
    for _ in range(5):
        n = np.round((tp - T0) / P)
        A = np.vstack([np.ones_like(n), n]).T
        (T0, P), *_ = np.linalg.lstsq(A, tp, rcond=None)
    n = np.round((tp - T0) / P); resid = tp - (T0 + n * P)
    return dict(P=float(P), T0=float(T0), rms_min=float(np.std(resid) * 1440), n_ecl=int(len(tp)),
                ncycles=int(n.max() - n.min()))


def analyse(tic, P, t0, dur, star, outdir):
    os.makedirs(outdir, exist_ok=True)
    secs = stage2.load_all(tic)
    res = dict(tic=tic, star=star, sectors=[s["sector"] for s in secs])
    if not secs: return res
    ev = eclipse_times(secs)
    try:
        eph = timing_ephemeris(ev, P) if len(ev) >= 4 else None
    except Exception:
        eph = None
    res["eclipses"] = ev; res["timing_ephemeris"] = eph
    if eph and eph["rms_min"] < 30:
        # timing ephemeris of the deeper eclipse; let the BLS refine only very locally around it
        P1, t01, d1 = refine_period(secs, eph["P"], dur, rel=0.0005)
        if abs(P1 - eph["P"]) / eph["P"] > 2e-4: P1, t01 = eph["P"], eph["T0"]
    else:
        P1, t01, d1 = refine_period(secs, P, dur, rel=0.01)
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    # test whether the true period is 2P: compare depths of alternate eclipses
    best = None
    for k in (1, 2):
        Pk = P1 * k
        ph, prof, c = fold_profile(t, f, Pk, t01 - 0.25 * d1, nb=int(max(200, Pk / (d1 / 6))))
        res[f"profile_k{k}"] = None
        if k == 2:
            # depth at phase of the eclipse (bin containing t01) and half a period later
            ep = np.round((t - t01) / P1).astype(int)
            intr = np.abs(t - t01 - ep * P1) < 0.3 * d1
            d_even = 1 - np.mean(f[intr & (ep % 2 == 0)]); d_odd = 1 - np.mean(f[intr & (ep % 2 == 1)])
            res["alt_depths"] = (float(d_even), float(d_odd))
            ratio = min(d_even, d_odd) / max(d_even, d_odd)
            res["true_P"] = Pk if ratio < 0.85 else P1
    if eph and eph["rms_min"] < 30:
        res["true_P"] = P1  # timing ephemeris is built from one eclipse family, so P1 is already the orbital period
    Pt = res["true_P"]
    # primary/secondary at the true period
    ep = np.round((t - t01) / Pt).astype(int); dt = t - t01 - ep * Pt
    prim = np.abs(dt) < 0.3 * d1
    d_prim = 1 - np.mean(f[prim])
    # secondary search across phase (exclude primary)
    ph = ((t - t01) / Pt) % 1.0
    grid = np.linspace(0.05, 0.95, 361); sd = []
    for p0 in grid:
        m = np.abs(ph - p0) < 0.3 * d1 / Pt
        sd.append(1 - np.mean(f[m]) if m.sum() > 5 else np.nan)
    sd = np.array(sd); j = int(np.nanargmax(sd))
    oot = (np.abs(dt) > 1.5 * d1) & (np.abs(((ph - grid[j] + 0.5) % 1) - 0.5) * Pt > 1.5 * d1)
    noise = np.std(f[oot]) / np.sqrt(max(1, np.sum(np.abs(ph - grid[j]) < 0.3 * d1 / Pt)))
    res.update(P=float(Pt), t0=float(t01), dur_h=float(d1 * 24), depth_prim=float(d_prim), depth_sec=float(sd[j]),
               phase_sec=float(grid[j]), sec_sig=float(sd[j] / noise), oot_std=float(np.std(f[oot])))
    res["ecosw"] = float(np.pi / 2 * (grid[j] - 0.5)) if res["sec_sig"] > 5 else None
    # novelty
    try:
        res["vsx_etc_150"] = known.vizier_cone(star["ra"], star["dec"], 150)
        res["vizier_variability_15"], res["n_catalogs_15"] = known.vizier_all(star["ra"], star["dec"], 15)
        res["simbad_30"] = known.simbad(star["ra"], star["dec"], 30)
    except Exception as e:
        res["novelty_error"] = str(e)[:200]
    res["known_signals"] = known.known_signals(tic, Pt)[0]
    # localisation in the most recent FFI sector available
    loc = []
    for s in sorted(secs, key=lambda s: -s["sector"])[:2]:
        try:
            L = localise(star["ra"], star["dec"], s["sector"], Pt, t01, d1)
        except Exception:
            L = None
        if L: loc.append({k: v for k, v in L.items() if k not in ("diff_img", "direct_img", "wcs")})
    res["localisation"] = loc
    # plot
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    phf = ((t - t01) / Pt + 0.25) % 1.0 - 0.25
    ax[0].plot(phf, f, ",", color="0.6"); pb, prof, _ = fold_profile(t, f, Pt, t01 - 0.25 * Pt, nb=300)
    ax[0].plot(pb - 0.25, prof, "-", lw=0.8); ax[0].set_title(f"TIC {tic}  P={Pt:.6f} d (all sectors)", fontsize=9)
    for a, c0 in ((ax[1], 0.0), (ax[2], grid[j])):
        x = (((t - t01) / Pt - c0 + 0.5) % 1.0 - 0.5) * Pt * 24
        m = np.abs(x) < 4 * d1 * 24
        a.plot(x[m], f[m], ".", ms=1, color="0.6")
        hb = np.linspace(-4 * d1 * 24, 4 * d1 * 24, 60); idx = np.digitize(x[m], hb)
        bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if (idx == q).sum() > 3]
        bv = [f[m][idx == q].mean() for q in range(1, len(hb)) if (idx == q).sum() > 3]
        a.plot(bc, bv, "o-", ms=2, lw=0.8); a.set_xlabel("hours", fontsize=8)
    ax[1].set_title(f"primary, depth {d_prim * 100:.2f}%", fontsize=9)
    ax[2].set_title(f"secondary at phase {grid[j]:.3f}, depth {sd[j] * 100:.3f}% ({res['sec_sig']:.0f} sigma)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{outdir}/eb_fold.png", dpi=100); plt.close(fig)
    json.dump(res, open(f"{outdir}/eb.json", "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    return res


if __name__ == "__main__":
    tic = int(sys.argv[1]); sector = int(sys.argv[2])
    c = pd.read_csv(f"{W}/results/new_s{sector:04d}.csv"); r = c[c.tic == tic].iloc[0].to_dict()
    star = {k: (float(r[k]) if pd.notna(r[k]) else None) for k in ("Tmag", "Teff", "rad", "mass", "d", "ra", "dec", "contratio")}
    res = analyse(tic, r["P"], r["t0"], r["dur"], star, f"{W}/vet/EB_TIC{tic}")
    print(json.dumps({k: res.get(k) for k in ("P", "dur_h", "depth_prim", "depth_sec", "phase_sec", "sec_sig", "alt_depths", "known_signals", "vizier_variability_15", "simbad_30", "localisation")}, indent=1, default=str))
