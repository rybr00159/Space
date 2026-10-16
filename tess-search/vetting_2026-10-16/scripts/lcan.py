"""Light-curve vetting for one signal: all products and apertures, ephemeris, depths per sector and product,
odd/even, secondary scan, phase curve, transit-shape fit and implied stellar density, per-event depths.
Usage: python lcan.py TIC [P_override]"""
import os, sys, json, warnings
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/user/work/src")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astropy.io import fits
from wotan import flatten
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import stage2
from bls import _bls, robust_std
from vet import model as bmodel, rho_star

HERE = os.path.dirname(os.path.abspath(__file__))
TIC = int(sys.argv[1])
seed = json.load(open(f"{HERE}/seed.json"))[str(TIC)]
P0, T00, DUR0 = seed["P"], seed["t0"], seed["dur"]
if len(sys.argv) > 2: P0 = float(sys.argv[2])
if len(sys.argv) > 3: T00 = float(sys.argv[3])
if len(sys.argv) > 4: DUR0 = float(sys.argv[4])
FIXED = len(sys.argv) > 2 and os.environ.get('LCAN_FIX', '0') == '1'
stars = pd.read_csv(f"{HERE}/stars.csv"); st = stars[stars.ID == TIC].iloc[0]
OUT = f"{HERE}/out/{TIC}"; os.makedirs(OUT, exist_ok=True)


def products():
    p = stage2.products(TIC)
    return p.sort_values(["sector", "prov", "exptime"])


def read(path, prov):
    """Return dict of flux columns (normalised later), time, quality-good mask, background, centroids."""
    with fits.open(path) as f:
        d = f[1].data; names = d.columns.names; h1 = f[1].header
        t = np.asarray(d["TIME"], float); q = np.asarray(d["QUALITY"])
        cols = {}
        for c in ("SAP_FLUX", "PDCSAP_FLUX", "DET_FLUX", "KSPSAP_FLUX"):
            if c in names: cols[c] = np.asarray(d[c], float)
        bkg = np.asarray(d["SAP_BKG"], float) if "SAP_BKG" in names else None
        aper = h1.get("APERTURE", h1.get("APER", ""))
        crowd = f[1].header.get("CROWDSAP", None)
    return t, q, cols, bkg, aper, crowd


def bin10(t, f, exptime):
    if exptime > 500: return t, f
    k = np.floor((t - t[0]) / (600 / 86400)).astype(int)
    _, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    tb = np.bincount(inv, t) / cnt; fb = np.bincount(inv, f) / cnt
    return tb[cnt >= 2], fb[cnt >= 2]


def in_mask(t, P, t0, dur, k=1.0):
    return np.abs(((t - t0) / P + 0.5) % 1 - 0.5) * P < k * dur


def detrend(t, f, P, t0, dur, win=None):
    win = win or max(0.6, 3.5 * dur)
    m = ~in_mask(t, P, t0, dur, 0.75)
    tr = flatten(t[m], f[m], method="biweight", window_length=win, return_trend=True, break_tolerance=0.4)[1]
    ok = np.isfinite(tr)
    trend = np.interp(t, t[m][ok], tr[ok])
    return f / trend


def load_all(P, t0, dur):
    secs = []
    for _, r in products().iterrows():
        path = stage2.fetch(r.uri)
        if not path: continue
        t, q, cols, bkg, aper, crowd = read(path, r.prov)
        for c, fl in cols.items():
            ok = (q == 0) & np.isfinite(t) & np.isfinite(fl) & (fl > 0)
            if ok.sum() < 200: continue
            tt, ff = bin10(t[ok], fl[ok] / np.nanmedian(fl[ok]), r.exptime)
            # split at gaps > 1 d and detrend per segment
            seg = np.concatenate([[0], np.cumsum(np.diff(tt) > 1.0)])
            fd = np.concatenate([detrend(tt[seg == s], ff[seg == s], P, t0, dur) for s in np.unique(seg)])
            oot = ~in_mask(tt, P, t0, dur, 1.0)
            sig = robust_std(fd[oot])
            keep = fd < 1 + 5 * sig  # upward outliers only
            secs.append(dict(sector=int(r.sector), prov=r.prov, col=c, exptime=float(r.exptime), t=tt[keep], f=fd[keep],
                             sigma=sig, aper=str(aper), crowd=crowd))
    return secs


def best_per_sector(secs):
    order = {("SPOC", "PDCSAP_FLUX"): 0, ("TESS-SPOC", "PDCSAP_FLUX"): 1, ("QLP", "SAP_FLUX"): 2}
    out = {}
    for s in secs:
        k = order.get((s["prov"], s["col"]))
        if k is None: continue
        if s["sector"] not in out or k < out[s["sector"]][0]: out[s["sector"]] = (k, s)
    return [v[1] for k, v in sorted(out.items())]


def box_depth(t, f, sigma, P, t0, dur, phase=0.0, half=0.5):
    ph = ((t - t0) / P - phase + 0.5) % 1 - 0.5
    inn = np.abs(ph * P) < half * dur
    oot = (np.abs(ph * P) > 1.0 * dur) & (np.abs(ph * P) < 3.0 * dur)
    if inn.sum() < 3 or oot.sum() < 10: return np.nan, np.nan, int(inn.sum())
    d = np.median(f[oot]) - np.mean(f[inn])
    e = sigma * np.sqrt(1 / inn.sum() + 1 / oot.sum())
    return d, e, int(inn.sum())


def refine(secs, P, dur):
    """Fine BLS on the combined light curve around P; also per-sector periods."""
    j = stage2.joint_fit(secs, P, dur, rel=0.004)
    per = []
    for s in secs:
        if s["t"].max() - s["t"].min() < 2.2 * P: continue
        jj = stage2.joint_fit([s], P, dur, rel=0.02)
        per.append(dict(sector=s["sector"], prov=s["prov"], P=jj["P"], snr=jj["snr"]))
    return j, per


def trap(t, t0, depth, T14, T23):
    a = np.abs(t - t0)
    T23 = min(T23, T14 * 0.999)
    y = np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
    return 1 - depth * y


def shape_fit(secs, P, t0, dur):
    tt, ff, ww = [], [], []
    for s in secs:
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5
        m = np.abs(ph * P) < 3 * dur
        tt.append(ph[m] * P); ff.append(s["f"][m]); ww.append(np.full(m.sum(), 1 / s["sigma"]))
    x = np.concatenate(tt); y = np.concatenate(ff); w = np.concatenate(ww)
    def res(p):
        return (y - trap(x, p[0], p[1], p[2], p[2] * p[3])) * w
    best = None
    for fr in (0.1, 0.5, 0.8):
        r = least_squares(res, [0.0, max(1e-4, 1 - np.median(y[np.abs(x) < dur / 3])), dur, fr],
                          bounds=([-dur / 2, 0, dur / 4, 0.0], [dur / 2, 0.5, 3 * dur, 1.0]))
        if best is None or r.cost < best.cost: best = r
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * best.cost / max(1, len(y) - 4)); e = np.sqrt(np.diag(cov))
    # chi2 profile of the flat-bottom fraction: V shape (T23/T14 = 0) vs best
    def prof(fr):
        r = least_squares(lambda p: (y - trap(x, p[0], p[1], p[2], p[2] * fr)) * w, best.x[:3],
                          bounds=([-dur / 2, 0, dur / 4], [dur / 2, 0.5, 3 * dur]))
        return 2 * r.cost
    return dict(dt0=best.x[0], depth=best.x[1], depth_err=e[1], T14_h=best.x[2] * 24, T14_err_h=e[2] * 24,
                flat_frac=best.x[3], flat_frac_err=e[3], dchi2_V=prof(0.0) - 2 * best.cost, chi2=2 * best.cost, n=len(y))


def batman_fit(secs, P, t0, dur, depth):
    w = []
    for s in secs:
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5
        m = np.abs(ph * P) < 3 * dur
        if m.sum() > 10: w.append(dict(s, t=s["t"][m], f=s["f"][m]))
    def resid(x):
        t0_, rp, la, b = x
        return np.concatenate([(s["f"] - bmodel(s["t"], P, t0_, rp, 10 ** la, b, max(600.0, s["exptime"]))) / s["sigma"] for s in w])
    a0 = max(2.0, P / (np.pi * dur)); best = None
    for b0 in (0.1, 0.5, 0.85, 1.0):
        r = least_squares(resid, [t0, np.sqrt(max(depth, 1e-5)), np.log10(a0), b0],
                          bounds=([t0 - dur / 2, 1e-3, 0.0, 0.0], [t0 + dur / 2, 0.6, 2.5, 1.5]), x_scale="jac")
        if best is None or r.cost < best.cost: best = r
    J = best.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * best.cost / max(1, len(best.fun) - 4)); e = np.sqrt(np.diag(cov))
    t0_, rp, la, b = best.x
    return dict(t0=t0_, rp=rp, rp_err=e[1], aR=10 ** la, aR_err=10 ** la * np.log(10) * e[2], b=b, b_err=e[3],
                rho=rho_star(10 ** la, P), chi2=2 * best.cost, n=len(best.fun))


def secondary_scan(secs, P, t0, dur):
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    sg = np.concatenate([np.full(len(s["t"]), s["sigma"]) for s in secs])
    out = []
    for ph0 in np.arange(0.06, 0.94 + 1e-9, max(0.005, dur / P / 4)):
        ph = ((t - t0) / P - ph0 + 0.5) % 1 - 0.5
        inn = np.abs(ph * P) < 0.5 * dur
        if inn.sum() < 5: continue
        wv = 1 / sg[inn] ** 2
        d = np.sum(wv * (1 - f[inn])) / wv.sum(); e = 1 / np.sqrt(wv.sum())
        out.append((ph0, d, e))
    a = np.array(out)
    i = int(np.argmax(a[:, 1] / a[:, 2])); j = int(np.argmin(np.abs(a[:, 0] - 0.5)))
    snr = a[:, 1] / a[:, 2]
    return dict(max_phase=a[i, 0], max_depth=a[i, 1], max_err=a[i, 2], max_snr=a[i, 1] / a[i, 2],
                at05_depth=a[j, 1], at05_err=a[j, 2], snr_rms=float(np.std(snr)), n_trials=len(a)), a


def phase_curve(secs, P, t0, dur):
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    sg = np.concatenate([np.full(len(s["t"]), s["sigma"]) for s in secs])
    ph = ((t - t0) / P) % 1
    m = ~in_mask(t, P, t0, dur, 1.0) & ~in_mask(t, P, t0 + P / 2, dur, 1.0)
    A = np.column_stack([np.ones(m.sum()), np.cos(2 * np.pi * ph[m]), np.sin(2 * np.pi * ph[m]), np.cos(4 * np.pi * ph[m])])
    W = 1 / sg[m]
    c, *_ = np.linalg.lstsq(A * W[:, None], (f[m] - 1) * W, rcond=None)
    cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None])); e = np.sqrt(np.diag(cov))
    return dict(cos1=c[1], cos1_err=e[1], sin1=c[2], sin1_err=e[2], cos2=c[3], cos2_err=e[3])


def events(secs, P, t0, dur):
    out = []
    for s in secs:
        n = np.round((s["t"] - t0) / P)
        for k in np.unique(n):
            tc = t0 + k * P
            m = np.abs(s["t"] - tc) < 3 * dur
            if m.sum() < 8: continue
            d, e, nin = box_depth(s["t"][m], s["f"][m], s["sigma"], P, t0, dur)
            if np.isfinite(d) and nin >= max(2, 0.4 * dur / (max(600, s["exptime"]) / 86400)):
                out.append(dict(sector=s["sector"], epoch=int(k), tc=tc, depth=d, err=e))
    return pd.DataFrame(out)


if __name__ == "__main__":
    allsecs = load_all(P0, T00, DUR0)
    best = best_per_sector(allsecs)
    j, per = refine(best, P0, DUR0)
    P, t0, dur = j["P"], j["t0"], j["dur"]
    if FIXED: P, t0, dur = P0, T00, DUR0
    # move t0 to the middle of the data span to minimise P-t0 covariance
    tm = np.median(np.concatenate([s["t"] for s in best])); t0 = t0 + np.round((tm - t0) / P) * P
    res = dict(tic=TIC, P=P, t0=t0, dur_h=dur * 24, joint_snr=j["snr"], second_alias_snr=j["second"], per_sector_P=per)
    # depths per sector and product (apertures differ between products)
    tab = []
    for s in allsecs:
        d, e, nin = box_depth(s["t"], s["f"], s["sigma"], P, t0, dur)
        dd, de, _ = box_depth(s["t"], s["f"], s["sigma"], 2 * P, t0 + P, dur)  # odd events only (alternate)
        tab.append(dict(sector=s["sector"], prov=s["prov"], col=s["col"], aper=s["aper"][:24], crowd=s["crowd"],
                        depth_ppm=round(d * 1e6), err_ppm=round(e * 1e6), snr=round(d / e, 1) if e == e and e > 0 else np.nan))
    res["depths"] = tab
    # odd / even on best products
    t = np.concatenate([s["t"] for s in best]); f = np.concatenate([s["f"] for s in best])
    sg = np.concatenate([np.full(len(s["t"]), s["sigma"]) for s in best])
    de, ee, _ = box_depth(t, f, np.median(sg), 2 * P, t0, dur)
    do, eo, _ = box_depth(t, f, np.median(sg), 2 * P, t0 + P, dur)
    res["even"] = (de, ee); res["odd"] = (do, eo); res["oddeven_sigma"] = abs(de - do) / np.hypot(ee, eo)
    res["shape"] = shape_fit(best, P, t0, dur)
    res["batman"] = batman_fit(best, P, t0, dur, res["shape"]["depth"])
    sec, scan = secondary_scan(best, P, t0, dur); res["secondary"] = sec
    res["phase_curve"] = phase_curve(best, P, t0, dur)
    ev = events(best, P, t0, dur); res["n_events"] = len(ev)
    if len(ev):
        w = 1 / ev.err ** 2; mu = np.sum(w * ev.depth) / w.sum()
        res["event_chi2"] = float(np.sum(((ev.depth - mu) / ev.err) ** 2)); res["event_dof"] = len(ev) - 1
        ev.to_csv(f"{OUT}/events.csv", index=False)
    # stellar density for comparison
    R = st.radius_flame if np.isfinite(st.radius_flame) else st.rad; M = st.mass_flame if np.isfinite(st.mass_flame) else st.mass
    res["rho_star_catalogue"] = M / R ** 3 if np.isfinite(M) and np.isfinite(R) else np.nan
    json.dump(res, open(f"{OUT}/lcan.json", "w"), indent=1, default=float)
    # summary
    print(f"TIC {TIC}: P={P:.6f} d t0={t0:.4f} dur={dur*24:.2f} h joint S/N {j['snr']:.1f} (next alias {j['second']:.1f})")
    print("  per-sector P:", [(x['sector'], round(x['P'], 5), round(x['snr'], 1)) for x in per])
    print(pd.DataFrame(tab).to_string(index=False))
    print(f"  even {de*1e6:.0f}±{ee*1e6:.0f}  odd {do*1e6:.0f}±{eo*1e6:.0f}  ({res['oddeven_sigma']:.1f} sigma)")
    sh = res["shape"]; bm = res["batman"]
    print(f"  trapezoid: depth {sh['depth']*1e6:.0f}±{sh['depth_err']*1e6:.0f} ppm T14 {sh['T14_h']:.2f}±{sh['T14_err_h']:.2f} h flat fraction {sh['flat_frac']:.2f}±{sh['flat_frac_err']:.2f}  dchi2(V-shape) {sh['dchi2_V']:.1f}")
    print(f"  batman: Rp/R* {bm['rp']:.4f}±{bm['rp_err']:.4f} a/R* {bm['aR']:.2f}±{bm['aR_err']:.2f} b {bm['b']:.2f}±{bm['b_err']:.2f} rho {bm['rho']:.3f} (catalogue {res['rho_star_catalogue']:.3f}) rho_sun")
    s_ = res["secondary"]
    print(f"  secondary: at 0.5 {s_['at05_depth']*1e6:.0f}±{s_['at05_err']*1e6:.0f} ppm; strongest at phase {s_['max_phase']:.3f}: {s_['max_depth']*1e6:.0f}±{s_['max_err']*1e6:.0f} ppm ({s_['max_snr']:.1f} sigma; rms of scan S/N {s_['snr_rms']:.2f}, {s_['n_trials']} trials)")
    pc = res["phase_curve"]
    print(f"  phase curve (ppm): cos1 {pc['cos1']*1e6:.0f}±{pc['cos1_err']*1e6:.0f} sin1 {pc['sin1']*1e6:.0f}±{pc['sin1_err']*1e6:.0f} cos2 {pc['cos2']*1e6:.0f}±{pc['cos2_err']*1e6:.0f}")
    if len(ev): print(f"  events: {len(ev)}; depth chi2 {res['event_chi2']:.1f} for {res['event_dof']} dof")
    # figure
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for i, s in enumerate(best):
        ph = ((s["t"] - t0) / P + 0.5) % 1 - 0.5; m = np.abs(ph * P) < 3 * dur
        hb = np.linspace(-3 * dur * 24, 3 * dur * 24, 37); x = ph[m] * P * 24; idx = np.digitize(x, hb)
        bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if np.sum(idx == q) > 2]
        bv = [np.mean(s["f"][m][idx == q] - 1) * 1e6 for q in range(1, len(hb)) if np.sum(idx == q) > 2]
        ax[0, 0].plot(bc, np.array(bv) - 1.5 * i * max(300, sh["depth"] * 1e6), "o-", ms=3, lw=0.7, label=f"S{s['sector']} {s['prov']}")
    ax[0, 0].set_title("transit per sector (offset)"); ax[0, 0].legend(fontsize=6); ax[0, 0].set_xlabel("hours")
    for k, (lab, off) in enumerate((("even", 0), ("odd", P))):
        ph = ((t - t0 - off) / (2 * P) + 0.5) % 1 - 0.5; m = np.abs(ph * 2 * P) < 3 * dur
        hb = np.linspace(-3 * dur * 24, 3 * dur * 24, 37); x = ph[m] * 2 * P * 24; idx = np.digitize(x, hb)
        bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if np.sum(idx == q) > 2]
        bv = [np.mean(f[m][idx == q] - 1) * 1e6 for q in range(1, len(hb)) if np.sum(idx == q) > 2]
        ax[0, 1].plot(bc, bv, "o-", ms=3, lw=0.7, label=lab)
    ax[0, 1].legend(); ax[0, 1].set_title(f"odd/even ({res['oddeven_sigma']:.1f} sigma)")
    ax[0, 2].errorbar(scan[:, 0], scan[:, 1] * 1e6, scan[:, 2] * 1e6, fmt=".", ms=2); ax[0, 2].axhline(0, color="k", lw=0.5)
    ax[0, 2].set_title("secondary scan (box depth vs phase)")
    ph = ((t - t0) / P) % 1; hb = np.linspace(0, 1, 101); idx = np.digitize(ph, hb)
    bv = [np.mean(f[idx == q] - 1) * 1e6 for q in range(1, len(hb))]
    ax[1, 0].plot((hb[:-1] + hb[1:]) / 2, bv, "k.-", lw=0.5); ax[1, 0].set_title("full phase, 1% bins")
    lim = max(5 * sh["depth"] * 1e6 / 50, 3 * np.nanstd(bv[20:80]))
    ax[1, 1].plot((hb[:-1] + hb[1:]) / 2, bv, "k.-", lw=0.5); ax[1, 1].set_ylim(-lim, lim); ax[1, 1].set_title("full phase, zoomed")
    if len(ev):
        ax[1, 2].errorbar(ev.tc, ev.depth * 1e6, ev.err * 1e6, fmt="o", ms=3); ax[1, 2].set_title("per-event depth"); ax[1, 2].set_xlabel("BTJD")
    fig.suptitle(f"TIC {TIC}  P={P:.6f} d"); fig.tight_layout(); fig.savefig(f"{OUT}/lcan.png", dpi=90); plt.close(fig)
