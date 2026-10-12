"""Secondary eclipse, odd/even and phase-curve (ellipsoidal / beaming / reflection) tests on the PSF PAIR light curve
(F_target + F_neighbour from psf_photometry/lc/nominal; data/lc_pair.npy from s02_lc.py), converted to the
neighbour's own flux scale with the dilution D = F_t/F_n = 1/r (r = 0.1879 +- 0.0015, results/stellar.json).
The pair light curve does not depend on the target/neighbour flux split (registration-degenerate), and is 6-15x less
noisy than the neighbour-only PSF light curve, so limits are tighter than psf_photometry's neighbour-only numbers.

A. Re-detrending: wotan biweight (1.0 d) of f_raw with BOTH the transit (phase 0) and the secondary (phase 0.5)
   masked (0.30 d full width), segments as in s02.
B. Secondary at phase 0.5: template = occultation shape (batman, u = 0) of the best neighbour-host transit model
   (k, b, a/R* from mcmc/mcmc_pair_local.json); amplitude via nbtools.depth_with_errors (per-epoch quadratic
   baselines, error = max(white, epoch scatter, 50 control phases)).  Also a scan over phases 0.30-0.70.
C. Odd/even: transit template from the same fit; odd and even epochs separately (absolute epoch parity from T0).
D. Phase curve: f_raw per segment = quadratic polynomial x (1 + A_ell cos(4 pi phi) + A_beam sin(2 pi phi)
   + A_refl cos(2 pi phi)), eclipses masked (|phase| < 0.04 around 0 and 0.5), weighted LSQ with per-sector OOT rms.
   Errors: scatter of the same fit at 60 control periods (2.6-7.0 d, avoiding P, P/2, 2P by 5 %); recovery
   factor from injection of A_ell = 1.0 ppt (neighbour scale) at the true period.
Output: results/pair_tests.json, fig/pair_secondary_scan.png, fig/pair_phasecurve.png
"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbtools as nt
from wotan import flatten
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = nt.ROOT
st = json.load(open(os.path.join(ROOT, "results/stellar.json")))["adopted"]
D, D_err = st["D_pair"], st["D_pair_err"]
SC = 1 + D  # pair depth -> neighbour depth
mc = json.load(open(os.path.join(ROOT, "mcmc/mcmc_pair_local.json")))
pb = mc["params"]
P, T0 = pb["P"]["med"], pb["T0"]["med"]
k, b, a = pb["k"]["med"], pb["b"]["med"], pb["a"]["med"]
SHAPE = dict(k=k, b=b, a=a, u1=pb["u1"]["med"], u2=pb["u2"]["med"])
SHAPE_OCC = dict(k=k, b=b, a=a, u1=1e-4, u2=0.0)
nt.P, nt.T0 = P, T0
lc = nt.load_lc("pair")
res = dict(D=D, D_err=D_err, ephemeris=dict(P=P, T0=T0), shape=SHAPE)


def redetrend(v, win=1.0, maskw=0.30):
    t, f, seg = v["t"], v["f_raw"], v["seg"]
    ph = ((t - T0) / P) % 1
    mask = (np.minimum(ph, 1 - ph) * P < maskw / 2) | (np.abs(ph - 0.5) * P < maskw / 2)
    tr = np.full_like(f, np.nan)
    for s_ in np.unique(seg):
        m = seg == s_
        ff = f[m].copy(); ff[mask[m]] = np.nan
        if np.isfinite(ff).sum() < 5: continue
        x = flatten(t[m], ff, method="biweight", window_length=win, return_trend=True, edge_cutoff=0)[1]
        g = np.isfinite(x)
        if g.sum() >= 2: tr[m] = np.interp(t[m], t[m][g], x[g])
    ok = np.isfinite(tr)
    return t[ok], (f / tr)[ok]


det = {s: redetrend(lc[s]) for s in nt.SECTORS}


def amp(sector_list, phase, shape, epochs_parity=None):
    out = {}
    for s in sector_list:
        t, f = det[s]
        kw = dict(shape=shape)
        if epochs_parity is not None:
            ep = np.unique(np.round((t - (T0 + phase * P)) / P))
            ep_abs = ep + round((T0 - 1658.41325) / P)
            kw["epochs"] = [e for e, ea in zip(ep, ep_abs) if int(ea) % 2 == epochs_parity]
        r = nt.depth_with_errors(t, f, lc[s]["exptime"], phase=phase, ncontrol=50, **kw)
        if r.get("n", 0) > 0 and np.isfinite(r.get("err", np.nan)):
            out[s] = dict(A=r["A"], err=r["err"], n=r["n"])
    if not out: return dict(A=np.nan, err=np.nan)
    m, e, c2 = nt.wmean([v["A"] for v in out.values()], [v["err"] for v in out.values()])
    return dict(A=m, err=e, chi2=c2, dof=len(out) - 1, per_sector=out)


# ---------------------------------------------------------------- B: secondary
sec = amp(nt.SECTORS, 0.5, SHAPE_OCC)
pri = amp(nt.SECTORS, 0.0, SHAPE)
res["primary_pair"] = pri
res["secondary_pair"] = sec
res["secondary_nb"] = dict(depth_ppt=sec["A"] * SC * 1e3, err_ppt=sec["err"] * SC * 1e3,
                           upper2sig_ppt=(max(sec["A"], 0) + 2 * sec["err"]) * SC * 1e3,
                           upper2sig_ppt_conservative=(abs(sec["A"]) + 2 * sec["err"]) * SC * 1e3)
res["primary_nb"] = dict(depth_ppt=pri["A"] * SC * 1e3, err_ppt=pri["err"] * SC * 1e3)
print("primary", res["primary_nb"], "secondary", res["secondary_nb"], flush=True)
scan = []
for ph in np.arange(0.30, 0.7001, 0.01):
    r = amp(nt.SECTORS, ph, SHAPE_OCC)
    scan.append((float(ph), r["A"] * SC * 1e3, r["err"] * SC * 1e3))
scan = np.array(scan)
res["secondary_scan"] = dict(phase=scan[:, 0].tolist(), depth_ppt=scan[:, 1].tolist(), err_ppt=scan[:, 2].tolist(),
                             max_abs_snr=float(np.max(np.abs(scan[:, 1] / scan[:, 2]))),
                             phase_max=float(scan[np.argmax(np.abs(scan[:, 1] / scan[:, 2])), 0]),
                             upper2sig_any_ppt=float(np.max(np.maximum(scan[:, 1], 0) + 2 * scan[:, 2])))
print("scan max |S/N|", res["secondary_scan"]["max_abs_snr"], "at", res["secondary_scan"]["phase_max"], flush=True)

# ---------------------------------------------------------------- C: odd / even
odd = amp(nt.SECTORS, 0.0, SHAPE, epochs_parity=1)
even = amp(nt.SECTORS, 0.0, SHAPE, epochs_parity=0)
dd = odd["A"] - even["A"]; de = np.hypot(odd["err"], even["err"])
res["oddeven_nb"] = dict(odd_ppt=odd["A"] * SC * 1e3, odd_err=odd["err"] * SC * 1e3, even_ppt=even["A"] * SC * 1e3,
                         even_err=even["err"] * SC * 1e3, diff_ppt=dd * SC * 1e3, diff_err=de * SC * 1e3,
                         sigma=float(dd / de), ratio=float(odd["A"] / even["A"]))
print("odd/even", res["oddeven_nb"], flush=True)


# ---------------------------------------------------------------- D: phase curve
def pc_fit(Pc, T0c, inject=0.0, Ptrue=P, T0true=T0):
    rows_A, rows_y, rows_w = [], [], []
    ncol_seg = 0
    blocks = []
    for s in nt.SECTORS:
        v = lc[s]; t = v["t"]; y = v["f_raw"].copy()
        if inject:
            y = y * (1 - inject * np.cos(4 * np.pi * (t - T0true) / Ptrue))  # ellipsoidal convention: minima at phases 0, 0.5
        ph = ((t - T0true) / Ptrue) % 1
        ecl = (np.minimum(ph, 1 - ph) < 0.04) | (np.abs(ph - 0.5) < 0.04)
        # bin to 30 min to speed up; the per-segment quadratic is fitted on the same bins
        sig = v["sigma"]
        for sg in np.unique(v["seg"]):
            m = (v["seg"] == sg) & ~ecl
            if m.sum() < 20 or np.ptp(t[m]) < 1.0: continue
            tb = np.floor((t[m] - t[m][0]) / (1 / 48))
            u_, inv = np.unique(tb, return_inverse=True)
            n_ = np.bincount(inv)
            tt = np.bincount(inv, t[m]) / n_; yy = np.bincount(inv, y[m]) / n_
            blocks.append((tt, yy, sig / np.sqrt(n_), np.mean(t[m])))
    nseg = len(blocks)
    N = sum(len(b_[0]) for b_ in blocks)
    A = np.zeros((N, 3 * nseg + 3)); yv = np.zeros(N); wv = np.zeros(N)
    i0 = 0
    for j, (tt, yy, ee, tm) in enumerate(blocks):
        n = len(tt); x = (tt - tm)
        A[i0:i0 + n, 3 * j] = 1; A[i0:i0 + n, 3 * j + 1] = x; A[i0:i0 + n, 3 * j + 2] = x ** 2
        phc = (tt - T0c) / Pc
        A[i0:i0 + n, -3] = np.cos(4 * np.pi * phc); A[i0:i0 + n, -2] = np.sin(2 * np.pi * phc)
        A[i0:i0 + n, -1] = np.cos(2 * np.pi * phc)
        yv[i0:i0 + n] = yy; wv[i0:i0 + n] = 1 / ee
        i0 += n
    # the harmonic terms multiply the local baseline ~1, so fitting them additively is equivalent to first order
    c, *_ = np.linalg.lstsq(A * wv[:, None], yv * wv, rcond=None)
    base = A[:, :-3] @ c[:-3]
    return dict(A_ell=-c[-3] / np.mean(base), A_beam=c[-2] / np.mean(base), A_refl=-c[-1] / np.mean(base), nseg=nseg)


real = pc_fit(P, T0)
rng = np.random.default_rng(3)
ctrl = []
for Pc in np.exp(rng.uniform(np.log(2.6), np.log(7.0), 400)):
    if any(abs(Pc / q - 1) < 0.05 for q in (P, P / 2, 2 * P, 2 * P / 3, 1.5 * P)): continue
    ctrl.append(pc_fit(Pc, T0 + rng.uniform(0, Pc)))
    if len(ctrl) >= 60: break
ce = {kk: float(np.std([c_[kk] for c_ in ctrl])) for kk in ("A_ell", "A_beam", "A_refl")}
cm = {kk: float(np.mean([c_[kk] for c_ in ctrl])) for kk in ("A_ell", "A_beam", "A_refl")}
inj_amp = 1.0e-3 / SC  # 1 ppt on the neighbour
inj = pc_fit(P, T0, inject=inj_amp)
rec = float((inj["A_ell"] - real["A_ell"]) / inj_amp)
pc = {}
for kk in ("A_ell", "A_beam", "A_refl"):
    val = real[kk] / rec * SC * 1e3; err = ce[kk] / rec * SC * 1e3
    pc[kk] = dict(ppt_nb=float(val), err_ppt_nb=float(err), upper2sig_ppt_nb=float(abs(val) + 2 * err),
                  pair_ppm=float(real[kk] * 1e6), ctrl_sd_pair_ppm=float(ce[kk] * 1e6), ctrl_mean_pair_ppm=float(cm[kk] * 1e6))
res["phase_curve"] = dict(fit=pc, recovery_factor=rec, n_control=len(ctrl), nseg=real["nseg"],
                          note="A_ell: amplitude of -cos(4 pi phi) (minima at phases 0 and 0.5); A_beam: sin(2 pi phi) "
                               "(maximum at phase 0.25); A_refl: -cos(2 pi phi) (maximum at phase 0.5). "
                               "Neighbour-scale values = pair value / recovery x (1 + D).")
print("phase curve", json.dumps(res["phase_curve"], default=float), flush=True)
json.dump(res, open(os.path.join(ROOT, "results/pair_tests.json"), "w"), indent=1, default=float)

# ---------------------------------------------------------------- figures
fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))
ax.errorbar(scan[:, 0], scan[:, 1], yerr=scan[:, 2], fmt="o", ms=3)
ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("phase"); ax.set_ylabel("occultation depth on neighbour (ppt)")
ax.set_title("Secondary-eclipse scan, PSF pair light curve x (1 + D)")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/pair_secondary_scan.png"), dpi=110)
fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))
allp, ally = [], []
for s in nt.SECTORS:
    t, f = det[s]
    allp.append(((t - T0) / P) % 1); ally.append((f - 1) * SC * 1e3)
allp = np.concatenate(allp); ally = np.concatenate(ally)
be = np.linspace(0, 1, 101); idx = np.digitize(allp, be)
bm = [np.median(ally[idx == i]) if np.sum(idx == i) > 10 else np.nan for i in range(1, 101)]
ax.plot(0.5 * (be[1:] + be[:-1]), bm, "o", ms=3)
ax.set_ylim(-30, 8); ax.set_xlabel("phase"); ax.set_ylabel("neighbour-scale flux - 1 (ppt)")
ax.set_title("Pair LC (re-detrended, both eclipses masked), folded, 1 % phase bins")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/pair_phasefold.png"), dpi=110)
print("DONE")
