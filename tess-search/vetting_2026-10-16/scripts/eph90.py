"""Clean per-sector ephemeris for TIC 90086266: unmasked biweight (1 d), fine BLS 2.2-2.5 d, then per-event times."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
import numpy as np
from wotan import flatten
from scipy.optimize import least_squares
import stage2
from bls import _bls, robust_std
sys.argv = [sys.argv[0], "90086266"]
from lcan import read, bin10, trap
TIC = 90086266
p = stage2.products(TIC)
data = {}
for _, r in p.iterrows():
    key = (int(r.sector), r.prov)
    path = stage2.fetch(r.uri); t, q, cols, bkg, aper, crowd = read(path, r.prov)
    c = "PDCSAP_FLUX" if "PDCSAP_FLUX" in cols else "SAP_FLUX"
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[c]) & (cols[c] > 0)
    tt, ff = bin10(t[ok], cols[c][ok] / np.nanmedian(cols[c][ok]), r.exptime)
    fl = flatten(tt, ff, method="biweight", window_length=1.0, break_tolerance=0.4)
    m = np.isfinite(fl); data[key] = (tt[m], fl[m], float(r.exptime))
res = {}
for key, (t, f, ex) in data.items():
    s = robust_std(f); y = f - 1; w = np.full(len(t), 1 / s**2); y -= np.sum(w*y)/np.sum(w)
    fr = 1 / np.linspace(2.2, 2.5, 6000)
    pw, dp, t0, du = _bls(t - t.min(), y, w, fr, np.array([2.5, 3.0, 3.5, 4.0, 4.5]) / 24, 5 / 1440)
    i = int(np.argmax(pw)); P = 1 / fr[i]; T0 = t0[i] + t.min() + du[i] / 2 - du[i] / 2
    # individual events: fit a trapezoid (shape fixed from the fold) for the time of each event
    n = np.round((t - T0) / P); tcs = []
    for k in np.unique(n):
        tc = T0 + k * P; m = np.abs(t - tc) < 0.3
        if m.sum() < 0.6 * 0.6 / (max(ex, 600) / 86400) or np.sum(np.abs(t[m] - tc) < 0.06) < 3: continue
        def r_(x): return (f[m] - trap(t[m], x[0], x[1], x[2], x[2] * 0.2)) / s
        rr = least_squares(r_, [tc, 0.016, 0.164], bounds=([tc - 0.1, 0, 0.08], [tc + 0.1, 0.06, 0.3]))
        J = rr.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * rr.cost / (m.sum() - 3))
        tcs.append((int(k), rr.x[0], np.sqrt(cov[0, 0]), rr.x[1]))
    a = np.array(tcs)
    A = np.column_stack([np.ones(len(a)), a[:, 0]]); W = 1 / a[:, 2]
    c, *_ = np.linalg.lstsq(A * W[:, None], a[:, 1] * W, rcond=None)
    cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None])); resid = a[:, 1] - A @ c
    chi2 = np.sum((resid / a[:, 2]) ** 2)
    res[f"S{key[0]}_{key[1]}"] = dict(bls_P=P, n_ev=len(a), T0=c[0], T0_err=np.sqrt(cov[0, 0]), P=c[1], P_err=np.sqrt(cov[1, 1]),
                                      chi2=chi2, dof=len(a) - 2, rms_min=np.std(resid) * 1440,
                                      depths=np.round(a[:, 3] * 1e6).tolist(), tc=a[:, 1].tolist(), tc_err=a[:, 2].tolist(), ep=a[:, 0].tolist())
    print(f"S{key[0]} {key[1]}: BLS P {P:.5f}; linear fit P {c[1]:.5f}±{np.sqrt(cov[1,1]):.5f} T0 {c[0]:.4f}±{np.sqrt(cov[0,0]):.4f} "
          f"({len(a)} events, chi2 {chi2:.1f}/{len(a)-2}, rms {np.std(resid)*1440:.1f} min); depths {np.round(a[:,3]*1e6).astype(int).tolist()}")
json.dump(res, open("out/90086266/eph.json", "w"), indent=1)
# cycle-count aliases linking S13 (TESS-SPOC) and S104 (QLP)
a, b = res["S13_TESS-SPOC"], res["S104_QLP"]
for N in range(1060, 1090):
    P = (b["T0"] - a["T0"]) / N
    if abs(P - a["P"]) < 4 * a["P_err"] + 0.002 or abs(P - b["P"]) < 4 * b["P_err"] + 0.002:
        print(f"  N={N}: P={P:.6f}  (S13 {(P-a['P'])/a['P_err']:+.1f} sigma, S104 {(P-b['P'])/b['P_err']:+.1f} sigma)")
