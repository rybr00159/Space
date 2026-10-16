"""Per-event times of the 5.69 d dips (TIC 214311906 field) from the QLP light curves of the M-dwarf pair TIC 214311910
(highest S/N) and of the Mira; linear ephemeris and predicted dip times for 2027."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "214311906"]
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from wotan import flatten
from astropy.time import Time
import stage2
from lcan import read, bin10, in_mask, trap
P0, T0, T14, FL = 5.69149, 3840.913, 0.067, 0.78
out = {}
for tic in (214311910, 214311906):
    p = stage2.products(tic); ev = []
    for _, r in p[p.prov == "QLP"].sort_values("sector").iterrows():
        t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), "QLP")
        col = "SAP_FLUX"
        ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[col]) & (cols[col] > 0)
        tt, ff = bin10(t[ok], cols[col][ok] / np.nanmedian(cols[col][ok]), r.exptime)
        m = ~in_mask(tt, P0, T0, T14, 1.5)
        tr = flatten(tt[m], ff[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        k_ = np.isfinite(tr); f = ff / np.interp(tt, tt[m][k_], tr[k_])
        s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
        for k in np.unique(np.round((tt - T0) / P0)):
            tc = T0 + k * P0; w = np.abs(tt - tc) < 0.3
            if np.sum(np.abs(tt[w] - tc) < T14) < 4 or w.sum() < 15: continue
            grid = np.arange(-0.15, 0.1501, 0.003); c2 = []
            for dt in grid:
                y = 1 - trap(tt[w], tc + dt, 1.0, T14, T14 * FL); A = np.column_stack([np.ones(w.sum()), y])
                c, *_ = np.linalg.lstsq(A, f[w], rcond=None); c2.append(np.sum(((f[w] - A @ c) / s) ** 2))
            j = int(np.argmin(c2))
            fit = least_squares(lambda x: (f[w] - x[1] + x[2] * (1 - trap(tt[w], tc + x[0], 1.0, T14, T14 * FL))) / s,
                                [grid[j], 1.0, 0.01], bounds=([-0.16, 0.95, -0.2], [0.16, 1.05, 0.2]))
            J = fit.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * fit.cost / (w.sum() - 3))
            if fit.x[2] / np.sqrt(cov[2, 2]) > 3:
                ev.append(dict(sector=int(r.sector), k=int(k), tc=tc + fit.x[0], err=float(np.sqrt(cov[0, 0])), depth=fit.x[2], snr=fit.x[2] / np.sqrt(cov[2, 2])))
    ev = pd.DataFrame(ev)
    W = 1 / ev.err.values; A = np.column_stack([np.ones(len(ev)), ev.k.values])
    c, *_ = np.linalg.lstsq(A * W[:, None], ev.tc.values * W, rcond=None); cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None]))
    res = ev.tc.values - A @ c; chi2 = float(np.sum((res * W) ** 2)); dof = len(ev) - 2
    sc = max(1.0, chi2 / dof); cov = cov * sc
    out[tic] = dict(P=c[1], P_err=float(np.sqrt(cov[1, 1])), T0=c[0], T0_err=float(np.sqrt(cov[0, 0])), cov=cov.tolist(), chi2=chi2, dof=dof, n=len(ev),
                    events=ev.to_dict("records"))
    print(f"TIC {tic}: {len(ev)} events at S/N>3; P = {c[1]:.6f} ± {np.sqrt(cov[1,1]):.6f} d, T0 = {c[0]:.4f} ± {np.sqrt(cov[0,0]):.4f} BTJD; chi2 {chi2:.1f}/{dof}; "
          f"O-C rms {np.std(res)*1440:.1f} min; per-sector events {ev.groupby('sector').size().to_dict()}")
# predicted 2027 windows from the pair's ephemeris (epochs relative to T0)
e = out[214311910]; Pp, Tp = e["P"], e["T0"]; C = np.array(e["cov"])
start = Time("2027-03-01").jd - 2457000; stop = Time("2027-09-30").jd - 2457000
ks = np.arange(np.ceil((start - Tp) / Pp), np.floor((stop - Tp) / Pp) + 1)
pred = []
for k in ks:
    tc = Tp + k * Pp; sig = np.sqrt(C[0, 0] + k ** 2 * C[1, 1] + 2 * k * C[0, 1])
    pred.append(dict(epoch=int(k), btjd=round(tc, 4), utc=Time(tc + 2457000, format="jd", scale="tdb").utc.iso[:16], sigma_min=round(sig * 1440, 1)))
pd.DataFrame(pred).to_csv("out/214311906/predicted_dips_2027.csv", index=False)
print("2027 predicted dips:", len(pred), "first:", pred[0], "1-sigma timing", pred[0]["sigma_min"], "-", pred[-1]["sigma_min"], "min")
json.dump(out, open("out/214311906/ephemeris.json", "w"), indent=1, default=float)
