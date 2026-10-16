"""Ephemeris of the 5.69 d signal from the QLP SAP light curve of TIC 214311908 (strongest dips), with QLP S103/S104
times corrected by the missing barycentric delay (TESScut TIMECORR median for this field), then 2027 predictions."""
import sys, json, glob, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "214311906"]
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import least_squares
from wotan import flatten
import stage2
from lcan import read, bin10, in_mask, trap
P0, T00, T14, FL = 5.69149, 3840.913, 0.067, 0.78
shift = {}
for s in (103, 104):
    f = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{s:04d}-*_astrocut.fits")[0]
    shift[s] = float(np.nanmedian(fits.open(f)[1].data["TIMECORR"]))
print("QLP time corrections (d):", {k: round(v, 5) for k, v in shift.items()}, "=", {k: round(v * 1440, 2) for k, v in shift.items()}, "min")
ev = []
p = stage2.products(214311908)
for _, r in p[p.prov == "QLP"].sort_values("sector").iterrows():
    sec = int(r.sector)
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), "QLP")
    t = t + shift.get(sec, 0.0)
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols["SAP_FLUX"]) & (cols["SAP_FLUX"] > 0)
    tt, ff = bin10(t[ok], cols["SAP_FLUX"][ok] / np.nanmedian(cols["SAP_FLUX"][ok]), r.exptime)
    m = ~in_mask(tt, P0, T00, T14, 2.0)
    tr = flatten(tt[m], ff[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    k_ = np.isfinite(tr); f = ff / np.interp(tt, tt[m][k_], tr[k_]); s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
    for e in np.unique(np.round((tt - T00) / P0)):
        tc = T00 + e * P0; w = np.abs(tt - tc) < 0.3
        if np.sum(np.abs(tt[w] - tc) < T14) < 4 or w.sum() < 15: continue
        grid = np.arange(-0.12, 0.1201, 0.002); c2 = []
        for dt in grid:
            y = 1 - trap(tt[w], tc + dt, 1.0, T14, T14 * FL); A = np.column_stack([np.ones(w.sum()), y])
            c, *_ = np.linalg.lstsq(A, f[w], rcond=None); c2.append(np.sum(((f[w] - A @ c) / s) ** 2) if c[1] < 0 else np.inf)
        j = int(np.argmin(c2))
        if not np.isfinite(c2[j]): continue
        fit = least_squares(lambda x: (f[w] - x[1] + x[2] * (1 - trap(tt[w], tc + x[0], 1.0, T14, T14 * FL))) / s, [grid[j], 1.0, 0.01], bounds=([-0.13, 0.9, 0], [0.13, 1.1, 0.2]))
        J = fit.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * fit.cost / (w.sum() - 3))
        snr = fit.x[2] / np.sqrt(cov[2, 2])
        if snr > 3.5: ev.append(dict(sector=sec, k=int(e), tc=float(tc + fit.x[0]), err=float(np.sqrt(cov[0, 0])), depth=float(fit.x[2]), snr=float(snr)))
ev = pd.DataFrame(ev)
def linfit(d):
    W = 1 / d.err.values; A = np.column_stack([np.ones(len(d)), d.k.values]); c, *_ = np.linalg.lstsq(A * W[:, None], d.tc.values * W, rcond=None)
    cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None])); res = d.tc.values - A @ c; chi2 = float(np.sum((res * W) ** 2))
    return c, cov, res, chi2
c, cov, res, chi2 = linfit(ev)
# iterative 4-sigma clip of outlying events (e.g. events hit by systematics), reported
keep = np.abs(res / ev.err.values) < 4 * np.sqrt(max(1, chi2 / (len(ev) - 2)))
c, cov, res2, chi2 = linfit(ev[keep]); dof = int(keep.sum()) - 2; cov = cov * max(1, chi2 / dof)
print(f"{len(ev)} events (S/N>3.5), {int((~keep).sum())} clipped: P = {c[1]:.6f} ± {np.sqrt(cov[1,1]):.6f} d, T0 = {c[0]:.4f} ± {np.sqrt(cov[0,0]):.4f} BTJD; chi2 {chi2:.1f}/{dof}; O-C rms {np.std(res2)*1440:.1f} min")
print("events per sector:", ev[keep].groupby("sector").size().to_dict(), "; median depth by sector (%):", (ev[keep].groupby("sector").depth.median() * 100).round(2).to_dict())
Tp, Pp = c[0], c[1]
start = Time("2027-02-15").jd - 2457000; stop = Time("2027-10-15").jd - 2457000
ks = np.arange(np.ceil((start - Tp) / Pp), np.floor((stop - Tp) / Pp) + 1); pred = []
for k in ks:
    tc = Tp + k * Pp; sig = np.sqrt(cov[0, 0] + k ** 2 * cov[1, 1] + 2 * k * cov[0, 1])
    pred.append(dict(epoch=int(k), bjd_tdb=round(tc + 2457000, 4), utc_mid=Time(tc + 2457000, format="jd", scale="tdb").utc.iso[:16], sigma_min=round(sig * 1440, 1)))
pd.DataFrame(pred).to_csv("out/214311906/predicted_transits_2027.csv", index=False)
print(f"{len(pred)} predicted transits Feb-Oct 2027; timing 1-sigma {pred[0]['sigma_min']}-{pred[-1]['sigma_min']} min; first {pred[0]}")
json.dump(dict(P=Pp, P_err=float(np.sqrt(cov[1, 1])), T0=Tp, T0_err=float(np.sqrt(cov[0, 0])), cov=cov.tolist(), chi2=chi2, dof=dof,
               qlp_shift_days=shift, events=ev.to_dict("records"), clipped=ev[~keep].to_dict("records")),
          open("out/214311906/ephemeris_v2.json", "w"), indent=1, default=float)
