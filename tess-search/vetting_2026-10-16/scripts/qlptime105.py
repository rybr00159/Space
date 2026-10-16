"""Is QLP S105 (and S106) also missing the barycentric correction? Time TIC 71316387's EB eclipses in QLP SAP and in
TESScut aperture photometry of the same sector."""
import sys, glob, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from wotan import flatten
import stage2
P, T0, T14 = 1.30367228, 3165.670416, 0.0667
def trapv(t, tc, d, T):  # V-ish eclipse
    a = np.abs(t - tc); return 1 - d * np.clip(1 - a / (T / 2), 0, 1)
def times(t, f):
    ph = np.abs(((t - T0) / P + 0.5) % 1 - 0.5) * P; m = (ph > 0.1) & (np.abs(np.abs(((t - T0) / P) % 1 - 0.5) * P) > 0.1)
    tr = flatten(t[m], f[m], method="biweight", window_length=0.8, return_trend=True, break_tolerance=0.4)[1]
    k = np.isfinite(tr); f = f / np.interp(t, t[m][k], tr[k]); s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
    out = {}
    for e in np.unique(np.round((t - T0) / P)):
        tc = T0 + e * P; w = np.abs(t - tc) < 0.15
        if np.sum(np.abs(t[w] - tc) < T14 / 2) < 5: continue
        r = least_squares(lambda x: (f[w] - trapv(t[w], tc + x[0], x[1], T14)) / s, [0, 0.006], bounds=([-0.02, 0], [0.02, 0.05]))
        J = r.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * r.cost / (w.sum() - 2))
        if r.x[1] > 0.002: out[int(e)] = (tc + r.x[0], np.sqrt(cov[0, 0]))
    return out
p = stage2.products(71316387)
for sec in (105, 106):
    import io, requests, fetch_qlp
    resp = requests.get(fetch_qlp.url(71316387, sec), timeout=120)
    r = [1] if resp.status_code == 200 else []
    cut = glob.glob(f"/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut/s{sec}/*astrocut.fits")
    if not len(r) or not cut: print(sec, "missing", len(r), cut); continue
    with fits.open(io.BytesIO(resp.content)) as f:
        d = f[1].data; t = np.array(d["TIME"], float); fl = np.array(d["SAP_FLUX"], float); q = np.array(d["QUALITY"])
    ok = (q == 0) & np.isfinite(t) & np.isfinite(fl); A = times(t[ok], fl[ok] / np.nanmedian(fl[ok]))
    with fits.open(cut[0]) as f:
        d = f[1].data; tc = np.array(d["TIME"], float); cube = np.array(d["FLUX"], float); qq = np.array(d["QUALITY"]); tcorr = np.array(d["TIMECORR"], float)
    ny, nx = cube.shape[1:]; yy, xx = np.mgrid[0:ny, 0:nx]; ap = np.hypot(xx - nx // 2, yy - ny // 2) <= 1.5
    good = (qq == 0) & np.isfinite(tc) & np.all(np.isfinite(cube.reshape(len(tc), -1)), 1)
    med = np.median(cube[good], 0); bk = np.median(cube[good][:, med < np.percentile(med, 30)], 1)
    fc = cube[good][:, ap].sum(1) - bk * ap.sum(); B = times(tc[good], fc / np.median(fc))
    common = sorted(set(A) & set(B))
    dd = np.array([(A[e][0] - B[e][0]) * 1440 for e in common]); ss = np.array([np.hypot(A[e][1], B[e][1]) * 1440 for e in common])
    w = 1 / ss ** 2
    print(f"S{sec}: {len(common)} eclipses; QLP minus TESScut = {np.sum(w*dd)/w.sum():+.2f} ± {1/np.sqrt(w.sum()):.2f} min; TESScut TIMECORR median {np.nanmedian(tcorr)*1440:+.2f} min")
