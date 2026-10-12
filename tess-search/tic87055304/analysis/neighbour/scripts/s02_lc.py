"""Light curves for the neighbour-host analysis.

Inputs (read-only):
  * PSF photometry of psf_photometry/ (lc/<variant>/psf_s*.npz): per-cadence PRF-fit fluxes of the target (jt) and
    the neighbour (jn) with all other stars modelled.  Cadence masks exactly as psf_photometry/scripts/analyze_lc.py
    (background < median + 6 MAD, chi2 < 5 x median), plus the S104 scattered-light window BTJD 4178.80-4179.20.
  * QLP SAP flux (data/qlp_raw_columns.npy) for cross-checks.
Products (data/lc_<name>.npy, dict sector -> arrays):
  pair   : F_target + F_neighbour (PSF), normalised.  Insensitive to the target/neighbour flux split (which is
           degenerate with sub-0.05 px registration errors); contains all of the neighbour's light and the target's.
  nbcat  : neighbour-only PSF flux, normalised to the CATALOGUE-predicted neighbour flux
           F_n,cat = median(F_pair) * r / (1 + r), r = F_n/F_t in the TESS band (from s03; default 0.1897).
  qlp    : QLP SAP_FLUX (target file), same cuts as transit/prep.py (QUALITY == 0, +-1 h around momentum dumps removed).
Each: t, f (detrended), f_raw (normalised, not detrended), trend, sigma (OOT MAD), exptime, bkg, seg, and for the PSF
LCs the registration offsets dx, dy.
Detrending: wotan biweight, window 1.0 d, transit mask of full width 0.30 d on P=4.136195, T0=1658.41325, segments
split at gaps > 0.5 d and at momentum dumps (QLP QUALITY bit 32 times).  5-MAD clip of OOT outliers only.
Usage: python s02_lc.py [variant=nominal] [r=0.1897]
"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from wotan import flatten, transit_mask

ROOT = "/home/user/work/tic87/neighbour"
PSF = "/home/user/work/tic87/psf_photometry/lc"
DATA = "/home/user/work/tic87/data"
P0, T00 = 4.136195, 1658.41325
SECTORS = [13, 39, 66, 93, 104]
EXPT = {13: 1800.0, 39: 600.0, 66: 200.0, 93: 200.0, 104: 200.0}
SCAT = (4178.80, 4179.20)
GAP = 0.5
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "nominal"
R_NB = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1897
raw = np.load(f"{DATA}/qlp_raw_columns.npy", allow_pickle=True).item()


def dumps_of(s):
    q = np.asarray(raw[s]["QUALITY"], np.int64); t = np.asarray(raw[s]["TIME"], float)
    return np.unique(np.round(t[(q & 32) > 0], 3))


def segments(t, breaks):
    lab = np.zeros(len(t), int); k = 0
    for i in range(1, len(t)):
        if t[i] - t[i - 1] > GAP or np.any((breaks > t[i - 1]) & (breaks <= t[i])):
            k += 1
        lab[i] = k
    return lab


def detrend(t, f, breaks, win=1.0, maskw=0.30):
    lab = segments(t, breaks)
    mask = transit_mask(time=t, period=P0, duration=maskw, T0=T00)
    trend = np.full_like(f, np.nan)
    for k in np.unique(lab):
        m = lab == k
        if m.sum() < 5: continue
        ff = f[m].copy(); ff[mask[m]] = np.nan
        if np.isfinite(ff).sum() < 5: continue
        tr = flatten(t[m], ff, method="biweight", window_length=win, return_trend=True, edge_cutoff=0)[1]
        g = np.isfinite(tr)
        if g.sum() >= 2: tr = np.interp(t[m], t[m][g], tr[g])
        trend[m] = tr
    fd = f / trend
    good = np.isfinite(trend)
    oot = good & ~mask
    mad = 1.4826 * np.nanmedian(np.abs(fd[oot] - np.nanmedian(fd[oot])))
    clip = oot & (np.abs(fd - 1) > 5 * mad)
    keep = good & ~clip
    return fd, trend, mad, keep, lab


def psf_masks(d):
    bk = d["bkg"][:, 0]; med = np.median(bk); mad = 1.4826 * np.median(np.abs(bk - med))
    c2 = d["chi2"]
    return (bk < med + 6 * mad) & (c2 < 5 * np.median(c2))


out = {"pair": {}, "nbcat": {}, "tgt": {}, "qlp": {}}
info = {}
for s in SECTORS:
    d = np.load(f"{PSF}/{VARIANT}/psf_s{s}.npz")
    t = d["t"]; jt, jn = int(d["jt"]), int(d["jn"])
    Ft, Fn = d["flux"][:, jt], d["flux"][:, jn]
    m = psf_masks(d) & np.isfinite(Ft) & np.isfinite(Fn) & ~((t > SCAT[0]) & (t < SCAT[1]))
    o = np.argsort(t[m])
    t, Ft, Fn = t[m][o], Ft[m][o], Fn[m][o]
    bkg = d["bkg"][m][o][:, 0]; reg = d["d"][m][o]
    Fp = Ft + Fn
    medp = np.median(Fp)
    Fn_cat = medp * R_NB / (1 + R_NB)
    Ft_cat = medp / (1 + R_NB)
    series = {"pair": Fp / medp,
              "nbcat": 1 + (Fn - np.median(Fn)) / Fn_cat,
              "tgt": 1 + (Ft - np.median(Ft)) / Ft_cat}
    br = dumps_of(s)
    for name, f in series.items():
        fd, tr, mad, keep, lab = detrend(t, f, br)
        out[name][s] = dict(t=t[keep], f=fd[keep], f_raw=f[keep], trend=tr[keep], sigma=mad, exptime=EXPT[s],
                            bkg=bkg[keep], dx=reg[keep, 0], dy=reg[keep, 1], seg=lab[keep], dumps=br)
    info[s] = dict(n=int(len(t)), F_t_med=float(np.median(Ft)), F_n_med=float(np.median(Fn)), F_pair_med=float(medp),
                   F_n_cat=float(Fn_cat), psf_ratio=float(np.median(Fn) / np.median(Ft)),
                   sigma_pair_ppm=float(out["pair"][s]["sigma"] * 1e6), sigma_nbcat_ppm=float(out["nbcat"][s]["sigma"] * 1e6))
    # QLP SAP (as transit/prep.py)
    v = raw[s]
    tq = np.asarray(v["TIME"], float); fq = np.asarray(v["SAP_FLUX"], float); q = np.asarray(v["QUALITY"], np.int64)
    bq = np.asarray(v["SAP_BKG"], float)
    ok = (q == 0) & np.isfinite(tq) & np.isfinite(fq) & (fq > 0) & ~((tq > SCAT[0]) & (tq < SCAT[1]))
    if len(br):
        ok &= ~(np.min(np.abs(tq[:, None] - br[None, :]), axis=1) < 1 / 24)
    tq, fq, bq = tq[ok], fq[ok], bq[ok]
    o = np.argsort(tq); tq, fq, bq = tq[o], fq[o] / np.nanmedian(fq[o]), bq[o]
    fd, tr, mad, keep, lab = detrend(tq, fq, br)
    out["qlp"][s] = dict(t=tq[keep], f=fd[keep], f_raw=fq[keep], trend=tr[keep], sigma=mad, exptime=EXPT[s], bkg=bq[keep],
                         seg=lab[keep], dumps=br)
    info[s]["sigma_qlp_ppm"] = float(mad * 1e6)
    print(s, json.dumps({k: round(v, 3) for k, v in info[s].items()}), flush=True)

sfx = "" if VARIANT == "nominal" else f"_{VARIANT}"
for name, dd in out.items():
    if name == "qlp" and sfx: continue
    np.save(os.path.join(ROOT, f"data/lc_{name}{sfx}.npy"), dd, allow_pickle=True)
json.dump(dict(variant=VARIANT, r_nb=R_NB, sectors=info), open(os.path.join(ROOT, f"data/lc_info{sfx}.json"), "w"), indent=1)
