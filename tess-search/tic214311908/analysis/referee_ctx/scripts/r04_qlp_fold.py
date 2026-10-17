"""R04: independent QLP depth / odd-even / phase-0.5 measurement for the M3V, pair A and the Mira.

Estimator: joint fit, common trapezoid depth, quadratic baseline per +-0.3 d window (rlib.joint_fit), on QLP
SAP_FLUX and DET_FLUX (QUALITY == 0), S39-S104 (S12 30-min reported separately). S103/S104 times corrected with the
qlp_shift_days of ephemeris_v2.json (QLP lacks the barycentric correction there). Errors: the same estimator on fake
window sets (all real windows shifted by a common offset; offsets >= 0.6 d apart and away from phase 0 and 0.5 are
independent) -> empirical sd. Odd/even: separate depths for even/odd epochs; difference with empirical error from the
same fake sets (odd/even split of the fakes).
Output: out/r04_qlp_fold.json
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
import rlib as R

LC = "tic214311908/analysis/ctx/lc"
TIC = {"M3V": "0000000214311908", "pairA": "0000000214311910", "Mira": "0000000214311906"}
SECT = (39, 66, 93, 103, 104)


def load(tic, sec, col):
    f = glob.glob(f"{LC}/hlsp_qlp_tess_ffi_s{sec:04d}-{tic}_tess_v01_llc.fits")
    if not f: return None
    d = fits.getdata(f[0], 1)
    col = col if col in d.columns.names else "KSPSAP_FLUX"
    t = np.array(d["TIME"], float) + R.QLP_SHIFT.get(sec, 0.0); fl = np.array(d[col], float)
    ok = (d["QUALITY"] == 0) & np.isfinite(t) & np.isfinite(fl)
    return t[ok], fl[ok] / np.nanmedian(fl[ok])


def measure(data, offset=0.0, phase=0.0, parity=None):
    """Joint depth over all sectors; returns depth (fraction)."""
    rows = []
    T = np.concatenate([d[0] for d in data]); Fl = np.concatenate([d[1] for d in data]); EX = np.concatenate([np.full(len(d[0]), d[2]) for d in data])
    tcs, ks = R.centres(T.min(), T.max(), phase)
    tcs = tcs + offset
    keep = np.array([np.any(np.abs(T - c) < 0.05) for c in tcs])
    tcs, ks = tcs[keep], ks[keep]
    if parity is not None:
        sel = (ks % 2) == parity; tcs = tcs[sel]
    # joint fit across sectors: exposure differs by sector; use window-specific exposure
    blocks = []
    for c in tcs:
        m = np.abs(T - c) < 0.30
        if m.sum() < 15: continue
        s = R.shape(T[m], c, EX[m][0])
        if np.sum(s > 0.5) < 2 or np.sum(s == 0) < 10: continue
        blocks.append((m, s, c))
    if not blocks: return np.nan, 0
    nb = len(blocks); npar = 1 + 3 * nb; A_ = []; y = []
    for b, (m, s, c) in enumerate(blocks):
        x = (T[m] - c) / 0.3; A = np.zeros((m.sum(), npar)); A[:, 0] = -s; A[:, 1 + 3 * b] = 1; A[:, 2 + 3 * b] = x; A[:, 3 + 3 * b] = x * x
        A_.append(A); y.append(Fl[m])
    coef = np.linalg.lstsq(np.vstack(A_), np.concatenate(y), rcond=None)[0]
    return coef[0], nb


out = {}
for star, tic in TIC.items():
    for col in ("SAP_FLUX", "DET_FLUX"):
        data = []
        for sec in SECT:
            L = load(tic, sec, col)
            if L is not None: data.append((L[0], L[1], R.exptime(sec)))
        d0, nb = measure(data); d5, _ = measure(data, phase=0.5)
        de, _ = measure(data, parity=0); do, _ = measure(data, parity=1)
        # fake offsets
        offs = np.arange(0.6, R.P - 0.6, 0.05)
        offs = offs[(np.abs(offs - R.P / 2) > 0.6)]
        fk = np.array([[measure(data, o)[0], measure(data, o, parity=0)[0], measure(data, o, parity=1)[0]] for o in offs])
        ind = [0]
        for i in range(1, len(offs)):
            if offs[i] - offs[ind[-1]] >= 0.6: ind.append(i)
        e_all = np.nanstd(fk[:, 0]); e_ind = np.nanstd(fk[ind, 0])
        e_diff = np.nanstd(fk[:, 2] - fk[:, 1])
        out[f"{star}_{col}"] = dict(depth_ppm=d0 * 1e6, n_windows=nb, err_fake_ppm=e_all * 1e6, err_fake_indep_ppm=e_ind * 1e6,
                                    n_fake=len(offs), n_fake_indep=len(ind), snr=d0 / e_all,
                                    phase05_ppm=d5 * 1e6, even_ppm=de * 1e6, odd_ppm=do * 1e6,
                                    odd_minus_even_ppm=(do - de) * 1e6, odd_minus_even_err_fake_ppm=e_diff * 1e6,
                                    odd_even_sigma=(do - de) / e_diff)
        print(star, col, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in out[f"{star}_{col}"].items()}, flush=True)
# per-sector M3V DET depths
per = {}
for sec in SECT:
    L = load(TIC["M3V"], sec, "DET_FLUX"); d, nb = measure([(L[0], L[1], R.exptime(sec))])
    per[sec] = dict(depth_ppm=d * 1e6, n=nb)
out["M3V_DET_per_sector"] = per; print(per)
json.dump(out, open("out/r04_qlp_fold.json", "w"), indent=1, default=float)
