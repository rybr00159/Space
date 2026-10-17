"""s06: Lomb-Scargle periodograms of the public light curves (ctx/lc, read-only): TARS (TIC 214311908, 1-px aperture,
S12-S93), QLP SAP and DET (M3V 214311908, pair 214311910, Mira 214311906, TIC 1515709367), ELEANOR-LITE S12.
Questions: (1) is the ~8.5-9.2 d TARS peak present and in which products/stars; (2) is the 0.6485-d signal (s01)
present, with which amplitude (relative and in e-/s of the aperture) in each product.
The 5.69-d transit windows (+-0.1 d) are removed. Errors are not used (LS on 30-min bins, median-normalised per orbit
segment for DET; SAP only median-normalised per sector).
Output: out/public_periodograms.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from astropy.io import fits
from astropy.timeseries import LombScargle
import glob

LC = f"{BASE}/ctx/lc"
P065 = 0.64853


def bin30(t, f, w=1 / 48):
    k = np.floor((t - t[0]) / w).astype(int); u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    return np.bincount(inv, t) / n, np.bincount(inv, f) / n


def ls_summary(t, f, pmin=0.2, pmax=15.0):
    fr = np.linspace(1 / pmax, 1 / pmin, 60000)
    L = LombScargle(t, f); pw = L.power(fr)
    order = np.argsort(pw)[::-1]
    peaks = []
    for i in order:
        if len(peaks) >= 5: break
        if all(abs(fr[i] - fr[j]) > 0.02 * fr[j] for j in peaks): peaks.append(i)
    near = np.abs(1 / fr - P065) < 0.006
    j = np.argmax(np.where(near, pw, 0))
    # semi-amplitude of a sinusoid fit at the local best period near 0.6485
    X = np.column_stack([np.ones(len(t)), np.sin(2 * np.pi * fr[j] * t), np.cos(2 * np.pi * fr[j] * t)])
    c, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ c
    # empirical significance: power at 200 random frequencies in 1.2-2 c/d
    fr2 = RNG.uniform(1 / 0.8, 1 / 0.52, 300); pr = L.power(fr2)
    k = (np.abs(1 / fr - 9.0) < 1.0)
    jj = np.argmax(np.where(k, pw, 0))
    return dict(top=[(float(1 / fr[i]), float(pw[i])) for i in peaks],
                p065=dict(P=float(1 / fr[j]), power=float(pw[j]), semi_amp=float(np.hypot(c[1], c[2])),
                          power_rank_vs_random_1p2_2cd=float((pr >= pw[j]).mean())),
                p8_10=dict(P=float(1 / fr[jj]), power=float(pw[jj])), rms=float(np.std(f)))


RNG = np.random.default_rng(5)


def main():
    out = {}
    # TARS
    for f in sorted(glob.glob(f"{LC}/hlsp_tars_*214311908*.fits")):
        h = fits.open(f); s = int(f.split("_s00")[1][:2]) if "_s00" in f else None
        d = h[1].data; t = np.array(d["TIME"], float); fl = np.array(d["FLUX"], float)
        m = np.isfinite(t) & np.isfinite(fl) & ~transit_mask(t)
        t, fl = t[m], fl[m] / np.nanmedian(fl[m]) - 1
        tb, fb = bin30(t, fl)
        sec = int(f.split("ffi_s")[1][:4])
        out[f"TARS_M3V_S{sec}"] = ls_summary(tb, fb)
        out[f"TARS_M3V_S{sec}"]["header_periods"] = [h[0].header.get(f"PER{i}") for i in range(1, 6)]
    # QLP
    for f in sorted(glob.glob(f"{LC}/hlsp_qlp_*.fits")):
        h = fits.open(f); d = h[1].data
        bn = os.path.basename(f); sec = int(bn.split("ffi_s")[1][:4]); tic = int(bn.split("-")[1][:16])
        t = np.array(d["TIME"], float)
        if sec in (103, 104): t = t + (5.74 if sec == 103 else 7.27) / 1440.0      # referee_phot QLP time fix (irrelevant here)
        q = np.array(d["QUALITY"], int)
        for col in [c for c in ("SAP_FLUX", "DET_FLUX", "KSPSAP_FLUX") if c in d.columns.names]:
            fl = np.array(d[col], float)
            m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & ~transit_mask(t)
            if m.sum() < 200: continue
            tt, ff = t[m], fl[m]
            med = np.nanmedian(ff)
            ff = ff / med - 1
            # SAP: remove a slow trend (biweight 2 d) for the 0.65-d question only
            tb, fb = bin30(tt, ff)
            key = f"QLP_{tic}_S{sec}_{col}"
            out[key] = ls_summary(tb, fb)
            out[key]["median_flux_e_s"] = float(med)
            if col == "SAP_FLUX":
                from wotan import flatten
                fd = flatten(tb, fb + 1, method="biweight", window_length=2.0, break_tolerance=0.5) - 1
                k = np.isfinite(fd)
                out[key + "_flat2d"] = ls_summary(tb[k], fd[k])
                out[key + "_flat2d"]["median_flux_e_s"] = float(med)
    save_json(out, f"{ROOT}/out/public_periodograms.json")
    for k, v in out.items():
        print(f"{k:42s} top {v['top'][0][0]:7.3f} ({v['top'][0][1]:.3f}) {v['top'][1][0]:7.3f} ({v['top'][1][1]:.3f}) | "
              f"0.6485: P {v['p065']['P']:.4f} pow {v['p065']['power']:.3f} A {v['p065']['semi_amp']*1e3:.2f} ppt FAPemp {v['p065']['power_rank_vs_random_1p2_2cd']:.3f}"
              + (f" A_e/s {v['p065']['semi_amp']*v['median_flux_e_s']:.1f}" if 'median_flux_e_s' in v else "") +
              f" | 8-10d: {v['p8_10']['P']:.2f} ({v['p8_10']['power']:.3f})")


if __name__ == "__main__":
    main()
