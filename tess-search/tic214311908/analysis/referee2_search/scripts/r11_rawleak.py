"""r11: is the RRab signal present in the RAW pixels around the M3V (then any photometry would carry it), or only in
the phot model-B light curve (then it is a pipeline coupling, as the analyst proposes)?

Per sector, from the TESScut cube (QUALITY == 0): flux sums over (a) the 3x3 box centred on the M3V (contains the Mira
core, 0.85 px away), (b) the single pixel nearest the M3V, (c) a control 3x3 box 4 px from the M3V on the side away
from the Mira. Each series: 5.69-d windows removed, divided by its 1.5-d running median (relative), 3-harmonic fit at
the Gaia RRab period with phases referred to the Gaia epoch of maximum. Reported: first-harmonic semi-amplitude in
e-/s (relative amplitude x median flux) and in units of the M3V-B leak amplitude; the projection onto the RRab phase
(negative = inverted like the M3V-B leak); significance from 100 fake periods 0.5-0.9 d.
Variant: with and without subtracting a common mode (median of pixels > 6 px from the Mira).
Output: out/r11_rawleak.json
"""
import sys, glob
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter

CUBES = "/home/user/work/tesscut5/214311906"
M3V = (260.969674, -51.323821); MIRA = (260.9776003, -51.3237201)
rng = np.random.default_rng(1111)


def h1(t, y, Pp, Tmax):
    ph = 2 * np.pi * (t - Tmax) / Pp
    X = np.column_stack([np.ones(len(t))] + [f(h * ph) for h in (1, 2, 3) for f in (np.cos, np.sin)])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return c[1], c[2]          # cos (projection on 'maximum at Tmax'), sin


def main():
    g = pd.read_csv(f"{ME}/data/gaia_rrab.csv").iloc[0]
    Pp = float(g.pf); Tmax = float(g.epoch_g) + 2455197.5 - 2457000.0
    leak = json.load(open(f"{ME}/out/r01_rrab_leak.json"))["C_per_sector"]
    out = {}
    for sec in SECTORS:
        fn = glob.glob(f"{CUBES}/tess-s{sec:04d}-*_astrocut.fits")[0]
        with fits.open(fn, memmap=True) as h:
            tab = h[1].data; w = WCS(h[2].header)
            t = np.array(tab["TIME"], float); F = np.array(tab["FLUX"], np.float64); q = np.array(tab["QUALITY"], int)
        ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1) & (np.abs(phase_dist(t)) > 0.1)
        t, F = t[ok], F[ok]
        mx, my = w.all_world2pix([M3V], 0)[0]; Mx, My = w.all_world2pix([MIRA], 0)[0]
        ix, iy = int(round(mx)), int(round(my))
        yy, xx = np.mgrid[0:21, 0:21]
        cm = np.median(F[:, np.hypot(xx - Mx, yy - My) > 6], axis=1)
        sx = int(np.sign(mx - Mx)) or 1
        series = {"box3x3": F[:, iy - 1:iy + 2, ix - 1:ix + 2].sum((1, 2)), "pix": F[:, iy, ix],
                  "control3x3": F[:, iy - 1:iy + 2, ix + 4 * sx - 1:ix + 4 * sx + 2].sum((1, 2))}
        res = {}
        dt = np.median(np.diff(t))
        win = max(5, int(1.5 / dt) | 1)
        for nm, s0 in series.items():
            for cmv in (False, True):
                s = s0 - (9 if "3x3" in nm else 1) * cm if cmv else s0
                tr = median_filter(s, size=win, mode="nearest")
                y = s / tr - 1
                y = np.where(np.abs(y - np.median(y)) < 6 * rstd(y), y, np.nan); k = np.isfinite(y)
                a, b = h1(t[k], y[k], Pp, Tmax)
                amp = np.hypot(a, b) * np.median(tr)
                proj = a * np.median(tr)            # e-/s along the RRab max phase (+ = in phase with RRab light)
                fk = []
                for pf in rng.uniform(0.5, 0.9, 100):
                    if abs(pf - Pp) < 0.01 or abs(pf - Pp / 2) < 0.01: continue
                    aa, bb = h1(t[k], y[k], pf, Tmax); fk.append(np.hypot(aa, bb) * np.median(tr))
                fk = np.array(fk)
                res[f"{nm}{'_cm' if cmv else ''}"] = dict(amp1_e_s=float(amp), proj_on_RR_phase_e_s=float(proj), fake_median=float(np.median(fk)),
                                                          fake_p95=float(np.percentile(fk, 95)), p_emp=float(np.mean(fk >= amp)),
                                                          median_flux=float(np.median(tr)))
        res["M3V_B_leak_amp1_e_s"] = leak[str(sec)]["amp_ppt"][0] * FREF / 1e3
        out[sec] = res
        print(sec, json.dumps({k: ({kk: round(vv, 2) for kk, vv in v.items()} if isinstance(v, dict) else round(v, 2)) for k, v in res.items()}), flush=True)
        del F
    jdump(out, f"{ME}/out/r11_rawleak.json")


if __name__ == "__main__":
    main()
