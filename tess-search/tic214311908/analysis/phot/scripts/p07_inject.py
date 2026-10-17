"""Injection-recovery in the pixel data: cross-talk of a transit on one star into the deblended light curves of the
others, for each sector and each photometric model (B, Bp, A).

For each sector and host X in {M3V, pair, Mira}: a train of trapezoid transits (T14 = 0.067 d, flat fraction 0.78,
exposure-integrated) removing DELTA e-/s at mid-transit from star X is subtracted from the cutout cube using X's
image: the SPOC PRF at X's position for the M dwarfs (plus a 'mismatched' M3V variant: PRF blurred by 0.15 px and
displaced by (0.10, -0.10) px, i.e. a template error the fit does not know about), and the Mira's empirical image
(p05 b image, normalised to unit sum within 3.5 px) for the Mira. Fake ephemerides: P = 5.6915 d at phase offsets
0.2, 0.35, 0.5, 0.65, 0.8 (real transits avoided), plus P = 3.71 d and 8.27 d. The full pipeline (regression with the
fake windows masked, per-cadence fit, detrending) is rerun and the dip is measured on every light curve with the same
linear template fit as for the real signal. Response = recovered(injected) - recovered(same phases, no injection).
Outputs: out/inject_s{sector}.json and out/inject_summary.json (cross-talk matrices per model)
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p05_phot import setup, photometry
from p06_lc import detrend, FREF, BKG_EXC_MAX

DELTA = 15.0      # e/s removed at mid-transit (1.6% of the M3V's 940 e/s)
T14, FL = 0.067, 0.78
KEYS = (("B", "M3V"), ("B", "Mira"), ("Bp", "pair"), ("A", "M3V"), ("A", "pair"), ("A", "Mira"))


def shape(t, P, T0, exptime, nsub=7):
    """Exposure-averaged trapezoid (1 at full depth)."""
    off = (np.arange(nsub) + 0.5) / nsub - 0.5
    y = np.zeros(len(t))
    for o in off:
        tt = t + o * exptime
        ph = ((tt - T0) / P + 0.5) % 1 - 0.5; a = np.abs(ph * P); T23 = FL * T14
        y += np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
    return y / nsub


def windows(t, P, T0, half=0.1):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    return np.abs(ph * P) < half


def measure(R, P, T0, exptime, fref):
    """Dip (e/s at full depth) on each light curve at ephemeris (P, T0), with real + fake windows masked in detrending."""
    t = R["t"]
    chiB = R["FIT"]["B"]["chi2"]
    bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
    ok = R["qok"] & (bexc < BKG_EXC_MAX) & R["jit_ok"] & np.isfinite(chiB) & (chiB < 3 * np.nanmedian(chiB[R["q0"]]))
    mask = windows(t, P_EPH, T0_EPH) | windows(t, P, T0)
    y = shape(t, P, T0, exptime)
    out = {}
    for mk, nm in KEYS:
        names = R["FIT"][mk]["names"]; dF = R["FIT"][mk]["F"][:, names.index(nm)]
        raw = 1 + dF / fref[nm]
        tr = detrend(t, raw, ok, mask)
        d = raw / tr
        m = ok & np.isfinite(d) & ~windows(t, P_EPH, T0_EPH)        # never use real transits
        if np.sum(m & (y > 0.5)) < 3: out[f"{mk}:{nm}"] = (np.nan, np.nan); continue
        A = np.column_stack([np.ones(m.sum()), -y[m]]); c, *_ = np.linalg.lstsq(A, d[m], rcond=None)
        sd = robust_std(d[m & (y == 0)]); e = sd / np.sqrt(np.sum((y[m] - y[m].mean()) ** 2))
        out[f"{mk}:{nm}"] = (float(c[1] * fref[nm]), float(e * fref[nm]))
    return out


def run(sector):
    S = setup(sector); C, W, J, g, prf, rBA = S
    R0 = photometry(sector, S=S, verbose=False)
    fref = dict(FREF, Mira=float(R0["mu"]["F"] - 0.9 * (FREF["M3V"] + FREF["pair"])))
    t = C.t; exptime = C.exptime
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mx, my = R0["pos"]["Mira"]
    b = R0["b_img"] / R0["b_img"][np.hypot(xx - mx, yy - my) <= 3.5].sum()
    b = np.where(np.hypot(xx - mx, yy - my) <= 3.5, b, 0.0)
    imgs = dict(M3V=prf.image(*R0["pos"]["M3V"]), pair=prf.image(*R0["pos"]["pair"]), Mira=b,
                M3V_mismatch=prf.image(R0["pos"]["M3V"][0] + 0.10, R0["pos"]["M3V"][1] - 0.10, blur=np.hypot(prf.blur, 0.15)))
    ephs = [(P_EPH, T0_EPH + f * P_EPH, f"P5.69_ph{f}") for f in (0.2, 0.35, 0.5, 0.65, 0.8)] + [(3.71, T0_EPH + 0.3, "P3.71"), (8.27, T0_EPH + 1.7, "P8.27")]
    res = []
    for P, T0, lab in ephs:
        tt = t[np.isfinite(t)]
        fake_epochs = T0 + P * np.arange(np.floor((tt.min() - T0) / P) - 1, np.ceil((tt.max() - T0) / P) + 2)
        base = measure(R0, P, T0, exptime, fref)
        for host, img in imgs.items():
            amp = DELTA * shape(np.nan_to_num(t), P, T0, exptime)
            R = photometry(sector, S=S, inject=[dict(img=img, amp=amp)], extra_mask_times=fake_epochs, verbose=False)
            rec = measure(R, P, T0, exptime, fref)
            row = dict(sector=sector, eph=lab, host=host, P=P, T0=T0)
            for k in rec:
                row[f"rec_{k}"] = rec[k][0] - base[k][0]; row[f"err_{k}"] = rec[k][1]; row[f"base_{k}"] = base[k][0]
            res.append(row)
        for r in res[-len(imgs):]:
            print(f"S{sector} {lab} inj {r['host']:13s}: " + ", ".join(f"{mk}:{nm} {r['rec_'+mk+':'+nm]:6.1f}" for mk, nm in KEYS), flush=True)
    save_json(dict(sector=sector, delta=DELTA, fref=fref, rows=res), f"{ROOT}/out/inject_s{sector}.json")
    return res


if __name__ == "__main__":
    allr = []
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        allr += run(s)
