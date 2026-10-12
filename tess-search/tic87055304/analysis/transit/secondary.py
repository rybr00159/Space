"""Secondary-eclipse search (complementary; full-phase detrended SAP flux, S104 scattered-light event masked).
Occultation model: symmetric trapezoid with the T14 and ingress fraction of the best-fit transit trapezoid
(vshape.json), exposure-integrated; depth free; per-sector offsets; fit at phase 0.5 and scanned over phase
0.30-0.70 (circular orbit expects 0.5) in 0.002 steps. Errors: per-sector white noise (MAD of out-of-transit
flux) x beta (vshape.json). Look-elsewhere: the rms of the scan in units of its errors.
Also the expected secondary depth if the signal is a background EB on the 9" neighbour (an M-dwarf companion,
Teff 3300-3600 K, behind a 5570 K star, k from mcmc_nbhost), diluted to the target light curve.
Writes secondary.json and secondary.png."""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
lc = np.load(os.path.join(HERE, "lc_prepared.npy"), allow_pickle=True).item()
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
vs = json.load(open(os.path.join(HERE, "vshape.json")))
T0, P = mc["params"]["T0"]["med"], mc["params"]["P"]["med"]
T14 = vs["models"]["T"]["T14_h"] / 24; fin = vs["models"]["T"]["f_ingress"]
BETA = {int(k): v for k, v in vs["beta"].items()}
EXPT = {13: 1800.0, 39: 600.0, 66: 200.0, 93: 200.0, 104: 200.0}


def shape(dt, E):
    """unit-depth trapezoid, box-averaged over exposure E (s)."""
    sub = np.linspace(-0.5, 0.5, 15) * E / 86400
    half = T14 / 2; T12 = fin * half
    a = np.abs(dt[:, None] + sub[None, :])
    return np.clip((half - a) / T12, 0, 1).mean(axis=1)


def fit_phase(ph0):
    rows = []
    num = 0.0; den = 0.0
    for s, v in lc.items():
        t, f = v["t"], v["f"]
        tr = np.abs(((t - T0) / P + 0.5) % 1 - 0.5) * P < 0.20   # mask primary transits
        tt, ff = t[~tr], f[~tr]
        dt = ((tt - T0) / P - ph0 + 0.5) % 1 - 0.5
        dt *= P
        m = np.abs(dt) < 0.6
        if m.sum() < 20: continue
        sx = shape(dt[m], EXPT[s])
        if sx.sum() == 0: continue
        sig = v["sigma"] * BETA[s]
        X = np.vstack([np.ones(m.sum()), -sx]).T
        c, *_ = np.linalg.lstsq(X, ff[m], rcond=None)
        cov = np.linalg.inv(X.T @ X) * sig ** 2
        d, e = c[1] / c[0], np.sqrt(cov[1, 1])
        num += d / e ** 2; den += 1 / e ** 2
        rows.append((s, d, e))
    return num / den, 1 / np.sqrt(den), rows


d05, e05, rows = fit_phase(0.5)
phs = np.arange(0.30, 0.7001, 0.002)
scan = np.array([fit_phase(p)[:2] for p in phs])
z = scan[:, 0] / scan[:, 1]
out = dict(depth_ppm_phase05=d05 * 1e6, err_ppm_phase05=e05 * 1e6, per_sector={int(s): [d * 1e6, e * 1e6] for s, d, e in rows},
           scan_max_sigma=float(z.max()), scan_max_phase=float(phs[np.argmax(z)]), scan_min_sigma=float(z.min()),
           scan_rms_sigma=float(np.std(z)), T14_h=T14 * 24, ingress_fraction=fin)
out["err_ppm_phase05_scanscaled"] = out["err_ppm_phase05"] * max(1.0, out["scan_rms_sigma"])
out["sigma_phase05_scanscaled"] = out["depth_ppm_phase05"] / out["err_ppm_phase05_scanscaled"]
# expected secondary of a background EB on the neighbour
nb = json.load(open(os.path.join(HERE, "mcmc_nbhost.json")))
k = nb["params"]["k"]["med"]
from astropy import constants as c_
import astropy.units as u
def planck_tess(T):
    lam = np.linspace(600, 1000, 200) * 1e-9
    B = 2 * c_.h.value * c_.c.value ** 2 / lam ** 5 / (np.exp(c_.h.value * c_.c.value / (lam * c_.k_B.value * T)) - 1)
    return np.trapezoid(B * lam, lam)   # photon-weighted
exp = {}
for Tc in (3100, 3300, 3600):
    sb = planck_tess(Tc) / planck_tess(5570)
    occ = k ** 2 * sb   # occultation depth on the neighbour (total occultation)
    exp[Tc] = dict(on_neighbour_ppm=occ * 1e6, in_target_lc_ppm=occ / 5.35 * 1e6)
out["expected_secondary_background_EB"] = exp
json.dump(out, open(os.path.join(HERE, "secondary.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.errorbar(phs, scan[:, 0] * 1e6, scan[:, 1] * 1e6, fmt=".", ms=2, lw=0.5)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0.5, color="C3", lw=0.6)
ax.set_xlabel("phase"); ax.set_ylabel("occultation depth (ppm)")
ax.set_title(f"phase 0.5: {d05 * 1e6:.0f} +- {e05 * 1e6:.0f} ppm; scan max {z.max():.1f} sigma at {phs[np.argmax(z)]:.3f}", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "secondary.png"), dpi=120)
