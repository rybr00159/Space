"""Two-star SED test (RUWE 18.8 suggests an unresolved companion): primary + companion at the same distance and
extinction, CK04 atmospheres, F99 dust. Same photometry as the 'noapass' single-star fit (Gaia, SkyMapper, 2MASS,
WISE W1/W2; jitter fixed to the single-star best value).
For a grid of companion Teff (3500-6000 K, logg 4.6 for T2 < 5500 K else 4.4) and companion/primary flux ratio in the
TESS band (0-1.2), the primary (Teff1, [M/H], A_V, theta1; logg1 4.3) is optimised and chi^2 recorded.
Results: best two-star chi^2 vs single-star chi^2 (Delta chi^2, BIC), and the 2-sigma (Delta chi^2 = 4 above the
minimum along the flux-ratio axis, profile likelihood) upper limit on the companion's TESS-band flux ratio per T2.
Companion masses for each T2 are taken from MIST ([Fe/H] = -0.3, 5 Gyr main sequence).
Output: out/sed_binary.json, fig/sed_binary.png
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, pandas as pd
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import sedlib
ROOT = "/home/user/work/tic87/stellar"
ph = pd.read_csv(os.path.join(ROOT, "out/sed_photometry.csv"))
ph = ph[(ph.band != "W3") & ~ph.source.isin(["APASS DR9", "Tycho-2"])]
bands = list(ph.band); mobs = ph.mag.values
single = json.load(open(os.path.join(ROOT, "out/sed/noapass.json")))
sj = single["best"]["sj"]
sig2 = ph.sigma.values ** 2 + sj ** 2
allb = list(pd.read_csv(os.path.join(ROOT, "out/sed_photometry.csv")).band) + ["TESS", "U", "R", "I"]
M = sedlib.SEDModel(allb)
idx = np.array([allb.index(b) for b in bands]); iT = allb.index("TESS")


def flux(T, lg, mh, av, th):
    return 10 ** (-0.4 * M.mags(T, lg, mh, av, th))


def chi2_two(p, T2, fr):
    T1, mh, av, lth1 = p
    if not (4500 < T1 < 7500 and -1.5 < mh < 0.5 and 0 <= av < 1.5): return 1e9
    th1 = 10 ** lth1
    f1 = flux(T1, 4.3, mh, av, th1)
    if fr > 0:
        # companion angular diameter chosen so that its TESS flux = fr * primary TESS flux
        lg2 = 4.6 if T2 < 5500 else 4.4
        f2u = flux(T2, lg2, mh, av, 1.0)
        th2 = np.sqrt(fr * f1[iT] / f2u[iT])
        f2 = f2u * th2 ** 2
    else:
        f2 = 0 * f1
    m = -2.5 * np.log10(f1 + f2)
    if not np.all(np.isfinite(m[idx])): return 1e9
    return float(np.sum((mobs - m[idx]) ** 2 / sig2))


def best(T2, fr, x0):
    r = None
    for start in [x0, [5900, -0.4, 0.25, np.log10(0.03)], [6200, -0.3, 0.4, np.log10(0.028)], [5600, -0.6, 0.1, np.log10(0.031)]]:
        rr = minimize(chi2_two, start, args=(T2, fr), method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-4, maxiter=4000))
        if r is None or rr.fun < r.fun: r = rr
    return r


x0 = [5900, -0.4, 0.25, np.log10(0.0311)]
r0 = best(0, 0.0, x0); chi_single = r0.fun
print("single star chi2", chi_single, r0.x)
T2s = np.arange(3500, 6001, 250); frs = np.array([0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2])
C = np.full((len(T2s), len(frs)), np.nan); P = {}
for i, T2 in enumerate(T2s):
    xp = r0.x
    for j, fr in enumerate(frs):
        r = best(T2, fr, xp); C[i, j] = r.fun; xp = r.x; P[(T2, fr)] = r.x.tolist()
    print(T2, np.round(C[i] - chi_single, 2), flush=True)
cmin = np.nanmin(C)
# MIST masses for companion Teff on the 5 Gyr, [Fe/H]=-0.3 main sequence
z = np.load(os.path.join(ROOT, "cache/mist_compact.npz")); t = pd.DataFrame(z["tab"], columns=list(z["cols"]))
ms = t[(t.feh == -0.25) & (np.abs(t.logage - 9.7) < 0.01) & (t.phase == 0)].sort_values("mini")
M2 = np.interp(np.log10(T2s), ms.logT.values, ms.mass.values)
ul = []
for i, T2 in enumerate(T2s):
    prof = C[i] - cmin
    ok = np.where(prof <= 4.0)[0]
    j = ok.max() if len(ok) else 0
    # interpolate the flux ratio where profile crosses 4
    if j < len(frs) - 1 and prof[j + 1] > 4:
        fr_ul = np.interp(4.0, [prof[j], prof[j + 1]], [frs[j], frs[j + 1]])
    else:
        fr_ul = frs[j]
    ul.append(float(fr_ul))
k_extra = 1  # one extra free parameter (flux ratio) for given T2; T2 itself is a second
bic_single = chi_single; bic_two = cmin + 2 * np.log(len(bands))
ib = np.unravel_index(np.nanargmin(C), C.shape)
res = dict(n_bands=len(bands), jitter_fixed=sj, chi2_single=chi_single, chi2_two_min=float(cmin),
           delta_chi2=float(chi_single - cmin), best_two=dict(T2=float(T2s[ib[0]]), flux_ratio_TESS=float(frs[ib[1]])),
           delta_BIC_two_minus_single=float(bic_two - bic_single),
           upper_limit_flux_ratio_TESS_2sigma={str(int(T)): u for T, u in zip(T2s, ul)},
           companion_mass_MIST={str(int(T)): float(m) for T, m in zip(T2s, M2)},
           chi2_grid=C.tolist(), T2=T2s.tolist(), fr=frs.tolist(), primary_single=r0.x.tolist())
json.dump(res, open(os.path.join(ROOT, "out/sed_binary.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("chi2_grid",)}, indent=1))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
im = ax[0].pcolormesh(frs, T2s, C - cmin, vmin=0, vmax=20, shading="nearest", cmap="viridis_r")
ax[0].contour(frs, T2s, C - cmin, levels=[2.3, 4.0, 9.0], colors=["w", "y", "r"])
plt.colorbar(im, ax=ax[0], label="Delta chi^2")
ax[0].set_xlabel("companion/primary flux ratio (TESS band)"); ax[0].set_ylabel("companion Teff (K)")
ax[1].plot(T2s, ul, "o-"); ax[1].set_xlabel("companion Teff (K)"); ax[1].set_ylabel("2-sigma upper limit, TESS flux ratio")
ax2 = ax[1].twiny(); ax2.set_xlim(ax[1].get_xlim()); ax2.set_xticks(T2s[::2]); ax2.set_xticklabels([f"{m:.2f}" for m in M2[::2]]); ax2.set_xlabel("MIST mass (M_sun)")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/sed_binary.png"), dpi=110)
