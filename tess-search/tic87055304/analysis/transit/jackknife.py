"""Sector jackknife: compare each drop-one-sector MCMC (mcmc_dropS) with the full free-density fit (mcmc_free).
Shift significance reported two ways:
  (a) shift / sigma_full;
  (b) shift / sqrt(|sigma_drop^2 - sigma_full^2|) -- the appropriate scale for nested data sets (the variance of the
      difference between a subset estimate and the full-set estimate for Gaussian, independent data).
Writes jackknife.json and jackknife.png."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
full = json.load(open(os.path.join(HERE, "mcmc_free.json")))["params"]
keys = ["k", "b", "rho", "a", "inc", "depth_ppm", "T14_h", "T23_h", "P", "T0", "u1", "u2"]
out = {}
fig, axes = plt.subplots(1, 6, figsize=(16, 3.4))
for s in (13, 39, 66, 93, 104):
    p = os.path.join(HERE, f"mcmc_drop{s}.json")
    if not os.path.exists(p):
        continue
    js = json.load(open(p)); d = js["params"]
    row = dict(converged=js["converged"], n_over_tau=js["n_over_tau"], acceptance=js["acceptance"])
    for k in keys:
        sf = 0.5 * (full[k]["lo"] + full[k]["hi"]); sd = 0.5 * (d[k]["lo"] + d[k]["hi"])
        sh = d[k]["med"] - full[k]["med"]
        row[k] = dict(med=d[k]["med"], lo=d[k]["lo"], hi=d[k]["hi"], shift=sh, shift_over_sigma_full=sh / sf,
                      shift_over_sigma_nested=sh / np.sqrt(abs(sd ** 2 - sf ** 2)) if abs(sd - sf) > 1e-12 else np.nan)
    out[f"drop{s}"] = row
json.dump(out, open(os.path.join(HERE, "jackknife.json"), "w"), indent=1)
for i, k in enumerate(["k", "b", "rho", "depth_ppm", "T14_h", "T23_h"]):
    ax = axes[i]
    ax.axhspan(full[k]["med"] - full[k]["lo"], full[k]["med"] + full[k]["hi"], color="0.85")
    ax.axhline(full[k]["med"], color="k", lw=0.8)
    for j, (nm, row) in enumerate(out.items()):
        ax.errorbar(j, row[k]["med"], [[row[k]["lo"]], [row[k]["hi"]]], fmt="o", color=f"C{j}")
    ax.set_xticks(range(len(out))); ax.set_xticklabels([n.replace("drop", "-S") for n in out], fontsize=8)
    ax.set_title(k, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "jackknife.png"), dpi=120)
for nm, row in out.items():
    print(nm, {k: (round(row[k]["med"], 4), round(row[k]["shift_over_sigma_full"], 2), round(row[k]["shift_over_sigma_nested"], 2)) for k in ["k", "b", "rho", "depth_ppm", "T14_h", "T23_h", "P"]})
