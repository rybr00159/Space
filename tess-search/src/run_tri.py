"""TRICERATOPS false-positive probability for a candidate using our phase-folded light curve."""
import sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import triceratops.triceratops as tr

tic, sectors, P, t0, dur, lcfile, out = int(sys.argv[1]), [int(s) for s in sys.argv[2].split(",")], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6], sys.argv[7]
z = np.load(lcfile)
ph = ((z["t"] - t0) / P + 0.5) % 1 - 0.5
m = np.abs(ph * P) < 3 * dur
time_d, flux, ferr = ph[m] * P, z["f"][m], np.full(m.sum(), float(z["sigma"]))
order = np.argsort(time_d); time_d, flux, ferr = time_d[order], flux[order], ferr[order]
# bin to 5 minutes for speed
edges = np.arange(time_d.min(), time_d.max() + 5 / 1440, 5 / 1440); idx = np.digitize(time_d, edges)
tb, fb, eb = [], [], []
for k in np.unique(idx):
    s = idx == k
    if s.sum() >= 3: tb.append(time_d[s].mean()); fb.append(flux[s].mean()); eb.append(ferr[s].mean() / np.sqrt(s.sum()))
tb, fb, eb = map(np.array, (tb, fb, eb))
import os as _os
target = tr.target(ID=tic, sectors=sectors, search_radius=int(_os.environ.get("TRI_RADIUS", 10)))
import os
if os.environ.get("TRI_MASS"):  # fill missing TIC stellar parameters for the target (row 0)
    target.stars.loc[0, "mass"] = float(os.environ["TRI_MASS"])
if os.environ.get("TRI_RAD"):
    target.stars.loc[0, "rad"] = float(os.environ["TRI_RAD"])
if os.environ.get("TRI_PLX"):
    target.stars.loc[0, "plx"] = float(os.environ["TRI_PLX"])
print(target.stars.head(5)[["ID", "Tmag", "Teff", "rad", "mass", "sep (arcsec)"]].to_string(), flush=True)
apertures = None
target.calc_depths(tdepth=float(1 - np.min(np.convolve(fb, np.ones(3) / 3, "same"))), all_ap_pixels=apertures)
target.calc_probs(time=tb, flux_0=fb, flux_err_0=float(np.median(eb)), P_orb=P, parallel=False, N=200000, verbose=0)
res = dict(FPP=float(target.FPP), NFPP=float(target.NFPP))
df = target.probs
res["top"] = df.sort_values("prob", ascending=False).head(8)[["ID", "scenario", "prob"]].to_dict("records")
json.dump(res, open(out, "w"), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
