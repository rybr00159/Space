"""COPY of /home/user/work/src/run_tri.py (unchanged logic) with ONE addition: after calc_depths(), when env TRI_ONTARGET=1,
target.remove_nearby_stars() drops every star except the target, i.e. the FPP conditional on the signal being on
TIC 87055315 (as established by the pixel-level localisation; the analogue of an on-target ground-based detection).
Only the on-target scenarios remain: TP, EB, EBx2P, PTP, PEB, PEBx2P, STP, SEB, SEBx2P (+ DTP/DEB/BTP/BEB background
scenarios, which TRICERATOPS keeps for unresolved background stars). Same arguments and env variables as run_tri.py."""
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
# transit depth in the aperture light curve. Earlier versions used np.convolve(..., "same"), whose zero-padded
# edges gave ~0.33 regardless of the data; that made every nearby star need a >100 % eclipse, so TRICERATOPS
# dropped all nearby false-positive scenarios (NFPP = 0). Use the 'valid' smoothing, or an explicit TRI_DEPTH.
tdepth = float(os.environ["TRI_DEPTH"]) if os.environ.get("TRI_DEPTH") else float(1 - np.min(np.convolve(fb, np.ones(3) / 3, "valid")))
print("tdepth used:", tdepth, flush=True)
target.calc_depths(tdepth=tdepth, all_ap_pixels=apertures)
if os.environ.get("TRI_ONTARGET") == "1":
    target.remove_nearby_stars()
    print(target.stars[["ID", "Tmag", "fluxratio", "tdepth"]].to_string(), flush=True)
target.calc_probs(time=tb, flux_0=fb, flux_err_0=float(np.median(eb)), P_orb=P, parallel=False, N=200000, verbose=0)
res = dict(FPP=float(target.FPP), NFPP=float(target.NFPP))
df = target.probs
res["top"] = df.sort_values("prob", ascending=False).head(8)[["ID", "scenario", "prob"]].to_dict("records")
res["all"] = df.sort_values("prob", ascending=False)[["ID", "scenario", "prob"] + [c for c in ("M_s", "R_s", "P_orb", "inc", "b", "R_p", "M_EB", "R_EB", "fluxratio_EB", "fluxratio_comp") if c in df.columns]].to_dict("records")
json.dump(res, open(out, "w"), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
