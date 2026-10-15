"""Variant of the FIXED runner /home/user/work/src/run_tri.py (identical light-curve handling, folding, 5-min binning,
depth estimate and calc_depths/calc_probs calls). Additions, all optional via environment variables:
  TRI_DROP_NEARBY=1 : call target.remove_nearby_stars() after calc_depths (resolved neighbours excluded by TESS localisation)
  TRI_CC=<file>, TRI_CC_FILT=<TESS|Vis|J|H|K> : contrast curve passed to calc_probs
  TRI_N : number of MC draws (default 200000, as in the fixed runner); TRI_SEED : numpy seed
  TRI_TEFF : override target Teff
  TRI_DROP : comma-separated scenarios passed as drop_scenario (e.g. planet scenarios excluded by the odd/even eclipse difference)
Writes the FULL scenario table (all scenarios, lnZ, best-fit parameters) instead of only the top 8."""
import sys, json, warnings, os
import numpy as np
warnings.filterwarnings("ignore")
import triceratops.triceratops as tr
tic, sectors, P, t0, dur, lcfile, out = int(sys.argv[1]), [int(s) for s in sys.argv[2].split(",")], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6], sys.argv[7]
if os.environ.get("TRI_SEED"): np.random.seed(int(os.environ["TRI_SEED"]))
z = np.load(lcfile)
ph = ((z["t"] - t0) / P + 0.5) % 1 - 0.5
m = np.abs(ph * P) < 3 * dur
time_d, flux, ferr = ph[m] * P, z["f"][m], np.full(m.sum(), float(z["sigma"]))
order = np.argsort(time_d); time_d, flux, ferr = time_d[order], flux[order], ferr[order]
edges = np.arange(time_d.min(), time_d.max() + 5 / 1440, 5 / 1440); idx = np.digitize(time_d, edges)
tb, fb, eb = [], [], []
for k in np.unique(idx):
    s = idx == k
    if s.sum() >= 3: tb.append(time_d[s].mean()); fb.append(flux[s].mean()); eb.append(ferr[s].mean() / np.sqrt(s.sum()))
tb, fb, eb = map(np.array, (tb, fb, eb))
np.savez(out.replace('.json','_binned.npz'), t=tb, f=fb, e=eb)
target = tr.target(ID=tic, sectors=sectors, search_radius=int(os.environ.get("TRI_RADIUS", 10)))
if os.environ.get("TRI_MASS"): target.stars.loc[0, "mass"] = float(os.environ["TRI_MASS"])
if os.environ.get("TRI_RAD"): target.stars.loc[0, "rad"] = float(os.environ["TRI_RAD"])
if os.environ.get("TRI_PLX"): target.stars.loc[0, "plx"] = float(os.environ["TRI_PLX"])
if os.environ.get("TRI_TEFF"): target.stars.loc[0, "Teff"] = float(os.environ["TRI_TEFF"])
print(target.stars.head(5)[["ID", "Tmag", "Teff", "rad", "mass", "plx", "sep (arcsec)"]].to_string(), flush=True)
tdepth = float(os.environ["TRI_DEPTH"]) if os.environ.get("TRI_DEPTH") else float(1 - np.min(np.convolve(fb, np.ones(3) / 3, "valid")))
print("tdepth used:", tdepth, "flux_err_0:", float(np.median(eb)), "nbins:", len(tb), flush=True)
target.calc_depths(tdepth=tdepth, all_ap_pixels=None)
if os.environ.get("TRI_DROP_NEARBY") == "1": target.remove_nearby_stars()
kw = {}
if os.environ.get("TRI_CC"): kw = dict(contrast_curve_file=os.environ["TRI_CC"], filt=os.environ.get("TRI_CC_FILT", "TESS"))
N = int(os.environ.get("TRI_N", 200000))
if os.environ.get("TRI_DROP"): kw["drop_scenario"] = os.environ["TRI_DROP"].split(",")
target.calc_probs(time=tb, flux_0=fb, flux_err_0=float(np.median(eb)), P_orb=P, parallel=False, N=N, verbose=0, **kw)
df = target.probs.copy(); df["lnZ"] = target.lnZ; df["fluxratio_EB"] = target.fluxratio_EB; df["fluxratio_comp"] = target.fluxratio_comp
res = dict(FPP=float(target.FPP), NFPP=float(target.NFPP), tdepth=tdepth, N=N, P=P, sectors=sectors, lcfile=lcfile,
           env={k: v for k, v in os.environ.items() if k.startswith("TRI_")}, scenarios=df.to_dict("records"))
json.dump(res, open(out, "w"), indent=1, default=str)
print(df.sort_values("prob", ascending=False).head(12)[["ID", "scenario", "prob", "lnZ", "M_s", "R_s", "inc", "b", "M_EB", "R_EB"]].to_string())
print("FPP", target.FPP, "NFPP", target.NFPP)
