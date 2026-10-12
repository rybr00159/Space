"""Empirical extinction-vs-distance profile along the sightline of TIC 87055304 from neighbouring stars with good
parallaxes, using the SAME catalogue/method as the target's own extinction estimate (so method systematics cancel).
Catalogues: Zhang+23 (XP forward model, E ~ E(B-V)), SHBoost (Khalatyan+24, A_V), StarHorse2021 (A_V),
Gaia DR3 GSP-Phot (A0, from the 90" + wider Gaia query), Gaia DR3 PHOENIX-library GSP-Phot is not in gaia_source.
For each catalogue: neighbours within R arcmin with Gaia DR3 parallax_over_error > 10, RUWE < 1.4, d < 2.5 kpc.
A running median A(d) is computed and the target's A is converted into a 'reddening distance' range.
Also the total-column SFD / Planck reddening is fetched from IRSA.
Output: out/reddening/*.csv, out/reddening/summary.json, fig/reddening_distance.png
"""
import os, json, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
OUT = "/home/user/work/tic87/stellar/out/reddening"; os.makedirs(OUT, exist_ok=True)
FIG = "/home/user/work/tic87/stellar/fig"
c0 = SkyCoord(272.606477 * u.deg, -43.995466 * u.deg)
TID = 6720986204578421888
RAD = 30 * u.arcmin


def vq(cat, cols, radius=RAD):
    V = Vizier(columns=cols, row_limit=-1, timeout=600)
    for i in range(4):
        try:
            r = V.query_region(c0, radius=radius, catalog=cat)
            return r[0].to_pandas() if len(r) else pd.DataFrame()
        except Exception as e:
            print(cat, "retry", e); time.sleep(15)
    return pd.DataFrame()

cache = os.path.join(OUT, "gaia_30arcmin.csv")
if os.path.exists(cache):
    g = pd.read_csv(cache)
else:
    g = vq("I/355/gaiadr3", ["Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "RUWE", "Gmag", "BP-RP", "A0", "AG", "E(BP-RP)", "Dist", "Teff", "logg"])
    g.to_csv(cache, index=False)
print("gaia rows", len(g))
g = g.rename(columns={"Source": "src"})
good = (g["Plx"] / g["e_Plx"] > 10) & (g["RUWE"] < 1.4) & (g["Plx"] > 0.4)
g["dplx"] = 1000 / g["Plx"]

cats = {
    "Zhang23": ("J/MNRAS/524/1855/xpparams", ["GaiaDR3", "Ext", "e_Ext", "Quality"], "GaiaDR3", "Ext"),
    "SHBoost": ("V/160/shboost", ["GaiaDR3", "AV", "AVmean", "s_AVmean"], "GaiaDR3", "AV"),
    "StarHorse21": ("I/354/starhorse2021", ["Source", "AV50", "AV16", "AV84", "dist50"], "Source", "AV50"),
}
res = {}
fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
target_vals = {}
for ax, (name, (cat, cols, idc, acol)) in zip(axs, cats.items()):
    f = os.path.join(OUT, f"{name}_30arcmin.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
    else:
        d = vq(cat, cols); d.to_csv(f, index=False)
    if len(d) == 0:
        res[name] = "query failed"; continue
    d = d.rename(columns={idc: "src"})
    tv = d.loc[d["src"] == TID, acol]
    m = d.merge(g[good][["src", "dplx", "Plx", "e_Plx"]], on="src")
    if name == "Zhang23":
        m = m[m["Quality"] == 0] if "Quality" in m else m
    m = m[np.isfinite(m[acol]) & (m["dplx"] < 2500)]
    edges = np.array([0, 150, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500])
    mid, med, p16, p84, nn = [], [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = m[(m["dplx"] >= lo) & (m["dplx"] < hi)][acol]
        if len(s) >= 5:
            mid.append(0.5 * (lo + hi)); med.append(s.median()); p16.append(s.quantile(0.16)); p84.append(s.quantile(0.84)); nn.append(len(s))
    prof = pd.DataFrame(dict(d_mid=mid, med=med, p16=p16, p84=p84, n=nn))
    prof.to_csv(os.path.join(OUT, f"{name}_profile.csv"), index=False)
    tval = float(tv.iloc[0]) if len(tv) else np.nan
    target_vals[name] = tval
    # reddening distance: where the running median crosses the target value (+-1 sigma of target from the scatter)
    out = dict(target=tval, n_neigh=int(len(m)), profile=prof.to_dict(orient="list"))
    for lab, dd in [("A_at_356pc", 356), ("A_at_464pc", 464), ("A_at_843pc", 843), ("A_at_998pc", 998)]:
        out[lab] = dict(med=float(np.interp(dd, prof.d_mid, prof.med)), p16=float(np.interp(dd, prof.d_mid, prof.p16)),
                        p84=float(np.interp(dd, prof.d_mid, prof.p84)))
    # fraction of neighbours in each distance bin with A <= target (percentile of target within bin)
    pct = {}
    for lo, hi in [(300, 420), (420, 600), (750, 950), (900, 1100)]:
        s = m[(m["dplx"] >= lo) & (m["dplx"] < hi)][acol]
        pct[f"{lo}-{hi}"] = dict(n=int(len(s)), pct_below_target=float((s <= tval).mean()) if len(s) else None,
                                 median=float(s.median()) if len(s) else None)
    out["target_percentile_in_bins"] = pct
    res[name] = out
    ax.scatter(m["dplx"], m[acol], s=3, c="0.6", label=f"neighbours (N={len(m)})")
    ax.plot(prof.d_mid, prof.med, "k-"); ax.fill_between(prof.d_mid, prof.p16, prof.p84, color="k", alpha=0.15)
    ax.axhline(tval, color="r", label=f"target {acol}={tval:.3f}")
    for dd, ls in [(356, "--"), (998, ":")]:
        ax.axvline(dd, color="b", ls=ls)
    ax.set_xlabel("1/parallax (pc)"); ax.set_ylabel(acol); ax.set_title(name); ax.set_xlim(0, 2000); ax.legend(fontsize=7)

# Gaia DR3 GSP-Phot A0 for neighbours vs parallax distance (target's own A0 from MARCS lib = 0.166 at 998 pc,
# PHOENIX lib = 0.309 at 334 pc, A lib = 0.390 at 843 pc)
ax = axs[3]
m = g[good & np.isfinite(g["A0"])]
ax.scatter(m["dplx"], m["A0"], s=3, c="0.6")
edges = np.array([0, 150, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500])
mid = []; med = []; p16 = []; p84 = []
for lo, hi in zip(edges[:-1], edges[1:]):
    s = m[(m["dplx"] >= lo) & (m["dplx"] < hi)]["A0"]
    if len(s) >= 5: mid.append(0.5 * (lo + hi)); med.append(s.median()); p16.append(s.quantile(.16)); p84.append(s.quantile(.84))
ax.plot(mid, med, "k-"); ax.fill_between(mid, p16, p84, color="k", alpha=0.15)
for v, lab, c in [(0.1661, "MARCS 998pc", "r"), (0.3092, "PHOENIX 334pc", "m"), (0.3903, "A-lib 843pc", "orange")]:
    ax.axhline(v, color=c, label=lab)
ax.set_xlim(0, 2000); ax.set_title("Gaia GSP-Phot A0"); ax.set_xlabel("1/parallax (pc)"); ax.legend(fontsize=7)
res["GSPPhot_A0_profile"] = dict(d_mid=mid, med=med, p16=p16, p84=p84, n=int(len(m)))
fig.tight_layout(); fig.savefig(os.path.join(FIG, "reddening_distance.png"), dpi=110)

# total column from IRSA dust service
try:
    from astroquery.ipac.irsa.irsa_dust import IrsaDust
    t = IrsaDust.get_query_table(c0, section="ebv")
    res["IRSA_total_column"] = {k: float(t[k][0]) for k in t.colnames if ("SandF" in k or "SFD" in k) and ("ref" in k and "RA" not in k and "Dec" not in k and "sys" not in k or "mean" in k or "std" in k)}
except Exception as e:
    res["IRSA_total_column"] = f"failed: {e}"
json.dump(res, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "profile"}) for k, v in res.items()}, indent=1, default=float))
