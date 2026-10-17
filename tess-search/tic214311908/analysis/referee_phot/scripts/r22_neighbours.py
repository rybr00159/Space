"""Referee r22: census of every Gaia DR3 source within 40" of the M3V that could host the dip: TESS flux in PRF units
(TIC Tmag where available, else Stassun relation, x the PRF/TIC flux ratio 1.107 from r17), the eclipse depth it would
need (15.0 e-/s flux loss, referee host-test value), and its offset from the referee difference-image centroid
(r10: E +0.52 +- 0.89", N -1.36 +- 0.76" from the M3V), as chi2 (2 dof) with statistical errors and with 1.5" added
systematic per axis. Stars needing > 100% depth are impossible hosts. Positions propagated to 2023.
Output: out/neighbours.csv"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from scipy.stats import chi2 as C2
g = pd.read_csv(f"{R}/data/gaia330.csv")
tic = pd.read_csv(f"{R}/data/tic330.csv", dtype={"GAIA": str}); tic = tic[tic.GAIA.notna()]; TT = dict(zip(tic.GAIA.astype(np.int64), tic.Tmag))
D = json.load(open(f"{R}/out/diffimg.json")); K = json.load(open(f"{R}/out/invscale.json"))["prf_over_tic_overall"]
ep = 2023.0 - 2016.0
ra = g.ra + np.nan_to_num(g.pmra) * ep / 3.6e6 / np.cos(np.radians(g.dec)); de = g.dec + np.nan_to_num(g.pmdec) * ep / 3.6e6
m3 = g.source_id == GID["M3V"]
dE = (ra - float(ra[m3].iloc[0])) * np.cos(np.radians(float(de[m3].iloc[0]))) * 3600; dN = (de - float(de[m3].iloc[0])) * 3600
g = g.assign(dE=dE, dN=dN, sep=np.hypot(dE, dN), T=[TT.get(i, np.nan) for i in g.source_id])
g["T"] = np.where(np.isfinite(g["T"]), g["T"], gaia_T(g))
g = g[g.sep < 40].copy()
g["F_prf"] = tflux(g["T"]) * K
g["need_depth"] = 15.0 / g.F_prf
cE, cN, sE, sN = D["centroid_E"], D["centroid_N"], D["centroid_E_err"], D["centroid_N_err"]
g["chi2_stat"] = (g.dE - cE) ** 2 / sE ** 2 + (g.dN - cN) ** 2 / sN ** 2
g["chi2_sys15"] = (g.dE - cE) ** 2 / (sE ** 2 + 2.25) + (g.dN - cN) ** 2 / (sN ** 2 + 2.25)
g["p_sys15"] = C2.sf(g.chi2_sys15, 2)
g["possible"] = g.need_depth <= 1.0
g = g.sort_values("sep")
cols = ["source_id", "phot_g_mean_mag", "bp_rp", "T", "F_prf", "need_depth", "sep", "dE", "dN", "parallax", "ruwe", "chi2_stat", "chi2_sys15", "p_sys15", "possible"]
g[cols].to_csv(f"{R}/out/neighbours.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 220)
print(g[cols].round(3).to_string(index=False))
