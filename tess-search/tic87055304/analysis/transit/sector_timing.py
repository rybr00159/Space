"""Per-sector mean transit-time offsets from a global template fit of each sector's detrended data
(complements the per-transit O-C of ttv.py). For each sector: model = c * (1 - A * s(t - dt)), with dt on a grid
(0.2 min steps over +-40 min) and c, A linear; dt error from Delta chi^2 = 1 with chi^2 scaled by beta^2 (vshape.json).
Offsets are relative to the free-density MCMC ephemeris. Writes sector_timing.json."""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy import stats
import tmodel as tm
from transit_fits import Template

HERE = os.path.dirname(os.path.abspath(__file__))
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th); T0, P = th[0], th[1]
BETA = {int(k): v for k, v in json.load(open(os.path.join(HERE, "vshape.json")))["beta"].items()}
data = tm.load()
out = {}
grid = np.arange(-40, 40.01, 0.2) / 1440
for s in tm.SECTORS:
    d = data[s]
    tau = d["t"] - T0 - np.round((d["t"] - T0) / P) * P
    sig = mc["params"][f"sigma{s}_ppm"]["med"] * 1e-6 * BETA[s]
    chi = []
    for dt in grid:
        sx = 1 - tmpl(tau - dt, d["exptime"])
        X = np.vstack([np.ones_like(sx), -sx]).T
        c, *_ = np.linalg.lstsq(X, d["f"], rcond=None)
        chi.append(np.sum((d["f"] - X @ c) ** 2) / sig ** 2)
    chi = np.array(chi); i = int(np.argmin(chi))
    # parabola around min
    j = slice(max(i - 15, 0), min(i + 16, len(grid)))
    a, b, c0 = np.polyfit(grid[j], chi[j], 2)
    dt0 = -b / (2 * a); err = 1 / np.sqrt(a)
    out[s] = dict(dt_min=float(dt0 * 1440), err_min=float(err * 1440), beta=BETA[s])
v = np.array([out[s]["dt_min"] for s in tm.SECTORS]); e = np.array([out[s]["err_min"] for s in tm.SECTORS])
w = 1 / e ** 2; m = np.sum(v * w) / w.sum()
chi2 = float(np.sum((v - m) ** 2 * w))
# dof: 5 sector offsets relative to an ephemeris fitted to the same data (2 parameters) -> 3
out["chi2_about_ephemeris"] = float(np.sum((v / e) ** 2)); out["dof"] = 3
out["p_value"] = float(stats.chi2.sf(out["chi2_about_ephemeris"], 3))
json.dump(out, open(os.path.join(HERE, "sector_timing.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
