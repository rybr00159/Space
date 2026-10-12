"""Quadratic limb-darkening priors for TESS from Claret (2017, A&A 600, A30; VizieR J/A+A/600/A30, table25:
ATLAS models, quadratic law, least-squares method LSM).

Grid: Teff 5700-7300 K, logg 3.5-4.2 (4.2 interpolated linearly between the 4.0 and 4.5 grid points),
[M/H] = -0.5 and 0.0 (GSP-Phot gives -0.52; TIC assumes solar), microturbulence 2 km/s.
Prior: Gaussian on u1 and u2 with mean = grid mean and sigma = sqrt(grid_sd^2 + 0.05^2)
(0.05 added for the known ~0.05 systematic offset between model and empirical TESS LDCs; e.g. Patel & Espinoza 2022).
Writes ld_prior.json and ld_grid.csv.
"""
import json, os
import numpy as np
from astroquery.vizier import Vizier

HERE = os.path.dirname(os.path.abspath(__file__))
v = Vizier(row_limit=-1, columns=["**"])
t = v.get_catalogs("J/A+A/600/A30/table25")[0].to_pandas()
t.to_csv(os.path.join(HERE, "claret2017_table25_full.csv"), index=False)
sel = t[(t.Teff >= 5700) & (t.Teff <= 7300) & (t.Z.isin([-0.5, 0.0])) & (t.xi == 2.0) & (t.logg.isin([3.5, 4.0, 4.5]))]
rows = []
for (Z, Teff), g in sel.groupby(["Z", "Teff"]):
    g = g.set_index("logg")
    if not all(x in g.index for x in (3.5, 4.0, 4.5)):
        continue
    for lg in (3.5, 4.0, 4.2):
        if lg == 4.2:
            a = 0.6 * g.loc[4.0, "aLSM"] + 0.4 * g.loc[4.5, "aLSM"]; b = 0.6 * g.loc[4.0, "bLSM"] + 0.4 * g.loc[4.5, "bLSM"]
        else:
            a, b = g.loc[lg, "aLSM"], g.loc[lg, "bLSM"]
        rows.append(dict(Z=Z, Teff=Teff, logg=lg, u1=float(a), u2=float(b)))
import pandas as pd
grid = pd.DataFrame(rows)
grid.to_csv(os.path.join(HERE, "ld_grid.csv"), index=False)
res = dict(source="Claret 2017 J/A+A/600/A30 table25 (ATLAS, quadratic, LSM, xi=2 km/s)",
           Teff_range=[5700, 7300], logg_range=[3.5, 4.2], MH=[-0.5, 0.0], n_grid=len(grid),
           u1_mean=grid.u1.mean(), u1_sd=grid.u1.std(), u1_min=grid.u1.min(), u1_max=grid.u1.max(),
           u2_mean=grid.u2.mean(), u2_sd=grid.u2.std(), u2_min=grid.u2.min(), u2_max=grid.u2.max())
res["u1_prior_sigma"] = float(np.hypot(res["u1_sd"], 0.05)); res["u2_prior_sigma"] = float(np.hypot(res["u2_sd"], 0.05))
# also the two stellar solutions separately
for lab, T, lg in (("TIC_7276_logg4.2", 7276, 4.2), ("GSPPhot_5696_logg3.5", 5696, 3.5)):
    g = grid[(grid.logg == lg)]
    i = (g.Teff - T).abs().sort_values().index[:2]
    res[lab] = dict(u1=float(g.loc[i, "u1"].mean()), u2=float(g.loc[i, "u2"].mean()))
json.dump(res, open(os.path.join(HERE, "ld_prior.json"), "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
