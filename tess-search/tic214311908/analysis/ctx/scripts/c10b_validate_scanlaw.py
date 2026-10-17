"""Validate the scanning-law crossing times of c10 with real Gaia DR3 epoch photometry of the Mira (17.8" from the M3V,
so its FOV transits are the same events): match every DR3 G-band observation time to the nearest predicted crossing.
Output: out/gaia_scanlaw_validation.json, lc/gaia_epoch/mira_dr3_epoch_G.csv
"""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from astroquery.gaia import Gaia
dl = Gaia.load_data(ids=[5925373268528684288], data_release="Gaia DR3", retrieval_type="EPOCH_PHOTOMETRY", data_structure="DATAMODEL_STANDARD", format="votable")
key = list(dl.keys())[0]; t = dl[key][0].to_table()
df = t.to_pandas(); df.to_csv("lc/gaia_epoch/mira_dr3_epoch_G.csv", index=False)
g = df[np.isfinite(df.g_transit_time)].copy(); g["mag"] = g.g_transit_mag
obs_tcb = g.g_transit_time.values + 2455197.5
pred = pd.read_csv("out/gaia_fov_transits.csv")
dt = np.array([np.min(np.abs(pred.bjd_tcb.values - o)) * 86400 for o in obs_tcb])
dr3_end = 2457901.5
npred = int(np.sum(pred.bjd_tcb < dr3_end))
res = dict(n_dr3_G_obs=int(len(g)), n_pred_in_dr3_window=npred, match_dt_seconds_median=float(np.median(dt)),
           match_dt_seconds_max=float(np.max(dt)), n_obs_matched_within_60s=int(np.sum(dt < 60)),
           n_pred_without_dr3_obs=int(np.sum([np.min(np.abs(obs_tcb - p)) * 86400 > 60 for p in pred.bjd_tcb[pred.bjd_tcb < dr3_end]])),
           g_mag_range=[float(g.mag.min()), float(g.mag.max())], rejected_by_photometry=int(g.rejected_by_photometry.astype(bool).sum()))
json.dump(res, open("out/gaia_scanlaw_validation.json", "w"), indent=1, default=float); print(res)
