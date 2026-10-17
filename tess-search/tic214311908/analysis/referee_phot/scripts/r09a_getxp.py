"""Referee r09a: download Gaia DR3 XP sampled spectra (336-1020 nm) for the field stars of the r01b flux-scale sample,
the M3V, pair A, and the 3" background star. Output: data/ext/xp_sampled.csv (source_id, wavelength nm, flux W/m2/nm)."""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from astroquery.gaia import Gaia
d = pd.read_csv(f"{R}/out/fluxscale_stars.csv")
ids = list(d[d.iso & d.inside & (d.Ttic < 14.5)].source_id.unique()) + [GID["M3V"], GID["pA"], GID["bg3"]]
ids = sorted(set(int(i) for i in ids))
dl = Gaia.load_data(ids=ids, data_release="Gaia DR3", retrieval_type="XP_SAMPLED", data_structure="DATAMODEL_GAIA", format="csv", valid_data=False)
rows = []
for k, v in dl.items():
    for tab in v:
        t = tab.to_pandas() if hasattr(tab, "to_pandas") else tab
        for r in t.itertuples():
            wl = np.arange(336, 1021, 2.0)
            fl = np.array(str(r.flux).strip("()[]").split(","), float); fe = np.array(str(r.flux_error).strip("()[]").split(","), float)
            rows.append(pd.DataFrame(dict(source_id=int(r.source_id), wl=wl, flux=fl, flux_err=fe)))
X = pd.concat(rows); X.to_csv(f"{R}/data/ext/xp_sampled.csv", index=False)
print(len(ids), "requested;", X.source_id.nunique(), "spectra")
