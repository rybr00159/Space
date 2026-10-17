"""Known transiting planets with 4-7 R_earth and P < 10 d around M dwarfs (Teff < 3900 K), from the NASA Exoplanet
Archive Planetary Systems Composite table (pscomppars, TAP), plus the same for 3900-4000 K (K7/M0 boundary).
Output: out/known_mdwarf_neptunes.csv, out/known_mdwarf_neptunes.json (counts by Teff and period bins)."""
import io, json, requests
import pandas as pd, numpy as np
q = ("select pl_name,hostname,disc_year,disc_facility,pl_orbper,pl_rade,pl_radeerr1,pl_radeerr2,pl_bmasse,pl_bmasseerr1,pl_eqt,"
     "st_teff,st_rad,st_mass,sy_dist,sy_tmag,tran_flag,pl_rade_reflink "
     "from pscomppars where tran_flag=1 and st_teff<4000 and pl_rade>=3.5 and pl_rade<=8 and pl_orbper<15")
r = requests.get("https://exoplanetarchive.ipac.caltech.edu/TAP/sync", params={"query": q, "format": "csv"}, timeout=300)
r.raise_for_status()
d = pd.read_csv(io.StringIO(r.text)).sort_values("st_teff")
d["pl_refname"] = d.pl_rade_reflink.astype(str).str.replace(r"<[^>]+>", "", regex=True).str.strip()
d.to_csv("out/known_mdwarf_neptunes.csv", index=False)
core = d[(d.pl_rade >= 4) & (d.pl_rade <= 7) & (d.pl_orbper < 10)]
out = dict(query=q, n_rows=len(d),
           n_4to7_P10_teff_lt_3900=int((core.st_teff < 3900).sum()), n_4to7_P10_teff_lt_4000=int(len(core)),
           n_4to7_P3to10_teff_lt_3900=int(((core.st_teff < 3900) & (core.pl_orbper >= 3)).sum()),
           n_4to7_P10_teff_lt_3500=int((core.st_teff < 3500).sum()),
           n_4to7_P10_mass_lt_0p45=int((core.st_mass < 0.45).sum()),
           planets=core[["pl_name", "st_teff", "st_mass", "pl_rade", "pl_orbper", "pl_bmasse", "pl_eqt", "sy_dist", "disc_year", "pl_refname"]].round(3).to_dict("records"))
json.dump(out, open("out/known_mdwarf_neptunes.json", "w"), indent=1, default=str)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 40)
print(core[["pl_name", "st_teff", "st_mass", "pl_rade", "pl_orbper", "pl_bmasse", "pl_eqt", "sy_dist", "disc_year", "pl_refname"]].round(3).to_string(index=False))
print({k: v for k, v in out.items() if k.startswith("n_")})
