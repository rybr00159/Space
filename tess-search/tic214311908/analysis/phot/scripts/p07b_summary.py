"""Cross-talk matrices from the pixel-level injections (p07): mean and scatter (over 7 fake ephemerides) of the
recovered dip on each light curve divided by the injected dip, per sector and host, and the overall average."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
KEYS = ["B:M3V", "B:Mira", "Bp:pair", "A:M3V", "A:pair", "A:Mira"]
rows = []
for s in SECTORS:
    J = json.load(open(f"{ROOT}/out/inject_s{s}.json"))
    for r in J["rows"]:
        rows.append(dict(sector=s, host=r["host"], eph=r["eph"], **{k: r[f"rec_{k}"] / J["delta"] for k in KEYS}))
d = pd.DataFrame(rows)
g = d.groupby(["host", "sector"])[KEYS]
mean, std = g.mean(), g.std()
pd.set_option("display.width", 200)
print("mean response (recovered / injected):\n", mean.round(3).to_string())
print("scatter over fake ephemerides:\n", std.round(3).to_string())
allm = d.groupby("host")[KEYS].mean(); alls = d.groupby("host")[KEYS].std()
print("all sectors:\n", allm.round(3).to_string(), "\n", alls.round(3).to_string())
save_json(dict(per_sector_mean={f"{h}|{s}": mean.loc[(h, s)].to_dict() for h, s in mean.index},
               per_sector_std={f"{h}|{s}": std.loc[(h, s)].to_dict() for h, s in std.index},
               all_mean={h: allm.loc[h].to_dict() for h in allm.index}, all_std={h: alls.loc[h].to_dict() for h in alls.index}),
          f"{ROOT}/out/inject_summary.json")
d.to_csv(f"{ROOT}/out/inject_all.csv", index=False, float_format="%.4f")
