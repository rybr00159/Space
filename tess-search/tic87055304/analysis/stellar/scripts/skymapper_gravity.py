"""Empirical, parallax-independent luminosity-class test with SkyMapper DR4 u, v photometry (Balmer jump / Ca HK
region): compare the target's (u-v), (v-g), (u-g) at its Gaia BP-RP with field stars on the same sightline whose
luminosity class is known from good Gaia parallaxes (plx/err > 10, RUWE < 1.4):
  dwarfs     : M_G > 3.6           (main sequence / near turn-off)
  subgiants  : 1.4 < M_G < 3.0     (the GSP-Phot MARCS solution has M_G = 2.04)
Colours are fitted as quadratic functions of BP-RP within 0.70 < BP-RP < 1.05 per class; the target's offset is
evaluated as a likelihood with the class scatter + target error. Reddening along this sightline is nearly flat
beyond ~300 pc (reddening_distance.py), so the two classes are compared without individual dereddening; a variant
using Zhang+23 extinction (E) to restrict to 0.04 < E < 0.12 is also given.
Output: out/skymapper_gravity.json, fig/skymapper_gravity.png
"""
import os, json, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
OUT = "/home/user/work/tic87/stellar/out"; FIG = "/home/user/work/tic87/stellar/fig"
c0 = SkyCoord(272.606477 * u.deg, -43.995466 * u.deg)
TID = 6720986204578421888
RAD = 60 * u.arcmin


def vq(cat, cols, radius=RAD, extra=None):
    V = Vizier(columns=cols, row_limit=-1, timeout=900, column_filters=extra or {})
    for i in range(4):
        try:
            r = V.query_region(c0, radius=radius, catalog=cat)
            return r[0].to_pandas() if len(r) else pd.DataFrame()
        except Exception as e:
            print(cat, "retry", e); time.sleep(20)
    return pd.DataFrame()

f1 = os.path.join(OUT, "reddening/smss_dr4_60arcmin.csv")
if os.path.exists(f1): sm = pd.read_csv(f1)
else:
    sm = vq("II/379/smssdr4", ["ObjectId", "uPSF", "e_uPSF", "vPSF", "e_vPSF", "gPSF", "e_gPSF", "iPSF", "e_iPSF", "flags", "ClassStar", "Gaia3Id1", "Gaia3Dist1", "NGaia3-15"],
            extra={"gPSF": "<16"})
    sm.to_csv(f1, index=False)
f2 = os.path.join(OUT, "reddening/gaia_60arcmin_bright.csv")
if os.path.exists(f2): g = pd.read_csv(f2)
else:
    g = vq("I/355/gaiadr3", ["Source", "Plx", "e_Plx", "RUWE", "Gmag", "BP-RP", "A0", "AG"], extra={"Gmag": "<16"})
    g.to_csv(f2, index=False)
f3 = os.path.join(OUT, "reddening/zhang23_60arcmin.csv")
if os.path.exists(f3): z = pd.read_csv(f3)
else:
    z = vq("J/MNRAS/524/1855/xpparams", ["GaiaDR3", "Ext", "e_Ext", "logg", "Teff", "[Fe/H]", "Quality"], extra={})
    z.to_csv(f3, index=False)
print(len(sm), len(g), len(z))
sm["flags"] = pd.to_numeric(sm["flags"].astype(str), errors="coerce")
d = sm.merge(g, left_on="Gaia3Id1", right_on="Source").merge(z, left_on="Source", right_on="GaiaDR3", how="left")
d["u_v"] = d.uPSF - d.vPSF; d["v_g"] = d.vPSF - d.gPSF; d["u_g"] = d.uPSF - d.gPSF; d["g_i"] = d.gPSF - d.iPSF
d["MG"] = d.Gmag + 5 * np.log10(d.Plx / 100.)
tgt = d[d.Source == TID].iloc[0]
sel = (d.Plx / d.e_Plx > 10) & (d.RUWE < 1.4) & (d["BP-RP"] > 0.70) & (d["BP-RP"] < 1.05) & (d.e_uPSF < 0.04) & (d.e_vPSF < 0.04) \
      & (d["flags"] < 4) & (d.Gaia3Dist1 < 1.0) & (d.Source != TID) & (d.Gmag < 15)
classes = {"dwarf": sel & (d.MG > 3.6), "subgiant": sel & (d.MG > 1.4) & (d.MG < 3.0)}
res = dict(target=dict(BP_RP=float(tgt["BP-RP"]), u_v=float(tgt.u_v), v_g=float(tgt.v_g), u_g=float(tgt.u_g),
                       e_u=float(tgt.e_uPSF), e_v=float(tgt.e_vPSF), e_g=float(tgt.e_gPSF)))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, col in zip(axs, ["u_v", "v_g", "u_g"]):
    out = {}
    tv = tgt[col]
    te = {"u_v": np.hypot(tgt.e_uPSF, tgt.e_vPSF), "v_g": np.hypot(tgt.e_vPSF, tgt.e_gPSF), "u_g": np.hypot(tgt.e_uPSF, tgt.e_gPSF)}[col]
    for name, m in classes.items():
        s = d[m]
        for it in range(3):  # quadratic fit with 3-sigma clipping
            cf = np.polyfit(s["BP-RP"], s[col], 2); r = s[col] - np.polyval(cf, s["BP-RP"])
            sd = 1.4826 * np.median(np.abs(r - np.median(r))); s = s[np.abs(r) < 3 * sd]
        pred = np.polyval(cf, tgt["BP-RP"]); sc = sd
        lnl = -0.5 * (tv - pred) ** 2 / (sc ** 2 + te ** 2) - 0.5 * np.log(2 * np.pi * (sc ** 2 + te ** 2))
        out[name] = dict(n=int(len(s)), pred=float(pred), scatter=float(sc), offset_sigma=float((tv - pred) / np.hypot(sc, te)), lnL=float(lnl))
        ax.scatter(d[m]["BP-RP"], d[m][col], s=4, alpha=0.5, label=f"{name} (N={len(s)})")
        xx = np.linspace(0.7, 1.05, 50); ax.plot(xx, np.polyval(cf, xx), "-")
    out["lnL_dwarf_minus_subgiant"] = out["dwarf"]["lnL"] - out["subgiant"]["lnL"]
    out["class_separation_mag"] = out["dwarf"]["pred"] - out["subgiant"]["pred"]
    res[col] = out
    ax.errorbar(tgt["BP-RP"], tv, yerr=te, fmt="r*", ms=14, label="TIC 87055304")
    ax.set_xlabel("Gaia BP-RP"); ax.set_ylabel(col.replace("_", "-") + " (SkyMapper DR4 PSF)"); ax.legend(fontsize=7)
    ax.set_title(f"{col}: dlnL(dwarf-subgiant) = {out['lnL_dwarf_minus_subgiant']:.2f}")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "skymapper_gravity.png"), dpi=110)
json.dump(res, open(os.path.join(OUT, "skymapper_gravity.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
