"""Validate the synthetic photometry (sedlib.band_mean + SVO zero points) with CALSPEC Vega (alpha_lyr_stis_011)
and the solar reference spectrum (sun_reference_stis_002). Vega-system bands should give Vega ~ 0.00-0.03 mag;
solar colours are compared with Casagrande & VandenBerg (2018) / Willmer (2018) values.
Output: out/synphot_validation.csv"""
import numpy as np, pandas as pd
from astropy.io import fits
import sedlib
rows = []
for name, fn in [("Vega", "alpha_lyr_stis_011.fits"), ("Sun", "sun_reference_stis_002.fits")]:
    d = fits.open("/home/user/work/tic87/stellar/cache/calspec/" + fn)[1].data
    wl = np.array(d["WAVELENGTH"], float); fl = np.array(d["FLUX"], float)
    for b in sedlib.FILTERS:
        F = sedlib.FILT[b]
        if F["w"][-1] > wl.max() or F["w"][0] < wl.min(): mag = np.nan
        else: mag = -2.5 * np.log10(sedlib.band_mean(wl, fl, F) / F["zp_flam"])
        rows.append(dict(star=name, band=b, mag=mag))
df = pd.DataFrame(rows).pivot(index="band", columns="star", values="mag")
# solar absolute magnitudes at 10 pc: m(1 AU) + 31.572
df["Sun_abs"] = df["Sun"] + 31.5721
df.to_csv("/home/user/work/tic87/stellar/out/synphot_validation.csv")
print(df.round(3).to_string())
print("Sun BP-RP", df.loc["BP", "Sun"] - df.loc["RP", "Sun"], "(CV18: 0.82)  G-Ks", df.loc["G", "Sun"] - df.loc["Ks", "Sun"], "(CV18: 1.53 approx)",
      " J-Ks", df.loc["J", "Sun"] - df.loc["Ks", "Sun"], "(Willmer18: 0.37)")
