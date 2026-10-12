"""Assemble the observed photometry table for the SED fit (out/sed_photometry.csv).
Values copied from the VizieR/Gaia rows saved by photometry.py / catalogs.py (out/phot_raw, data/).
sigma_fit = sqrt(sigma_cat^2 + floor^2); floors reflect zero-point / system-transformation uncertainties.
"""
import numpy as np, pandas as pd
R = "/home/user/work/tic87/stellar/out/phot_raw/"
g = pd.read_csv("/home/user/work/tic87/data/gaia_dr3_90arcsec.csv").sort_values("sep").iloc[0]
sm = pd.read_csv(R + "SMSS_DR4.csv").iloc[0]
ap = pd.read_csv(R + "APASS9.csv").iloc[0]
ty = pd.read_csv(R + "Tycho2.csv").iloc[0]
tm = pd.read_csv(R + "2MASS.csv").iloc[0]
aw = pd.read_csv(R + "AllWISE.csv").iloc[0]
rows = []
def add(band, mag, err, floor, src, use=True):
    rows.append(dict(band=band, mag=float(mag), err=float(err), floor=floor, sigma=float(np.hypot(err, floor)), source=src, use=use))
# Gaia DR3 (flux errors -> mag errors)
for b, fl, fe, m in [("G", g.phot_g_mean_flux, g.phot_g_mean_flux_error, g.phot_g_mean_mag),
                     ("BP", g.phot_bp_mean_flux, g.phot_bp_mean_flux_error, g.phot_bp_mean_mag),
                     ("RP", g.phot_rp_mean_flux, g.phot_rp_mean_flux_error, g.phot_rp_mean_mag)]:
    add(b, m, 1.0857 * fe / fl, 0.01, "Gaia DR3")
for b in "uvgriz":
    add("SM_" + b, sm[f"{b}PSF"], sm[f"e_{b}PSF"], 0.03 if b in "uv" else 0.02, "SkyMapper DR4 PSF")
add("B", ap["Bmag"], ap["e_Bmag"], 0.05, "APASS DR9")
add("V", ap["Vmag"], ap["e_Vmag"], 0.05, "APASS DR9")
add("gp", ap["g'mag"], ap["e_g'mag"], 0.05, "APASS DR9")
add("rp", ap["r'mag"], ap["e_r'mag"], 0.05, "APASS DR9")
add("ip", ap["i'mag"], ap["e_i'mag"], 0.05, "APASS DR9")
add("BT", ty["BTmag"], ty["e_BTmag"], 0.02, "Tycho-2")
add("VT", ty["VTmag"], ty["e_VTmag"], 0.02, "Tycho-2")
add("J", tm["Jmag"], tm["e_Jmag"], 0.02, "2MASS")
add("H", tm["Hmag"], tm["e_Hmag"], 0.02, "2MASS")
add("Ks", tm["Kmag"], tm["e_Kmag"], 0.02, "2MASS")
add("W1", aw["W1mag"], aw["e_W1mag"], 0.02, "AllWISE")
add("W2", aw["W2mag"], aw["e_W2mag"], 0.02, "AllWISE")
add("W3", aw["W3mag"], aw["e_W3mag"], 0.03, "AllWISE")
df = pd.DataFrame(rows)
df.to_csv("/home/user/work/tic87/stellar/out/sed_photometry.csv", index=False)
print(df.to_string())
