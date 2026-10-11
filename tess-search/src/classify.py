"""Split stage-1 candidates into planet-like and eclipsing-binary-like tracks."""
import numpy as np, pandas as pd

G = 6.674e-11; RSUN = 6.957e8; MSUN = 1.989e30


def tmax_hours(P, rad, mass):
    """Duration of a central transit (k << 1) for a circular orbit, hours."""
    rad = np.where(np.isfinite(rad) & (rad > 0), rad, 1.0); mass = np.where(np.isfinite(mass) & (mass > 0), mass, 1.0)
    a = (G * mass * MSUN * (P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return P * 24 / np.pi * np.arcsin(np.clip(rad * RSUN / a, 0, 1))


def classify(c, tic):
    t = tic.copy(); t["ID"] = t.ID.astype("int64")
    for k in ("Tmag", "Teff", "rad", "mass", "d", "ra", "dec", "contratio"): t[k] = pd.to_numeric(t[k], errors="coerce")
    c = c.merge(t[["ID", "Tmag", "Teff", "rad", "mass", "d", "ra", "dec", "contratio"]], left_on="tic", right_on="ID", how="left")
    c["depth_ppm"] = c.m_depth * 1e6
    c["rp_re"] = np.sqrt(c.m_depth.clip(0)) * c.rad * 109.1
    c["tmax_h"] = tmax_hours(c.P.values, c.rad.values, c.mass.values)
    c["dur_ratio"] = c.dur * 24 / c.tmax_h
    eb = (c.m_oddeven_sig >= 3.5) | (c.m_sec_sig >= 3.5) | (c.depth_ppm > 30000) | (c.rp_re > 25)
    var = (c.m_sec_sig <= -4) | (c.red_ratio > 1.8) | (c.m_oot_rms_ratio > 1.6) | ((c.m_oot_dip_sig > 4) & (c.m_oot_dip_ratio > 0.4))  # brightening at phase 0.5 / red-noise dominated: pulsator or spotted star
    planet = ~eb & ~var & (c.dur_ratio < 2.5) & (c.depth_ppm >= 80)
    c["track"] = np.where(planet, "planet", np.where(eb, "eb", np.where(var, "variable", "other")))
    return c
