"""NEOWISE-R (2013-2024) and AllWISE multi-epoch (2010-2011) single-exposure photometry of TIC 214311908 (W1/W2), folded on
the TESS ephemeris. The Mira (W1 ~ 4-5, saturated) is 17.8" away = ~3 W1 FWHM, so contamination flags are reported.
Times converted to BJD_TDB for the M3V. Output: out/neowise_fold.json, lc/neowise_m3v.csv, lc/allwise_mep_m3v.csv
"""
import os, sys, json, io, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, requests
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import common as C

ra16, de16, pmra, pmde = 260.9696736991784, -51.323821211802624, -113.67115, 85.54355
def pos_at(year):
    dt = year - 2016.0
    return ra16 + pmra * dt / 3.6e6 / np.cos(np.radians(de16)), de16 + pmde * dt / 3.6e6
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
def tap(q):
    r = requests.get(TAP, params={"QUERY": q, "FORMAT": "csv"}, timeout=600); r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))
out = {}
ra, de = pos_at(2019.0)
q1 = (f"SELECT ra, dec, mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, w1rchi2, w2rchi2, cc_flags, ph_qual, qual_frame, qi_fact, saa_sep, moon_masked, nb, na, sso_flg "
      f"FROM neowiser_p1bs_psd WHERE CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{de},{6/3600}))=1")
q2 = (f"SELECT ra, dec, mjd, w1mpro_ep, w1sigmpro_ep, w2mpro_ep, w2sigmpro_ep, w1rchi2_ep, w2rchi2_ep, nb, na, qi_fact, saa_sep, moon_masked, cc_flags "
      f"FROM allwise_p3as_mep WHERE CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{pos_at(2010.5)[0]},{pos_at(2010.5)[1]},{6/3600}))=1")
for name, q, w1, e1, w2, e2 in (("neowise", q1, "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro"), ("allwise_mep", q2, "w1mpro_ep", "w1sigmpro_ep", "w2mpro_ep", "w2sigmpro_ep")):
    d = tap(q)
    if len(d) == 0:
        out[name] = "no detections"; continue
    t = Time(d.mjd.values, format="mjd", scale="utc")
    yr = 2000 + (t.jd - 2451545.0) / 365.25
    rr, dd = pos_at(np.median(yr))
    # keep only detections near the M3V at that epoch (and not the Mira, 17.8" away)
    sep = SkyCoord(d.ra.values, d.dec.values, unit="deg").separation(SkyCoord(rr, dd, unit="deg")).arcsec
    d["sep"] = sep; d = d[sep < 3.0].copy(); t = Time(d.mjd.values, format="mjd", scale="utc")
    # barycentric correction for the M3V (WISE is in low-Earth orbit: geocentric is adequate at the 20-ms level)
    tgt = SkyCoord(rr, dd, unit="deg")
    ltt = t.light_travel_time(tgt, kind="barycentric", location=EarthLocation.from_geocentric(0, 0, 0, unit=u.m))
    d["btjd"] = (t.tdb + ltt).jd - 2457000.0
    d.to_csv(f"lc/{name}_m3v.csv", index=False)
    res = dict(n=len(d), years=[float(yr.min()), float(yr.max())], median_sep_arcsec=float(np.median(d.sep)),
               cc_flags_counts=d.cc_flags.astype(str).value_counts().to_dict())
    tmpl = C.template(d.btjd.values, 7.7 / 86400)        # WISE exposures 7.7 s (W1)
    ep = C.epoch(d.btjd.values); sigT = np.sqrt(C.T0_ERR**2 + (ep * C.P_ERR)**2 + 2 * ep * C.EPH["cov"][0][1])
    ph = ((d.btjd.values - C.T0) / C.P + 0.5) % 1 - 0.5
    res["n_in_transit"] = int(np.sum(tmpl > 0.5)); res["timing_sigma_min_at_epochs"] = [float(sigT.min() * 1440), float(sigT.max() * 1440)]
    for band, m, e in (("W1", w1, e1), ("W2", w2, e2)):
        ok = np.isfinite(d[m]) & np.isfinite(d[e])
        mag = d[m].values[ok]; err = d[e].values[ok]; tm = tmpl[ok]
        # per-visit (6-month season) median removal: the M3V is constant to first order, but this absorbs season-to-season offsets
        season = np.round((d.btjd.values[ok] - d.btjd.values[ok].min()) / 182.6).astype(int)
        flux = 10 ** (-0.4 * mag); rel = np.empty_like(flux)
        for s in np.unique(season):
            k = season == s; rel[k] = flux[k] / np.median(flux[k][tm[k] < 0.5]) if np.any(tm[k] < 0.5) else np.nan
        good = np.isfinite(rel)
        sc = C.robust_std(rel[good & (tm < 0.5)])
        dep, de_ = C.fit_depth(rel[good], np.full(good.sum(), sc), tm[good])
        res[band] = dict(n=int(good.sum()), n_in=int(np.sum(tm[good] > 0.5)), scatter_frac=float(sc), median_mag=float(np.median(mag)),
                         depth=float(dep), depth_err=float(de_),
                         in_transit_points=[dict(btjd=float(b), phase_min=float(p * C.P * 1440), rel=float(r_), tmpl=float(q)) for b, p, r_, q in
                                            zip(d.btjd.values[ok][good][tm[good] > 0.05], ph[ok][good][tm[good] > 0.05], rel[good][tm[good] > 0.05], tm[good][tm[good] > 0.05])])
    out[name] = res
    print(name, json.dumps(res, indent=1, default=float)[:3000])
json.dump(out, open("out/neowise_fold.json", "w"), indent=1, default=float)
