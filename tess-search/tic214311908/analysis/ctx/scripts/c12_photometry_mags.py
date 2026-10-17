"""Magnitudes of the M3V, the M-dwarf pair and the Mira in follow-up bands: SkyMapper DR4 (u v g r i z, PSF mags),
2MASS (J H Ks), VHS DR5 (J Ks), APASS9 (B V g r i), ASAS-SN variable-star catalogue entry of the Mira (V mean, amplitude),
Gaia DR3 (G BP RP) and Gaia-based V, R_C, I_C transformations (Riello et al. 2021 / Gaia EDR3 documentation 5.5.1, Table 5.9).
Output: out/mags.json
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

stars = {"M3V": (260.9696737, -51.3238212, 14.281043, 2.7561493), "pairA": (260.9748734, -51.3240552, 14.654037, 2.8778114),
         "Mira": (260.9776003, -51.3237201, 12.660598, 4.6837482)}
cats = {"smss4": "II/379/smssdr4", "2mass": "II/246/out", "vhs5": "II/367/vhs_dr5", "apass9": "II/336/apass9",
        "asassn_var": "II/366/catv2021", "vsx": "B/vsx/vsx", "gaia_lpv": "I/358/vlpv"}
v = Vizier(columns=["**"], row_limit=20)
out = {}
for lab, (ra, dec, G, bprp) in stars.items():
    c = SkyCoord(ra, dec, unit="deg"); out[lab] = {"gaia": {"G": G, "bp_rp": bprp}}
    # Gaia -> Johnson-Cousins (EDR3 doc Table 5.9, valid -0.5 < BP-RP < 5.0 for V, I; 'R' -0.3..4.0)
    x = bprp
    out[lab]["gaia_transform"] = {
        "V": G - (-0.02704 + 0.01424 * x - 0.2156 * x**2 + 0.01426 * x**3),
        "R_C": G - (-0.02275 + 0.3961 * x - 0.1243 * x**2 - 0.01396 * x**3 + 0.003775 * x**4),
        "I_C": G - (0.01753 + 0.76 * x - 0.0991 * x**2)}
    for key, cat in cats.items():
        try:
            r = v.query_region(c, radius=(4 if lab != "Mira" else 6) * u.arcsec, catalog=cat)
        except Exception as e:
            out[lab][key] = f"error {e}"; continue
        if len(r) == 0:
            out[lab][key] = None; continue
        t = r[0]
        # nearest row
        rac = [k for k in t.colnames if k in ("RAJ2000", "RA_ICRS", "RAICRS", "_RAJ2000", "RAdeg", "RA")] or [k for k in t.colnames if k.upper().startswith("RA")]
        dec_ = [k for k in t.colnames if k in ("DEJ2000", "DE_ICRS", "DEICRS", "_DEJ2000", "DEdeg", "DE", "Dec")] or [k for k in t.colnames if k.upper().startswith("DE")]
        try:
            cc = SkyCoord(t[rac[0]], t[dec_[0]], unit="deg")
        except Exception:
            cc = SkyCoord(t[rac[0]], t[dec_[0]], unit=(u.hourangle, u.deg))
        i = int(np.argmin(cc.separation(c).arcsec))
        row = {k: (float(t[k][i]) if np.issubdtype(type(t[k][i]), np.number) and not np.ma.is_masked(t[k][i]) else str(t[k][i])) for k in t.colnames}
        row["_sep_arcsec"] = float(cc.separation(c).arcsec[i]); row["_n_within"] = len(t)
        out[lab][key] = row
json.dump(out, open("out/mags.json", "w"), indent=1, default=str)
for lab in out:
    print(lab, {k: round(v_, 2) for k, v_ in out[lab]["gaia_transform"].items()})
    for key in cats:
        r = out[lab].get(key)
        if isinstance(r, dict):
            keep = {k: r[k] for k in r if any(s in k for s in ("mag", "Jmag", "Hmag", "Kmag", "Ksmag", "Vmag", "Bmag", "gmag", "rmag", "imag", "zmag", "Amp", "Period", "Type", "max", "min", "Name", "_sep", "Per", "Amag", "flags", "Qflg", "Cflg", "ngood"))}
            print("   ", key, json.dumps(keep, default=str)[:600])
        else:
            print("   ", key, r)
