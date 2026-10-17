"""R07 (test ctx did not run): rotation / spot modulation of the M3V and synchronisation with the 5.69-d period.
TARS (rotation-survey HLSP, 1-px aperture on the M3V) and TGLC S39 psf_flux of the M3V, pair A and the Mira
(TGLC removes the static neighbour model, but a variable neighbour leaks), Lomb-Scargle 0.3-13 d with the transits
masked. A 5.69-d peak in the M3V but not in the Mira/pair would hint at a synchronised close binary.
Output: out/r07_rotation.json
"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
import rlib as R
LC = "tic214311908/analysis/ctx/lc"
out = {}
def ls(t, f, tag):
    ph = ((t - R.T0) / R.P + 0.5) % 1 - 0.5; m = np.abs(ph * R.P) > 0.1
    t, f = t[m], f[m]; f = f / np.nanmedian(f) - 1
    # remove a linear trend per sector-length chunk (Mira's 294-d variability leaks as slow trends)
    p = np.polyfit(t - t.mean(), f, 2); f = f - np.polyval(p, t - t.mean())
    fr = np.linspace(1 / 13.0, 1 / 0.3, 20000); pw = LombScargle(t, f).power(fr)
    i = np.argsort(pw)[::-1]; peaks = []
    for j in i:
        if all(abs(1 / fr[j] - q) / q > 0.1 for q, _ in peaks): peaks.append((1 / fr[j], pw[j]))
        if len(peaks) == 4: break
    k = np.argmin(np.abs(1 / fr - R.P)); amp = np.sqrt(np.mean(f ** 2))
    out[tag] = dict(top_peaks=[(round(a, 3), round(b, 4)) for a, b in peaks], power_at_P=float(pw[k]),
                    fap_top=float(LombScargle(t, f).false_alarm_probability(peaks[0][1], minimum_frequency=fr[0], maximum_frequency=fr[-1])),
                    rms_ppt=float(amp * 1e3), n=len(t))
    print(tag, out[tag])
for f in sorted(glob.glob(f"{LC}/hlsp_tars_*214311908*")):
    d = fits.getdata(f, 1); sec = int(f.split("_s")[1][:4]); ok = np.isfinite(d["FLUX"])
    ls(np.array(d["TIME"][ok], float), np.array(d["FLUX"][ok], float), f"TARS_M3V_S{sec}")
for name, sid in (("M3V", 5925373272850784768), ("pairA", 5925373272850782208), ("Mira", 5925373268528684288), ("ref1", 5925373204131295872)):
    f = glob.glob(f"{LC}/tglc/hlsp_tglc_tess_ffi_gaiaid-{sid}-s0039-*.fits")[0]; d = fits.getdata(f, 1)
    ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0) & np.isfinite(d["psf_flux"])
    t = np.array(d["time"][ok], float); fl = np.array(d["psf_flux"][ok], float)
    fl = fl - np.nanmedian(fl) + np.nanmedian(fl) if name != "Mira" else fl
    ls(t, fl, f"TGLC_S39_{name}")
json.dump(out, open("out/r07_rotation.json", "w"), indent=1, default=float)
