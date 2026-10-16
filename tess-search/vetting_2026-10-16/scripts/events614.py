"""Single-event search for TIC 61141412: box dips 0.5-4 h in each product (detrended, unmasked), with quality context."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "61141412"]
import numpy as np
from astropy.io import fits
from wotan import flatten
import stage2
from lcan import read, bin10
p = stage2.products(61141412); allev = []
for _, r in p.iterrows():
    path = stage2.fetch(r.uri)
    t, q, cols, bkg, aper, crowd = read(path, r.prov)
    c = "PDCSAP_FLUX" if "PDCSAP_FLUX" in cols else "SAP_FLUX"
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[c]) & (cols[c] > 0)
    tt, ff = bin10(t[ok], cols[c][ok] / np.nanmedian(cols[c][ok]), r.exptime)
    fl = flatten(tt, ff, method="biweight", window_length=0.75, break_tolerance=0.4)
    m = np.isfinite(fl); tt, fl = tt[m], fl[m]
    s = 1.4826 * np.median(np.abs(fl - np.median(fl)))
    dt = np.median(np.diff(tt)); found = []
    for dur_h in (0.75, 1.5, 3.0):
        w = max(1, int(round(dur_h / 24 / dt)))
        cs = np.convolve(fl - 1, np.ones(w) / w, mode="same")
        snr = -cs / (s / np.sqrt(w))
        for i in np.argsort(snr)[::-1][:20]:
            if snr[i] < 5: break
            if any(abs(tt[i] - f[0]) < 0.2 for f in found): continue
            # flag cadences within 0.3 d that had nonzero quality (momentum dumps, scattered light) in the raw file
            near = (np.abs(t - tt[i]) < 0.15) & (q != 0)
            found.append((tt[i], dur_h, -cs[i] * 1e6, snr[i], int(near.sum())))
    for f in sorted(found):
        allev.append(dict(sector=int(r.sector), prov=r.prov, col=c, t=round(f[0], 3), dur_h=f[1], depth=round(f[2]), snr=round(f[3], 1), n_flag=f[4]))
        print(f"S{int(r.sector):3d} {r.prov:9s} {c:11s} t={f[0]:.3f} dur {f[1]:.2f} h depth {f[2]:5.0f} ppm S/N {f[3]:4.1f} flagged-cadences-nearby {f[4]}")
json.dump(allev, open("out/61141412/events_single.json", "w"), indent=1)
