"""Target-specific barycentric correction for the TESScut cutouts: TESScut TIME is barycentred for the FFI reference
direction; we compute the Roemer delay towards TIC 214311908 from JPL Horizons positions of TESS and the difference
to TIMECORR. Validation: compare TIMECORR itself with our Roemer delay toward the cutout centre (should agree to the
offset between the FFI reference direction and the cutout, i.e. small). Output: out/timecorr.json."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
out = {}
for s in SECTORS:
    C = Cutout(s); m = np.isfinite(C.t) & np.isfinite(C.timecorr)
    add = roemer_target(s, C.t[m], C.timecorr[m])
    add_c = roemer_target(s, C.t[m], C.timecorr[m], ra=260.977650, dec=-51.323690)
    out[s] = dict(add_s_mean=float(add.mean()), add_s_min=float(add.min()), add_s_max=float(add.max()), timecorr_min_mean=float(C.timecorr[m].mean() * 1440),
                  add_s_mean_cutout_centre=float(add_c.mean()))
    print(f"S{s}: TIMECORR mean {C.timecorr[m].mean()*1440:.3f} min; add for M3V {add.mean():+.2f} s (range {add.min():+.2f}..{add.max():+.2f}); for cutout centre {add_c.mean():+.2f} s")
save_json(out, f"{ROOT}/out/timecorr.json")
