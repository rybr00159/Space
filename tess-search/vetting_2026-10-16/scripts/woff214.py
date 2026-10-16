"""WCS correction used for TIC 214311906: offsets from fitting the Mira (very red) vs fitting a normal-colour bright star
(Gaia DR3 5925373204131295872, T 11.1, BP-RP 0.48, 41") in the direct image."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from loc5 import Cut, component_images, direct_offset, fit_psf
cfg = json.load(open("loccfg.json"))["214311906"]
for sec in (39, 66, 93, 103, 104):
    C = Cut(214311906, sec, cfg); psf = C.prf(); comps, O = component_images(C, cfg)
    w_mira = direct_offset(C, psf, O)
    r = C.g[C.g.source_id == 5925373204131295872].iloc[0]
    g = C.g; nb = g[(np.hypot(g.x - r.x, g.y - r.y) < 4) & (g.source_id != r.source_id) & (g["T"] < r["T"] + 3)]
    extra = [psf(q.x, q.y) for q in nb.itertuples()]
    F = fit_psf(O, np.sqrt(np.abs(O)) + 1, psf, r.x, r.y, rad=2.0, extra=extra, box=1.0)
    w_ref = np.array([F["x"] - r.x, F["y"] - r.y])
    print(f"S{sec}: offset from the Mira fit ({w_mira[0]:+.2f}, {w_mira[1]:+.2f}) px; from the blue reference star ({w_ref[0]:+.2f}, {w_ref[1]:+.2f}) px; difference {np.hypot(*(w_mira - w_ref)) * 21:.1f}\"")
