"""TIC 214311906: does the absolute eclipse amplitude (e-/s) stay constant while the Mira's flux changes?"""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from wotan import flatten
from loc5 import Cut, circle_mask, template
cfg = json.load(open("loccfg.json"))["214311906"]
rows = []
for sec in (12, 39, 66, 93, 103, 104):
    C = Cut(214311906, sec, cfg); g = C.good; t = C.t[g]; cube = C.flux[g]
    med = np.median(cube, 0); bpx = med < np.percentile(med, 25); bkg = np.median(cube[:, bpx], 1)
    ap = circle_mask((C.ny, C.nx), C.xt, C.yt, 2.5)
    f = cube[:, ap].sum(1) - bkg * ap.sum()
    F = np.median(f)
    tr = flatten(t, f, method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    ok = np.isfinite(tr); r_ = f[ok] - tr[ok]; y = template(t[ok], cfg)
    A = np.column_stack([np.ones(ok.sum()), y]); c, *_ = np.linalg.lstsq(A, r_, rcond=None)
    s = 1.4826 * np.median(np.abs(r_ - A @ c)); cov = np.linalg.pinv(A.T @ A) * s ** 2
    dF, e = -c[1], np.sqrt(cov[1, 1])
    rows.append((sec, C.camera, C.ccd, F, dF, e))
    print(f"S{sec:3d} cam{C.camera} ccd{C.ccd}: aperture flux {F:9.0f} e/s; eclipse amplitude {dF:6.1f} ± {e:5.1f} e/s; fractional {dF/F*1e6:6.0f} ± {e/F*1e6:4.0f} ppm")
a = np.array([r[3:] for r in rows]); F, dF, e = a[:, 0], a[:, 1], a[:, 2]
w = 1 / e ** 2; mu = np.sum(w * dF) / w.sum(); chi_const = np.sum(((dF - mu) / e) ** 2)
frac = dF / F; ef = e / F; wf = 1 / ef ** 2; muf = np.sum(wf * frac) / wf.sum(); chi_frac = np.sum(((frac - muf) / ef) ** 2)
print(f"Flux range {F.min():.0f}-{F.max():.0f} e/s (ratio {F.max()/F.min():.2f}).")
print(f"Constant absolute amplitude (eclipse on a constant star): chi2 {chi_const:.1f} for {len(F)-1} dof")
print(f"Constant fractional depth (eclipse scales with the Mira's flux): chi2 {chi_frac:.1f} for {len(F)-1} dof")
json.dump(dict(rows=rows, chi_const=chi_const, chi_frac=chi_frac), open("out/214311906/mira_test.json", "w"), indent=1, default=float)
