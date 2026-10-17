"""R02b: TGLC-faithful geometry. Each star's ePSF is placed where TGLC placed it (STAR_X/STAR_Y on the detector grid,
common WCS offset removed), the true source image at its proper-motion-propagated, calibrated position, optionally
displaced by TESS lateral colour (radially away from the S39 boresight, PA 165 deg as used by ctx; loc measured the
Mira's S39 displacement at PA ~173 deg, 7-8").
Compares with ctx's geometry and refits ctx's observed dF with covariance B (ctx-style, re-implemented in r02).
Output: out/r02b_tglc_geometry.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import leak as L
import r02_leak_gls as A      # reuses its data products (runs r02 quietly first)

F = A.F; names = L.NAMES
tp, d = L.tglc_positions(F)
out = {"common_offset_px": d.tolist()}
for k in names:
    x, y = F.xy(L.STARS[k]); out.setdefault("tglc_minus_true_arcsec", {})[k] = (F.J @ (np.array(tp[k]) - np.array([x, y]))).round(2).tolist()
print("TGLC ePSF position minus true (PM-propagated) position, arcsec E,N:", out["tglc_minus_true_arcsec"])
Wt = {k: L.weights_at(F, k, tp[k]) for k in names}
pa = np.radians(165.0); uEN = np.array([np.sin(pa), np.cos(pa)])
def src(xy, shift_arcsec=0.0):
    dd = F.Ji @ (uEN * shift_arcsec); return (xy[0] + dd[0], xy[1] + dd[1])
def cv(xy, W): img = F.stamp(*xy); return np.array([np.sum(W[k] * img) for k in names])
base_m3 = cv(F.xy(L.STARS["M3V"]), A.Wpsf)
res = {}
for sh in (0.0, 2.0, 3.5, 4.0):
    m3 = cv(src(F.xy(L.STARS["M3V"]), sh), Wt)
    res[f"M3V_colour{sh}"] = dict(coeffs=m3.round(4).tolist())
out["coeffs_baseline_M3V"] = base_m3.round(4).tolist(); out["coeffs_tglc_geometry"] = res
print("baseline M3V coeffs", base_m3.round(4))
for kk, v in res.items(): print(kk, v["coeffs"])
# noise-free: truth = TGLC-faithful geometry (no colour shift / colour 3.5"), fitted with the baseline model
Ci = np.linalg.inv(A.CB)
for kk in ("M3V_colour0.0", "M3V_colour3.5"):
    o = 13.9 * np.array(res[kk]["coeffs"])
    e = A.evaluate(o, A.CB, f"noise-free truth {kk}, baseline model, cov B")
    out[f"noisefree_{kk}_baseline_fit"] = dict(chi2_M3V=e["M3V"]["chi2"], chi2_free=e["free"]["chi2"], free_EN=e["free"]["EN_from_Mira"], chi2_G17_7=e["G17.7_5.5"]["chi2"])
# refit the observed ctx dF with the TGLC-faithful model (grid + hypotheses)
H0, G0 = A.HYPV, A.GRID
for sh in (0.0, 3.5):
    A.HYPV = {h: cv(src(v, sh if h in ("M3V", "pairA") else 0.0), Wt) for h, v in A.hyp.items()}
    A.GRID = {k: cv(src(F.en_to_xy(*k), sh), Wt) for k in G0}
    e1 = A.evaluate(A.o_ctx, A.CB, f"ctx dF, TGLC geometry, colour {sh}\" [covB]")
    e2 = A.evaluate(A.o_ref, A.CA, f"ref dF, TGLC geometry, colour {sh}\" [covA]")
    e3 = A.evaluate(A.o_ctx, np.diag(np.diag(A.CB)), f"ctx dF, TGLC geometry, colour {sh}\" [diag]")
    out[f"refit_tglc_geometry_colour{sh}"] = dict(ctx_covB={h: e1[h] for h in ("M3V", "G17.5_3.1", "G17.7_5.5", "pairA", "Mira", "free")},
                                                 ref_covA={h: e2[h] for h in ("M3V", "G17.5_3.1", "G17.7_5.5", "pairA", "Mira", "free")},
                                                 ctx_diag={h: e3[h] for h in ("M3V", "G17.7_5.5", "free")})
A.HYPV, A.GRID = H0, G0
json.dump(out, open("out/r02b_tglc_geometry.json", "w"), indent=1, default=float)
