"""R03: the same Sector 39 pixels that TGLC used, analysed directly (independent of ctx and loc code).

1. TESScut S39 21x21 (10-min). Cadences with QUALITY in {0, 2048}. Per cadence a background plane fitted to the
   faintest 30% of pixels is subtracted.
2. Per-pixel transit amplitude image A (e-/s lost) from a joint fit over the 4 S39 transit windows (quadratic baseline
   per +-0.3 d window, trapezoid T14 1.63 h). Times: TESScut TIME (FFI-direction barycentring; target correction
   < 1 min, negligible for 0.6-d windows and a 98-min transit).
3. Noise: the same estimator in single fake windows at independent centres (>= 0.6 d apart, away from transits and
   phase 0.5); the 4-window amplitude noise is emulated by averaging 4 random fake single-window images.
4. Direct fit A = a*PRF(source) + b over pixels within 3.5 px of the M3V; hypotheses M3V (with/without lateral colour
   shift), G17.5 (3.1"), G17.7 (5.5"), pair A, Mira (Gaia position and loc's measured apparent position), free grid.
5. Emulated TGLC: dF_k = sum(w_k * A) for the 6 TGLC stars (weights from leak.py) vs the real TGLC dF.
6. Injection: PRF source with a = observed amplitude on the M3V (or on G17.7 / G17.5) + fake noise images ->
   distribution of chi2(G17.7) - chi2(M3V), chi2(G17.5) - chi2(M3V) and of the free-fit offset; and the same
   injections pushed through the TGLC leak-model GLS (covariance from the fake images' projections).
Output: out/r03_pixels_s39.json, figs/r03_amp_s39.png
"""
import os, sys, json, glob, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
import leak as L
import rlib as R

F = L.Field(39)
with fits.open(F.file) as h:
    d = h[1].data; t = np.array(d["TIME"], float); fl = np.array(d["FLUX"], float); q = d["QUALITY"]
ok = np.isfinite(t) & ((q == 0) | (q == 2048)) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
t, fl = t[ok], fl[ok]
# background plane from the faintest 30% pixels (median image ranking)
med = np.median(fl, axis=0); yy, xx = np.mgrid[0:21, 0:21]
faint = med <= np.percentile(med, 30)
Ab = np.column_stack([np.ones(faint.sum()), xx[faint], yy[faint]])
coef = np.linalg.lstsq(Ab, fl[:, faint].T, rcond=None)[0]           # 3 x ncad
bkg = coef[0][:, None, None] + coef[1][:, None, None] * xx + coef[2][:, None, None] * yy
fl = fl - bkg
ex = R.exptime(39)
tcs, _ = R.centres(t.min(), t.max()); tcs = tcs[(tcs > t.min() + 0.3) & (tcs < t.max() - 0.3)]   # complete windows only (same 4 transits as TGLC)
t5, _ = R.centres(t.min(), t.max(), 0.5); avoid = np.concatenate([tcs, t5])
half = 0.30


def amp_image(centres):
    """Joint amplitude image over windows `centres` (vectorised over pixels)."""
    blocks = []; ys = []
    nb = len(centres); npar = 1 + 3 * nb; rows = []
    for b, tc in enumerate(centres):
        m = np.abs(t - tc) < half
        s = R.shape(t[m], tc, ex); x = (t[m] - tc) / half
        A = np.zeros((m.sum(), npar)); A[:, 0] = -s
        A[:, 1 + 3 * b] = 1; A[:, 2 + 3 * b] = x; A[:, 3 + 3 * b] = x * x
        rows.append(A); ys.append(fl[m].reshape(m.sum(), -1))
    A = np.vstack(rows); Y = np.vstack(ys)
    c = np.linalg.lstsq(A, Y, rcond=None)[0]
    return c[0].reshape(21, 21)


Areal = amp_image(tcs)
# fake single windows
cs = np.arange(t.min() + half, t.max() - half, 0.02)
cs = cs[np.min(np.abs(cs[:, None] - avoid[None, :]), axis=1) > 2 * half]
sel = [cs[0]]
for c in cs[1:]:
    if c - sel[-1] >= 2 * half: sel.append(c)
fakes = []
for c in sel:
    m = np.abs(t - c) < half
    if m.sum() < 60: continue
    fakes.append(amp_image([c]))
fakes = np.array(fakes); nf = len(fakes)
var1 = np.var(fakes, axis=0, ddof=1)
var4 = var1 / len(tcs)
sig = np.sqrt(var4)
print(f"{nf} independent fake single windows; real transits {len(tcs)}")

x3, y3 = F.xy(L.STARS["M3V"])
fitmask = np.hypot(xx - x3, yy - y3) < 3.5
mira_pos = json.load(open("tic214311908/analysis/loc/out/mira_pos.json"))
def fit(img, src, w=None):
    w = 1 / var4 if w is None else w
    P = F.stamp(*src)[fitmask]; y = img[fitmask]; ww = w[fitmask]
    A = np.column_stack([P, np.ones(len(P))])
    Aw = A * ww[:, None]; c = np.linalg.solve(A.T @ Aw, Aw.T @ y)
    r = y - A @ c; return c[0], float(np.sum(r * r * ww))

pa = np.radians(165.0); uEN = np.array([np.sin(pa), np.cos(pa)])
def shifted(xy, arcsec, pa_deg=165.0):
    u = np.array([np.sin(np.radians(pa_deg)), np.cos(np.radians(pa_deg))]); dd = F.Ji @ (u * arcsec); return (xy[0] + dd[0], xy[1] + dd[1])
mo = None
try:
    mp = mira_pos["39"] if "39" in mira_pos else mira_pos.get("S39")
    mo = mp.get("offset_EN_arcsec") or mp.get("offset_arcsec") or None
except Exception:
    pass
hyp = {"M3V": (x3, y3), "M3V_colour3.5": shifted((x3, y3), 3.5), "G17.5_3.1": F.xy(L.EXTRA["G17.5_3.1"]),
       "G17.7_5.5": F.xy(L.EXTRA["G17.7_5.5"]), "pairA": F.xy(L.STARS["pairA"]), "Mira_gaia": F.xy(L.STARS["Mira"]),
       "Mira_shift7.5": shifted(F.xy(L.STARS["Mira"]), 7.5, 173.0)}
res = {"n_fake_windows": nf, "fit_radius_px": 3.5, "n_pix": int(fitmask.sum())}
real = {h: fit(Areal, v) for h, v in hyp.items()}
gE = np.arange(-35, 5.1, 0.5); gN = np.arange(-15, 15.1, 0.5)
GKEYS = [(de, dn) for dn in gN for de in gE]
GP = {k: F.en_to_xy(*k) for k in GKEYS}
PG = np.array([F.stamp(*GP[k])[fitmask] for k in GKEYS])          # ngrid x npix (precomputed once)
def free(img, w=None):
    """Vectorised weighted fit a*PRF + b at every grid point; returns (chi2_min, key)."""
    w = (1 / var4)[fitmask] if w is None else w[fitmask]
    y = img[fitmask]
    S11 = PG ** 2 @ w; S1 = PG @ w; S0 = w.sum(); Sy = PG @ (w * y); Y0 = np.sum(w * y); YY = np.sum(w * y * y)
    det = S11 * S0 - S1 ** 2
    a = (Sy * S0 - S1 * Y0) / det; b = (S11 * Y0 - S1 * Sy) / det
    chi = YY - 2 * a * Sy - 2 * b * Y0 + a * a * S11 + 2 * a * b * S1 + b * b * S0
    i = int(np.argmin(chi)); return float(chi[i]), GKEYS[i]
fb = free(Areal)
m3en = F.xy_to_en(x3, y3)
res["real"] = {h: dict(amp=float(a), chi2=float(c2)) for h, (a, c2) in real.items()}
res["real_free"] = dict(chi2=fb[0], EN_from_Mira=list(fb[1]), offset_from_M3V=float(np.hypot(fb[1][0] - m3en[0], fb[1][1] - m3en[1])))
print("REAL S39 pixel fit (chi2, n_pix=%d):" % fitmask.sum(), {h: (round(v[0], 2), round(v[1], 1)) for h, v in real.items()})
print("   free best", fb, "M3V at", m3en.round(1))
# chi2 calibration of the real fit
chi_fake = [fit(fakes[i] / 1.0, hyp["M3V"], 1 / var1)[1] for i in range(nf)]
res["fake_single_window_chi2_median"] = float(np.median(chi_fake))
print("   chi2 of fake single windows (M3V model, own var):", np.median(chi_fake).round(1), "for", fitmask.sum() - 2, "dof")

# emulated TGLC
W = {k: F.weights(k, "psf") for k in L.NAMES}
Wap = {k: F.weights(k, "ap") for k in L.NAMES}
emu = {k: float(np.sum(W[k] * Areal)) for k in L.NAMES}
emu_ap = {k: float(np.sum(Wap[k] * Areal)) for k in L.NAMES}
proj_f = np.array([[np.sum(W[k] * fk) for k in L.NAMES] for fk in fakes])
emu_err = (np.std(proj_f, axis=0, ddof=1) / np.sqrt(len(tcs))).round(2)
r01 = json.load(open("out/r01_tglc_dF.json"))["S39_psf_flux_h0.3"]["dF"]
r01a = json.load(open("out/r01_tglc_dF.json"))["S39_aperture_flux_h0.3"]["dF"]
res["emulated_tglc_psf_dF"] = emu; res["emulated_tglc_psf_err"] = dict(zip(L.NAMES, emu_err.tolist()))
res["emulated_tglc_aperture_dF"] = emu_ap
res["real_tglc_psf_dF_referee"] = r01; res["real_tglc_aperture_dF_referee"] = r01a
print("emulated TGLC psf dF:", {k: round(v, 2) for k, v in emu.items()}, "+-", emu_err)
print("real TGLC psf dF    :", r01)
print("emulated TGLC ap dF :", {k: round(v, 2) for k, v in emu_ap.items()})
print("real TGLC ap dF     :", r01a)

# leak-model GLS on emulated dF, covariance from fake projections
Cproj = np.cov(proj_f.T) / len(tcs); Ci = np.linalg.inv(Cproj)
def gls(o, m):
    a = m @ Ci @ m; dd = m @ Ci @ o / a; r = o - dd * m; return float(r @ Ci @ r)
def gls_min(o):
    a = np.einsum('ij,jk,ik->i', GVM, Ci, GVM); bb = GVM @ (Ci @ o); return float(o @ Ci @ o - np.max(bb * bb / a))
def cvec(p): img = F.stamp(*p); return np.array([np.sum(W[k] * img) for k in L.NAMES])
HV = {h: cvec(v) for h, v in hyp.items()}
GV = {k: cvec(GP[k]) for k in GKEYS if (k[0] % 1 == 0 and k[1] % 1 == 0)}
GVK = list(GV); GVM = np.array([GV[k] for k in GVK])
o_emu = np.array([emu[k] for k in L.NAMES])
lk = {h: gls(o_emu, m) for h, m in HV.items()}; kb = min(GV, key=lambda k: gls(o_emu, GV[k]))
res["leak_model_on_emulated"] = dict(chi2=lk, free_EN=list(kb), free_chi2=gls(o_emu, GV[kb]))
print("leak-model GLS on emulated dF:", {h: round(v, 2) for h, v in lk.items()}, "free", kb, round(gls(o_emu, GV[kb]), 2))

# injections
rng = np.random.default_rng(42)
NI = 400
inj = {}
for truth in ("M3V", "M3V_colour3.5", "G17.7_5.5", "G17.5_3.1"):
    a0 = real[truth][0]
    d77, d75, offs, lkd, lkoff = [], [], [], [], []
    for i in range(NI):
        idx = rng.choice(nf, len(tcs), replace=False)
        noise = fakes[idx].mean(axis=0)
        img = noise + a0 * F.stamp(*hyp[truth])
        c = {h: fit(img, hyp[h])[1] for h in ("M3V", "G17.7_5.5", "G17.5_3.1")}
        d77.append(c["G17.7_5.5"] - c["M3V"]); d75.append(c["G17.5_3.1"] - c["M3V"])
        if i < 150:
            fbi = free(img); offs.append(float(np.hypot(fbi[1][0] - m3en[0], fbi[1][1] - m3en[1])))
            o = np.array([np.sum(W[k] * img) for k in L.NAMES])
            lkd.append(gls(o, HV["M3V"]) - gls_min(o))
    inj[truth] = dict(dchi2_G17_7_minus_M3V_pct=np.percentile(d77, [2.5, 16, 50, 84, 97.5]).round(1).tolist(),
                      dchi2_G17_5_minus_M3V_pct=np.percentile(d75, [2.5, 16, 50, 84, 97.5]).round(1).tolist(),
                      free_offset_from_M3V_pct=np.percentile(offs, [50, 68, 95]).round(1).tolist(),
                      leak_model_dchi2_M3V_minus_free_pct=np.percentile(lkd, [50, 84, 95, 99]).round(1).tolist(),
                      frac_leak_dchi2_ge_13=float(np.mean(np.array(lkd) >= 13)),
                      frac_dchi2_G17_7_ge_observed=float(np.mean(np.array(d77) >= real["G17.7_5.5"][1] - real["M3V"][1])),
                      frac_dchi2_G17_5_ge_observed=float(np.mean(np.array(d75) >= real["G17.5_3.1"][1] - real["M3V"][1])),
                      frac_leak_dchi2_ge_emulated_observed=float(np.mean(np.array(lkd) >= lk["M3V"] - gls_min(o_emu))))
    print("INJ", truth, inj[truth])
res["injections"] = inj
json.dump(res, open("out/r03_pixels_s39.json", "w"), indent=1, default=float)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a_, im, ttl in ((ax[0], Areal, "S39 amplitude image (e-/s lost)"), (ax[1], Areal / sig, "amplitude / sigma"),
                    (ax[2], (Areal - real["M3V"][0] * F.stamp(*hyp["M3V"])) / sig, "residual after M3V PRF fit / sigma")):
    p = a_.imshow(im, origin="lower", cmap="RdBu_r", vmin=-np.nanmax(np.abs(im[fitmask])), vmax=np.nanmax(np.abs(im[fitmask])))
    for h, (x, y) in hyp.items():
        if h in ("M3V", "G17.7_5.5", "G17.5_3.1", "pairA", "Mira_gaia"): a_.plot(x, y, "k+" if h != "M3V" else "g*", ms=8)
    a_.set_xlim(5.5, 15.5); a_.set_ylim(5.5, 15.5); a_.set_title(ttl, fontsize=9); plt.colorbar(p, ax=a_)
plt.tight_layout(); plt.savefig("figs/r03_amp_s39.png", dpi=90)
