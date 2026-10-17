"""R02: how robust is ctx's claim that TGLC S39 psf_flux prefers a source ~7" N of the M3V (dchi2 ~ 13)?

1. Observed dF vectors: (a) ctx's values (tglc_leak.json), (b) referee's local-baseline values (r01).
2. Covariances: (A) single-window fake depths at independent centres (>= 0.6 d apart) / n_transits,
   (B) ctx-style: 1-d biweight detrend + template amplitude at 400 random phases (re-implemented here),
   (C) diagonal of (B), (D) Ledoit-Wolf shrinkage of (B).
3. Hypotheses: M3V, pairA, Mira, G17.5 (3.1"), G17.7 (5.5"), free position grid (1" steps).
4. Leave-one-transit-out (3 transits in S39).
5. Model-mismatch test: truth = source on the M3V computed with a perturbed PSF/geometry (lateral-colour shift of the
   source image, PSF broadened by a Gaussian, different PRF reference position, per-star TGLC position errors),
   noise-free; fitted with the baseline model and covariance (A)/(B): what dchi2(M3V - best) does mismatch alone make?
Output: out/r02_leak_gls.json
"""
import os, sys, json, glob, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.io import fits
from wotan import flatten
import leak as L
import rlib as R

OUT = {}
F = L.Field(39)
names = L.NAMES
Wpsf = {k: F.weights(k, "psf") for k in names}


def cvec(src, Ws=Wpsf):
    img = F.stamp(*src); return np.array([np.sum(Ws[k] * img) for k in names])


def gls(o, Ci, m):
    a = m @ Ci @ m; d = (m @ Ci @ o) / a; r = o - d * m; return d, 1 / np.sqrt(a), float(r @ Ci @ r)


hyp = {"M3V": F.xy(L.STARS["M3V"]), "pairA": F.xy(L.STARS["pairA"]), "Mira": F.xy(L.STARS["Mira"]),
       "G17.5_3.1": F.xy(L.EXTRA["G17.5_3.1"]), "G17.7_5.5": F.xy(L.EXTRA["G17.7_5.5"])}
m3en = F.xy_to_en(*hyp["M3V"])
gE = np.arange(-30, 0.1, 1.0); gN = np.arange(-12, 18.1, 1.0)
GRID = {(de, dn): cvec(F.en_to_xy(de, dn)) for dn in gN for de in gE}
HYPV = {h: cvec(v) for h, v in hyp.items()}


def evaluate(o, C, tag):
    Ci = np.linalg.inv(C)
    res = {h: dict(zip(("dF", "err", "chi2"), map(float, gls(o, Ci, m)))) for h, m in HYPV.items()}
    c2 = {k: gls(o, Ci, m)[2] for k, m in GRID.items()}; kb = min(c2, key=c2.get)
    res["free"] = dict(chi2=float(c2[kb]), EN_from_Mira=[float(kb[0]), float(kb[1])],
                       offset_from_M3V_arcsec=float(np.hypot(kb[0] - m3en[0], kb[1] - m3en[1])))
    res["dchi2_M3V_minus_free"] = res["M3V"]["chi2"] - res["free"]["chi2"]
    res["dchi2_M3V_minus_G17.7"] = res["M3V"]["chi2"] - res["G17.7_5.5"]["chi2"]
    print(f"{tag:40s} chi2: M3V {res['M3V']['chi2']:6.2f} G17.5 {res['G17.5_3.1']['chi2']:6.2f} G17.7 {res['G17.7_5.5']['chi2']:6.2f} "
          f"pairA {res['pairA']['chi2']:6.2f} Mira {res['Mira']['chi2']:6.2f} free {res['free']['chi2']:6.2f} at {res['free']['EN_from_Mira']}  "
          f"d(M3V-free) {res['dchi2_M3V_minus_free']:5.2f}", flush=True)
    return res


# ---------------- data
def load(sid, col="psf_flux"):
    f = glob.glob(f"{L.TGLC}/hlsp_tglc_tess_ffi_gaiaid-{sid}-s0039-*.fits")[0]; d = fits.getdata(f, 1)
    ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0) & np.isfinite(d[col])
    return np.array(d["time"][ok], float), np.array(d[col][ok], float)

T = {k: load(sid) for k, sid in L.STARS.items()}
t0 = T["M3V"][0]; ex = R.exptime(39)
tcs, _ = R.centres(t0.min(), t0.max()); tcs = tcs[(tcs > t0.min()) & (tcs < t0.max())]
t5, _ = R.centres(t0.min(), t0.max(), 0.5); avoid = np.concatenate([tcs, t5])
print("S39 transit centres", tcs)

ctx = json.load(open("tic214311908/analysis/ctx/out/tglc_leak.json"))["psf_flux"]["S39"]
o_ctx = np.array([ctx["observed_dF"][k] for k in names])
r01 = json.load(open("out/r01_tglc_dF.json"))["S39_psf_flux_h0.3"]
o_ref = np.array([r01["dF"][k] for k in names])

# covariance A: single-window fakes at independent centres
half = 0.30
cs = np.arange(t0.min() + half, t0.max() - half, 0.02)
cs = cs[np.min(np.abs(cs[:, None] - avoid[None, :]), axis=1) > 2 * half]
sel = [cs[0]]
for c in cs[1:]:
    if c - sel[-1] >= 2 * half: sel.append(c)
FW = []
for c in sel:
    v = [R.joint_fit(*T[k], [c], ex, half)[0] for k in names]
    if np.all(np.isfinite(v)): FW.append(v)
FW = np.array(FW); CA = np.cov(FW.T) / len(tcs)
print("cov A from", len(FW), "independent single windows")

# covariance B: ctx-style random phases on a 1-d biweight-detrended light curve (re-implemented)
def detr(t, f, mask):
    off = 10 * np.nanmax(np.abs(f)) + 1
    tr = flatten(t[~mask], (f[~mask] + off) / off, method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    k = np.isfinite(tr); return f - (np.interp(t, t[~mask][k], tr[k]) * off - off)

def tmpl(t, ph0):
    ph = ((t - R.T0) / R.P - ph0 + 0.5) % 1 - 0.5
    return R.shape(ph * R.P, 0.0, ex)

RES = {}
for k, (t, f) in T.items():
    ph = ((t - R.T0) / R.P + 0.5) % 1 - 0.5; ph5 = ((t - R.T0) / R.P) % 1 - 0.5
    mask = (np.abs(ph * R.P) < R.T14) | (np.abs(ph5 * R.P) < R.T14)
    r = detr(t, f, mask); s = 1.4826 * np.median(np.abs(r[~mask] - np.median(r[~mask])))
    good = ~((~mask) & (np.abs(r) > 5 * s)); RES[k] = (t[good], r[good])
rng = np.random.default_rng(3); phs = []
while len(phs) < 400:
    p = rng.uniform(0, 1)
    if min(abs(p), abs(p - 0.5), abs(1 - p)) > 2.5 * R.T14 / R.P: phs.append(p)
RB = np.array([[-np.sum((q := tmpl(RES[k][0], p)) * RES[k][1]) / np.sum(q * q) for p in phs] for k in names])
o_B = np.array([-np.sum((q := tmpl(RES[k][0], 0.0)) * RES[k][1]) / np.sum(q * q) for k in names])
CB = np.cov(RB)
# effective number of independent phases: template support is T14 wide
n_eff_B = int(round(R.P * (1 - 4 * 2.5 * R.T14 / R.P) / R.T14))
lw_shr = None
try:
    from sklearn.covariance import LedoitWolf
    CD = LedoitWolf().fit(RB.T).covariance_; lw_shr = "sklearn"
except Exception:
    # simple Ledoit-Wolf-type shrinkage towards the diagonal with intensity 0.2
    CD = 0.8 * CB + 0.2 * np.diag(np.diag(CB)); lw_shr = "fixed 0.2 towards diagonal"
CC = np.diag(np.diag(CB))
OUT["inputs"] = dict(o_ctx=o_ctx.tolist(), o_ref_localbaseline=o_ref.tolist(), o_ctxstyle_reimpl=o_B.round(2).tolist(),
                     errA=np.sqrt(np.diag(CA)).round(2).tolist(), errB=np.sqrt(np.diag(CB)).round(2).tolist(),
                     corrA=np.round(np.corrcoef(FW.T), 3).tolist(), corrB=np.round(np.corrcoef(RB), 3).tolist(),
                     eigB=np.linalg.eigvalsh(CB).round(4).tolist(), eigA=np.linalg.eigvalsh(CA).round(4).tolist(),
                     n_single_windows_A=len(FW), n_random_phases_B=400, n_independent_phases_B_approx=n_eff_B, shrinkage=lw_shr)
print("o_ctx", o_ctx, "\no_ref", o_ref, "\no_B(reimpl)", o_B.round(2))
print("errA", np.sqrt(np.diag(CA)).round(2), "errB", np.sqrt(np.diag(CB)).round(2))
print("eig B", np.linalg.eigvalsh(CB).round(4), " eig A", np.linalg.eigvalsh(CA).round(4))

E = {}
E["ctx_dF__covB"] = evaluate(o_ctx, CB, "ctx dF, ctx-style cov (B)")
E["reimplB_dF__covB"] = evaluate(o_B, CB, "reimpl ctx-style dF, cov B")
E["ref_dF__covA"] = evaluate(o_ref, CA, "referee dF, single-window cov (A)")
E["ref_dF__covB"] = evaluate(o_ref, CB, "referee dF, cov B")
E["ctx_dF__diag"] = evaluate(o_ctx, CC, "ctx dF, diagonal cov")
E["ref_dF__diag"] = evaluate(o_ref, CC, "referee dF, diagonal cov")
E["ctx_dF__shrunk"] = evaluate(o_ctx, CD, f"ctx dF, shrunk cov ({lw_shr})")
# drop pairB (TGLC pair B fit is nearly a duplicate of pair A: corr ~1.00)
keep = [0, 1, 3, 4, 5]
def evaluate_sub(o, C, tag, idx):
    global HYPV, GRID
    H0, G0 = HYPV, GRID
    HYPV = {h: m[idx] for h, m in H0.items()}; GRID = {k: m[idx] for k, m in G0.items()}
    try:
        return evaluate(o[idx], C[np.ix_(idx, idx)], tag)
    finally:
        HYPV, GRID = H0, G0
E["ctx_dF__covB__no_pairB"] = evaluate_sub(o_ctx, CB, "ctx dF, cov B, without pairB", keep)
E["ctx_dF__covB__no_Mira_no_pairB"] = evaluate_sub(o_ctx, CB, "ctx dF, cov B, without Mira & pairB", [0, 1, 4, 5])
OUT["evaluations"] = E

# ---------------- leave-one-transit-out (referee estimator, cov A scaled to 2 transits)
LOO = {}
for i in range(len(tcs)):
    use = np.delete(tcs, i)
    o = np.array([R.joint_fit(*T[k], use, ex, half)[0] for k in names])
    LOO[f"drop_{tcs[i]:.2f}"] = dict(dF=o.round(2).tolist(), **{kk: vv for kk, vv in evaluate(o, CA * len(tcs) / len(use), f"LOO drop {tcs[i]:.2f}").items() if kk in ("M3V", "G17.7_5.5", "free", "dchi2_M3V_minus_free")})
for i in range(len(tcs)):
    o = np.array([R.joint_fit(*T[k], [tcs[i]], ex, half)[0] for k in names])
    LOO[f"only_{tcs[i]:.2f}"] = dict(dF=o.round(2).tolist(), **{kk: vv for kk, vv in evaluate(o, CA * len(tcs), f"single transit {tcs[i]:.2f}").items() if kk in ("M3V", "G17.7_5.5", "free", "dchi2_M3V_minus_free")})
OUT["leave_one_out"] = LOO

# ---------------- model mismatch (noise-free truth on the M3V with a perturbed model)
dF0 = 13.9
mm = {}
x3, y3 = hyp["M3V"]
pa = np.radians(165.0); uEN = np.array([np.sin(pa), np.cos(pa)])   # S39: away from boresight (ctx's value)
def truth_shift(arcsec):
    d = F.Ji @ (uEN * arcsec); return dF0 * cvec((x3 + d[0], y3 + d[1]))
def truth_blur(sig):
    img = gaussian_filter(F.stamp(x3, y3), sig); return dF0 * np.array([np.sum(Wpsf[k] * img) for k in names])
def truth_prfpos(dcol, drow):
    import PRF
    p2 = PRF.TESS_PRF(F.cam, F.ccd, 39, F.col0 + F.xm + dcol, F.row0 + F.ym + drow, localdatadir=L.PRFDIR)
    pad = 8; img = p2.locate(x3 + pad, y3 + pad, (37, 37))[pad:pad + 21, pad:pad + 21]
    W2 = {}
    for k in names:   # TGLC's own ePSF also from that PRF (field-averaged ePSF)
        c = F.cut[k]; x, y = F.xy(L.STARS[k]); full = p2.locate(x + pad, y + pad, (37, 37))[pad:pad + 21, pad:pad + 21]
        cx, cy = int(np.round(x)), int(np.round(y)); s11 = p2.locate(x - cx + 5, y - cy + 5, (11, 11)).sum()
        PX = full[c["oy"]:c["oy"] + 5, c["ox"]:c["ox"] + 5] / s11; yy, xx = np.mgrid[0:5, 0:5] - 2.0
        A = np.column_stack([PX.ravel(), xx.ravel(), yy.ravel(), np.ones(25)]); w = np.zeros((21, 21))
        w[c["oy"]:c["oy"] + 5, c["ox"]:c["ox"] + 5] = np.linalg.pinv(A)[0].reshape(5, 5); W2[k] = w
    # truth: real TGLC uses W2 (its ePSF = this PRF) and the source image = this PRF; model uses baseline PRF
    return dF0 * np.array([np.sum(W2[k] * img) for k in names])
def truth_starjitter(seed, sig_px=0.05):
    rr = np.random.default_rng(seed); W2 = {k: F.weights(k, "psf", star_shift=tuple(rr.normal(0, sig_px, 2))) for k in names}
    img = F.stamp(x3, y3); return dF0 * np.array([np.sum(W2[k] * img) for k in names])

cases = {f"source_colour_shift_{s}arcsec_PA165": truth_shift(s) for s in (-4.0, -2.0, 2.0, 4.0)}
cases.update({f"psf_blur_gauss_{s}px": truth_blur(s) for s in (0.2, 0.35, 0.5)})
cases.update({f"prf_from_ccd_offset_{dc}_{dr}px": truth_prfpos(dc, dr) for dc, dr in ((-75, -75), (75, 75), (-75, 75), (75, -75), (200, 0))})
cases.update({f"per_star_tglc_position_err_0.05px_seed{s}": truth_starjitter(s) for s in range(6)})
cases.update({f"per_star_tglc_position_err_0.10px_seed{s}": truth_starjitter(100 + s, 0.10) for s in range(6)})
for cname, o in cases.items():
    mm[cname] = {}
    for cn, C in (("covB", CB), ("covA", CA)):
        e = evaluate(o, C, f"MISMATCH {cname} [{cn}]")
        mm[cname][cn] = dict(chi2_M3V=e["M3V"]["chi2"], chi2_free=e["free"]["chi2"], free_EN=e["free"]["EN_from_Mira"],
                             dchi2=e["dchi2_M3V_minus_free"], chi2_G17_7=e["G17.7_5.5"]["chi2"], chi2_pairA=e["pairA"]["chi2"])
OUT["model_mismatch_noise_free"] = mm
json.dump(OUT, open("out/r02_leak_gls.json", "w"), indent=1, default=float)
