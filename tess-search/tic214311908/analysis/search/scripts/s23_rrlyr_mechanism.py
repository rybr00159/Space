"""s23: how does the RRab (Gaia DR3 5925372615693467904, 200" S of the M3V) enter the M3V light curve?
Per sector: (1) is the RRab's 3x3 pixel neighbourhood inside the faint-pixel set (lowest 15 % of the sector-median
image) used by phot p05 for its per-cadence background plane? (2) significance of the 0.648528-d signal (3 harmonics on
top of a 1.5-d B-spline) in each per-cadence regressor of the phot pipeline: the faint-pixel background B(t) and the
field-star pointing offsets dx(t), dy(t) (phot/data/jitter_s*.npz). Significance = first-harmonic amplitude relative to
the same fit at 100 fake periods (0.55-0.80 d). A pointing regressor that carries the RRab's variability imprints it,
through the Mira's large pointing response, on the pixels around the Mira and hence on the M3V fit.
Output: out/rrlyr_mechanism.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/search/scripts")
sys.path.insert(1, "tic214311908/analysis/phot/scripts")
from slib import *
from s05_pix_periodic import bspline_basis
from s07_periodic_model import harm_X
import photlib
RR = 5925372615693467904
RNG = np.random.default_rng(23)
pm = json.load(open(f"{ROOT}/out/periodic_model.json"))


def amp(tt, v, Bs, P):
    X = np.column_stack([Bs, harm_X(tt, P)[:, 1:]])
    c, *_ = np.linalg.lstsq(X, v, rcond=None); h = c[Bs.shape[1]:]
    return float(np.hypot(h[0], h[1])), float(np.arctan2(h[1], h[0]))


out = {}
for s in SECTORS:
    C = photlib.Cutout(s)
    g = C.gaia(19.0); W = json.load(open(f"{PHOT}/out/wcs_s{s}.json"))
    r = g[g.source_id == RR].iloc[0]; x, y = r.x + W["global_shift"][0], r.y + W["global_shift"][1]
    q0 = C.fin & (C.q == 0)
    med = np.nanmedian(C.flux[q0], 0); faint = med < np.percentile(med, 15)
    ix, iy = int(round(x)), int(round(y))
    nb = [(j, i) for j in range(iy - 1, iy + 2) for i in range(ix - 1, ix + 2) if 0 <= i < C.nx and 0 <= j < C.ny]
    J = np.load(f"{PHOT}/data/jitter_s{s}.npz")
    res = dict(rr_xy=[float(x), float(y)], n_faint_in_3x3=int(sum(faint[j, i] for j, i in nb)), n_3x3=len(nb),
               M3V_B_leak_amp1_e_s=pm["per_sector"][str(s)]["amp_e_s"][0])
    for key in ("bkg", "dx", "dy"):
        v = J[key]; k = q0 & np.isfinite(v) & np.isfinite(C.t)
        tt, vv = C.t[k], v[k]
        # 30-min bins
        kb = np.floor((tt - tt[0]) * 48).astype(int); u, inv = np.unique(kb, return_inverse=True); n = np.bincount(inv)
        tb, vb = np.bincount(inv, tt) / n, np.bincount(inv, vv) / n
        Bs = bspline_basis(tb)
        a, ph = amp(tb, vb, Bs, 0.648528)
        fk = np.array([amp(tb, vb, Bs, p)[0] for p in RNG.uniform(0.55, 0.80, 100) if abs(p - 0.648528) > 0.01])
        res[key] = dict(amp1=a, phase=ph, fake_median=float(np.median(fk)), p_emp=float(np.mean(fk >= a)), ratio=float(a / np.median(fk)))
    out[s] = res
    print(s, f"faint 3x3 {res['n_faint_in_3x3']}/{res['n_3x3']}; " + "; ".join(f"{k} amp {res[k]['amp1']:.4g} (x{res[k]['ratio']:.1f} fake median, p {res[k]['p_emp']:.2f})" for k in ("bkg", "dx", "dy")), flush=True)
save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/rrlyr_mechanism.json")
