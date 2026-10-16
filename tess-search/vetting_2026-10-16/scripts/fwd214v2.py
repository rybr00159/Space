"""Forward model of the aperture test: for each hypothesised host (M3V, pair, Mira), place the TESS PRF at that star
and predict the fraction of a unit dip falling in each r = 0.75 px aperture (M3V, pair, Mira, mirror) per sector.
Fit one free total dip amplitude per hypothesis to the measured aperture dips and compare chi2."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from wotan import flatten
from loc5 import fit_psf, Cut, circle_mask, direct_offset, component_images

def ref_offset(C, psf, O, ref=5925373204131295872):
    """WCS correction from a normal-colour bright star (the Mira's own fit is biased by its extreme colour)."""
    r = C.g[C.g.source_id == ref].iloc[0]; g = C.g
    nb = g[(np.hypot(g.x - r.x, g.y - r.y) < 4) & (g.source_id != r.source_id) & (g["T"] < r["T"] + 3)]
    F = fit_psf(O, np.sqrt(np.abs(O)) + 1, psf, r.x, r.y, rad=2.0, extra=[psf(q.x, q.y) for q in nb.itertuples()], box=1.0)
    return np.array([F["x"] - r.x, F["y"] - r.y])
E = json.load(open("out/214311906/ephemeris_v2.json")); P, T0 = E["P"], E["T0"]
cfg = json.load(open("loccfg.json"))["214311906"]; cfg = dict(cfg, P=P, t0=T0)
IDS = {"M3V": 5925373272850784768, "pair": 5925373272850782208, "Mira": 5925373268528684288}
APS = ("M3V", "pair", "Mira", "mirror")
def trap(t):
    tc = T0 + np.round((t - T0) / P) * P; a = np.abs(t - tc); T14, T23 = 0.067, 0.067 * 0.78
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
obs, err, pred = [], [], {h: [] for h in IDS}
for sec in (39, 66, 93, 103, 104):
    C = Cut(214311906, sec, cfg); psf = C.prf(); g = C.good; t = C.t[g]; cube = C.flux[g]
    comps, O = component_images(C, cfg); woff = ref_offset(C, psf, O)
    med = np.median(cube, 0); bk = np.median(cube[:, med < np.percentile(med, 25)], 1)
    pos = {k: (float(C.g[C.g.source_id == v].x.iloc[0]), float(C.g[C.g.source_id == v].y.iloc[0])) for k, v in IDS.items()}
    pos["mirror"] = (2 * pos["Mira"][0] - pos["M3V"][0], 2 * pos["Mira"][1] - pos["M3V"][1])
    masks = {a: circle_mask((C.ny, C.nx), *pos[a], 0.75) for a in APS}
    m = trap(t) == 0; y_ = trap(t)
    for a in APS:
        f = cube[:, masks[a]].sum(1) - bk * masks[a].sum()
        tr = flatten(t[m], f[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        k_ = np.isfinite(tr); r_ = f - np.interp(t, t[m][k_], tr[k_])
        A = np.column_stack([np.ones(len(t)), y_]); c, *_ = np.linalg.lstsq(A, r_, rcond=None)
        s = 1.4826 * np.median(np.abs(r_[m] - np.median(r_[m]))); obs.append(-c[1]); err.append(s / np.sqrt(np.sum(y_ ** 2)))
        for h in IDS:
            img = psf(pos[h][0] + woff[0], pos[h][1] + woff[1]); img = img / img.sum()
            pred[h].append(float(img[masks[a]].sum()))
    print(f"S{sec}: aperture pixel counts {[int(masks[a].sum()) for a in APS]}; predicted fractions " +
          "; ".join(f"{h}: " + "/".join(f"{pred[h][-4 + i]:.2f}" for i in range(4)) for h in IDS), flush=True)
obs, err = np.array(obs), np.array(err)
print("observed dips (e/s):", np.round(obs, 1).tolist())
res = {}
for h in IDS:
    x = np.array(pred[h]); w = 1 / err ** 2
    F = np.sum(w * x * obs) / np.sum(w * x * x); chi2 = float(np.sum(((obs - F * x) / err) ** 2))
    res[h] = dict(F_total=F, F_err=float(1 / np.sqrt(np.sum(w * x * x))), chi2=chi2, dof=len(obs) - 1)
    print(f"host = {h:5s}: best total dip {F:.1f} ± {1/np.sqrt(np.sum(w*x*x)):.1f} e/s; chi2 {chi2:.1f} for {len(obs)-1} dof")
json.dump(dict(obs=obs.tolist(), err=err.tolist(), pred=pred, fits=res), open("out/214311906/forward_apertures_v2.json", "w"), indent=1)
