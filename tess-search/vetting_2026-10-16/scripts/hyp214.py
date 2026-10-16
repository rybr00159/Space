"""Which star hosts the 5.69 d dips near TIC 214311906? For each sector: the observed free-fit source position, and
the positions recovered when the same signal is injected (at other phases) at the Mira, the M-dwarf pair TIC 214311910
or the M3V TIC 214311908. Per-hypothesis likelihood from the scatter of the injections."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from loc5 import Cut, component_images, fit_psf_multi, direct_offset, template
cfg = json.load(open("loccfg.json"))["214311906"]
cfg = dict(cfg, P=5.691471, t0=3840.9101)
H = {"Mira": 5925373268528684288, "pair TIC 214311910": 5925373272850782208, "M3V TIC 214311908": 5925373272850784768}
rng = np.random.default_rng(1)
rows = []
for sec in (39, 66, 93, 103, 104, 12):
    C = Cut(214311906, sec, cfg); psf = C.prf()
    comps, O = component_images(C, cfg); woff = direct_offset(C, psf, O)
    q = (comps[0][0] / comps[0][1]); yy, xx = np.unravel_index(np.nanargmax(q), q.shape)
    F = fit_psf_multi(comps, psf, float(xx), float(yy))
    obs = C.J @ (np.array([F["x"], F["y"]]) - woff - np.array([C.xt, C.yt]))
    # signal amplitude (e/s at PRF peak scale) from a fit at the observed position, used for injections
    A0 = comps[0][0]; amp = np.nansum(A0[np.hypot(*np.meshgrid(np.arange(C.nx) - F["x"], np.arange(C.ny) - F["y"])) < 3])
    row = dict(sector=sec, obs_dE=obs[0], obs_dN=obs[1])
    for name, gid in H.items():
        r = C.g[C.g.source_id == gid].iloc[0]
        rec = []
        for k in range(8):
            cfg2 = dict(cfg); cfg2["t0"] = cfg["t0"] + cfg["P"] * (0.2 + 0.6 * rng.random())
            fl0 = C.flux.copy(); tpl = template(C.t, cfg2)
            img = psf(r.x + woff[0], r.y + woff[1]); img = img / img.sum()
            C.flux = fl0 - amp * tpl[:, None, None] * img[None]
            c2, _ = component_images(C, cfg2); C.flux = fl0
            q2 = c2[0][0] / c2[0][1]; y2, x2 = np.unravel_index(np.nanargmax(q2), q2.shape)
            F2 = fit_psf_multi(c2, psf, float(x2), float(y2))
            rec.append(C.J @ (np.array([F2["x"], F2["y"]]) - woff - np.array([C.xt, C.yt])))
        rec = np.array(rec); mu = rec.mean(0); cov = np.cov(rec.T) + np.eye(2) * 1.0
        d = obs - mu; chi2 = float(d @ np.linalg.solve(cov, d))
        row[f"{name}: mean"] = np.round(mu, 1).tolist(); row[f"{name}: sd"] = np.round(np.sqrt(np.diag(cov)), 1).tolist(); row[f"{name}: chi2"] = round(chi2, 2)
    rows.append(row)
    print(f"S{sec}: observed ({obs[0]:+.1f}, {obs[1]:+.1f})\" | " + " | ".join(f"{n}: rec {row[n + ': mean']} sd {row[n + ': sd']} chi2 {row[n + ': chi2']}" for n in H), flush=True)
df = pd.DataFrame(rows)
for n in H: print(f"{n}: total chi2 {df[n + ': chi2'].sum():.1f} for {2 * len(df)} dof")
df.to_csv("out/214311906/host_hypotheses.csv", index=False)
