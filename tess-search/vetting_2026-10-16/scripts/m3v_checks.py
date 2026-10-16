"""TIC 214311908 (M3V) QLP SAP light curve, S103/S104 times corrected: phase-0.5 depth, secondary scan, odd/even,
and the half-period test (transits at phase 0.5 would mean P/2)."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "214311906"]
import numpy as np
from wotan import flatten
import stage2
from lcan import read, bin10, in_mask, trap
E = json.load(open("out/214311906/ephemeris_v2.json")); P, T0 = E["P"], E["T0"]; shift = {int(k): v for k, v in E["qlp_shift_days"].items()}
T14, FL = 0.067, 0.78
T, F, S = [], [], []
p = stage2.products(214311908)
for _, r in p[p.prov == "QLP"].sort_values("sector").iterrows():
    sec = int(r.sector); t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), "QLP"); t = t + shift.get(sec, 0.0)
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols["SAP_FLUX"]) & (cols["SAP_FLUX"] > 0)
    tt, ff = bin10(t[ok], cols["SAP_FLUX"][ok] / np.nanmedian(cols["SAP_FLUX"][ok]), r.exptime)
    m = ~in_mask(tt, P, T0, T14, 1.0) & ~in_mask(tt, P, T0 + P / 2, T14, 1.0)
    tr = flatten(tt[m], ff[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    k_ = np.isfinite(tr); f = ff / np.interp(tt, tt[m][k_], tr[k_]); s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
    T.append(tt); F.append(f); S.append(np.full(len(tt), s))
t, f, s = np.concatenate(T), np.concatenate(F), np.concatenate(S)
def depth_at(phase0, per=P, mask_half=0.5):
    ph = ((t - T0) / per - phase0 + 0.5) % 1 - 0.5; inn = np.abs(ph * per) < mask_half * T14
    w = 1 / s[inn] ** 2; return np.sum(w * (1 - f[inn])) / w.sum(), 1 / np.sqrt(w.sum()), int(inn.sum())
d0 = depth_at(0.0); d5 = depth_at(0.5)
print(f"primary depth {d0[0]*1e6:.0f} ± {d0[1]*1e6:.0f} ppm; phase 0.5 {d5[0]*1e6:.0f} ± {d5[1]*1e6:.0f} ppm (ratio {d5[0]/d0[0]:.2f})")
ep = np.round((t - T0) / P)
for par, lab in ((0, "even"), (1, "odd")):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5; inn = (np.abs(ph * P) < 0.5 * T14) & (ep % 2 == par)
    w = 1 / s[inn] ** 2; print(f"  {lab}: {np.sum(w*(1-f[inn]))/w.sum()*1e6:.0f} ± {1/np.sqrt(w.sum())*1e6:.0f} ppm")
scan = []
for ph0 in np.arange(0.05, 0.951, 0.0025):
    d, e, n = depth_at(ph0)
    if n > 5: scan.append((ph0, d / e, d))
a = np.array(scan); i = int(np.argmax(a[:, 1]))
print(f"secondary scan: max S/N {a[i,1]:.1f} at phase {a[i,0]:.3f} (depth {a[i,2]*1e6:.0f} ppm); rms of scan S/N {np.std(a[:,1]):.2f}; n trials {len(a)}")
json.dump(dict(primary=d0, half=d5, scan_max=a[i].tolist()), open("out/214311906/m3v_checks.json", "w"), indent=1, default=float)

# is the phase-0.825 feature coherent across sectors and events?
sec_of = np.concatenate([np.full(len(x), k) for k, x in zip(range(len(T)), T)])
names = [int(r.sector) for _, r in p[p.prov == "QLP"].sort_values("sector").iterrows()]
for ph0 in (a[i, 0], 0.0):
    ph = ((t - T0) / P - ph0 + 0.5) % 1 - 0.5; inn = np.abs(ph * P) < 0.5 * T14
    print(f"phase {ph0:.3f}: per-sector depth (ppm):", end=" ")
    for k, nm in enumerate(names):
        sel = inn & (sec_of == k)
        if sel.sum() > 3:
            w = 1 / s[sel] ** 2; print(f"S{nm} {np.sum(w*(1-f[sel]))/w.sum()*1e6:.0f}±{1/np.sqrt(w.sum())*1e6:.0f}", end="; ")
    print()
    ev = np.round((t[inn] - T0) / P - ph0); dep = []
    for e in np.unique(ev):
        sel = ev == e; w = 1 / s[inn][sel] ** 2; dep.append((int(e), np.sum(w * (1 - f[inn][sel])) / w.sum(), 1 / np.sqrt(w.sum())))
    dep = np.array(dep); top = dep[np.argsort(-dep[:, 1] / dep[:, 2])][:3]
    print(f"   {len(dep)} events; strongest single events (epoch, ppm, S/N): " + "; ".join(f"{int(x[0])} {x[1]*1e6:.0f} {x[1]/x[2]:.1f}" for x in top))
