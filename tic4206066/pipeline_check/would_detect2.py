"""Stage 1 on S6/S32 (QLP and SPOC-family products), then classification and stage 2 for the S98 detection."""
import sys, json, warnings, io
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tess-search", "src"))  # the search pipeline, unchanged
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import fetch_qlp, stage2
stage2.CACHE = os.environ.get("LC_CACHE", os.path.expanduser("~/.cache/tess_lc"))  # keep downloads out of the repo
from bls import detrend, search_robust, clean, diagnostics
from classify import classify

TIC = 4206066
EPH = json.load(open(os.path.join(HERE, "..", "sources", "ephemeris_final.json")))
A, B = EPH["A"], EPH["B"]

def stage1(t, f, label):
    t, f = detrend(t, f)
    best, _ = search_robust(t, f, max_mask=2)
    pre = best["snr_red"] >= 8.5 and best["passed"]
    post = best["snr_red"] >= 9 and best["n_tr"] >= 3 and best["max_event_frac"] <= 0.5 and best["depth"] > 0
    gaps = np.diff(t); big = np.where(gaps > 1.0)[0]
    print(f"\n== {label}: {len(t)} pts, BTJD {t.min():.1f}-{t.max():.1f}, gaps>1 d: "
          f"{[(round(t[i],1), round(t[i+1],1)) for i in big]}, sigma {best['sigma']*1e6:.0f} ppm")
    print(f"   winner: P={best['P']:.5f} d depth={best['depth']*1e6:.0f} ppm dur={best['dur']*24:.2f} h n_tr={best['n_tr']} "
          f"snr_red={best['snr_red']:.2f} SDE={best['sde']:.1f} alias_of={best.get('alias_of')}")
    print(f"   passes stage-1 cuts: {pre and post}")
    tc, fc, s = clean(t, f)
    d = diagnostics(tc, fc, A["P"], A["T0"], A["T14"], s)
    print(f"   at published 3.18 d ephemeris: depth={d['depth']*1e6:.0f}±{d['depth_err']*1e6:.0f} ppm n_tr={d['n_tr']} snr_red={d['snr_red']:.2f}")
    return best

# S98 detection again (kept for classification / stage 2)
t, f, n, meta = fetch_qlp.load(TIC, 98, binsize=3)
b98 = stage1(t, f, "S98 QLP DET_FLUX")

# S6 and S32 via the same product loader stage 2 uses (SPOC 2-min > TESS-SPOC FFI > QLP)
secs = stage2.load_all(TIC)
print("\nproducts used by stage 2:", [(s["sector"], s["prov"], s["exptime"]) for s in secs])
p = stage2.products(TIC)
print(p[["sector", "prov", "exptime"]].to_string(index=False))
for prov in ("QLP", "TESS-SPOC", "SPOC"):
    for _, r in p[p.prov == prov].sort_values("sector").iterrows():
        if r.sector == 98 and prov == "QLP": continue
        path = stage2.fetch(r.uri)
        lc = stage2.read_lc(path, prov) if path else None
        if lc is None: print(f"\nS{r.sector} {prov}: unreadable"); continue
        stage1(*lc, f"S{r.sector} {prov} {r.exptime:.0f}s (binned to >=10 min)")

# classification as in classify.py
tic = pd.DataFrame([dict(ID=TIC, Tmag=9.6817, Teff=4010.0, rad=0.648133, mass=0.63, d=35.6156, ra=80.8598, dec=-1.3273, contratio=0.0027)])
t2, f2 = detrend(t, f); tc, fc, s = clean(t2, f2)
d = diagnostics(tc, fc, b98["P"], b98["t0"], b98["dur"], s)
c = pd.DataFrame([dict(tic=TIC, P=b98["P"], dur=b98["dur"], red_ratio=b98["red_ratio"], **{f"m_{k}": v for k, v in d.items()})])
c = classify(c, tic)
print("\nclassification:", c[["track", "depth_ppm", "rp_re", "dur_ratio", "m_oddeven_sig", "m_sec_sig", "m_oot_rms_ratio",
                              "m_oot_dip_sig", "m_oot_dip_ratio", "red_ratio"]].round(2).to_dict("records"))

# stage 2: recurrence in the other sectors, discovery sector excluded
r = stage2._run(dict(tic=TIC, P=b98["P"], dur=b98["dur"], sector="98"))
print("\nstage 2: sectors", r["sectors"], " other sectors with FAP<=1%:", r["n_other_sig"], "/", r["n_other"],
      " combined p:", f"{r['p_comb']:.2g}")
for x in r["per"]:
    print(f"   S{x['sector']} {x['prov']}: narrow S/N {x['snr']:.1f}, FAP {x['fap']:.3f}, depth at joint ephemeris "
          f"{x['depth_joint']*1e6:.0f}±{x['depth_err']*1e6:.0f} ppm ({x['ntr_joint']} transits)")
j = r["joint"]
print(f"   joint: P={j['P']:.6f} d (published {A['P']:.6f}), S/N {j['snr']:.1f}, next-best alias S/N {j['second']:.1f}")
if r["others_only"]:
    o = r["others_only"]; print(f"   S6+S32 only: P={o['P']:.6f}, S/N {o['snr']:.1f}, next-best alias {o['second']:.1f}")
json.dump(dict(s98=b98, stage2=r), open(os.path.join(HERE, "would_detect2.json"), "w"), default=float, indent=1)
