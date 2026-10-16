"""Does the S6 QLP BLS winner (11.12 d, 2 events) line up with the published 11.13 d ephemeris?"""
import sys, json, warnings
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tess-search", "src"))  # the search pipeline, unchanged
warnings.filterwarnings("ignore")
import numpy as np
import stage2
stage2.CACHE = os.environ.get("LC_CACHE", os.path.expanduser("~/.cache/tess_lc"))  # keep downloads out of the repo
from bls import detrend, search_robust
B = json.load(open(os.path.join(HERE, "..", "sources", "ephemeris_final.json")))["B"]
A = json.load(open(os.path.join(HERE, "..", "sources", "ephemeris_final.json")))["A"]
p = stage2.products(4206066)
r = p[(p.sector == 6) & (p.prov == "QLP")].iloc[0]
t, f = stage2.read_lc(stage2.fetch(r.uri), "QLP")
t, f = detrend(t, f)
best, _ = search_robust(t, f, max_mask=2)
ev = best["t0"] + best["P"] * np.arange(-3, 4); ev = ev[(ev > t.min()) & (ev < t.max())]
print("S6 QLP winner events (BTJD):", np.round(ev, 3), "P", round(best["P"], 4), "dur h", round(best["dur"] * 24, 2))
for name, e in (("B", B), ("A", A)):
    k = np.round((ev - e["T0"]) / e["P"])
    pred = e["T0"] + k * e["P"]
    sig = np.sqrt(e["cov_T0_P"][0][0] + k ** 2 * e["cov_T0_P"][1][1] + 2 * k * e["cov_T0_P"][0][1]) if "cov_T0_P" in e else np.nan
    print(f"  nearest predicted {name} transits:", np.round(pred, 3), " offset (h):", np.round((ev - pred) * 24, 2), " 1-sigma pred (h):", np.round(sig * 24, 2))
