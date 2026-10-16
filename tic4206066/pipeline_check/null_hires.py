"""Stage-2 recurrence test for TIC 4206066 with finer null distributions (stage 2 used 150 random periods)."""
import sys, json, warnings
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tess-search", "src"))  # the search pipeline, unchanged
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import chi2
import stage2
stage2.CACHE = os.environ.get("LC_CACHE", os.path.expanduser("~/.cache/tess_lc"))  # keep downloads out of the repo

r = json.load(open(os.path.join(HERE, "would_detect2.json")))
P, dur = r["s98"]["P"], r["s98"]["dur"]
secs = {s["sector"]: s for s in stage2.load_all(4206066)}
out = {}
for nnull in (150, 3000):
    faps = {}
    for sec in (6, 32):
        rng = np.random.default_rng(4206066 + sec)
        t = stage2.sector_test(secs[sec], P, dur, nnull=nnull, rng=rng)
        faps[sec] = t["fap"]
    fisher = -2 * np.sum(np.log(list(faps.values())))
    p = chi2.sf(fisher, 4)
    out[nnull] = dict(faps=faps, p_comb=p)
    print(f"nnull={nnull}: FAP S6={faps[6]:.4f} S32={faps[32]:.4f}  combined p={p:.2g}  (vetting cut 1e-3)")
json.dump(out, open(os.path.join(HERE, "null_hires.json"), "w"), indent=1)
