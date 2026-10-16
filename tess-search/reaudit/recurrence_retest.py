"""Re-run the stage-2 recurrence test with a finer null for signals that missed p < 1e-3 only because of its resolution.

stage2.py compared each other sector with 150 random periods, so no sector could score below 1/151 and a star with one or
two other sectors could rarely reach the p < 1e-3 vetting cut. This re-tests the Sector 104 signals that missed the cut
with at least one sector at that floor, using 3,000 random periods per sector. The discovery sector is excluded, as in
stage 2. Usage: python recurrence_retest.py [nnull] [procs]
"""
import os, sys, json, warnings
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import chi2
import stage2
stage2.CACHE = os.environ.get("LC_CACHE", os.path.expanduser("~/.cache/tess_lc"))

# Sector 104 stage-2 inputs (period and box duration in days) from stage2_in_s0104.csv
CANDS = [
    dict(tic=90086266, P=2.360598141113784, dur=0.2707744926571693),
    dict(tic=61332759, P=0.7397897406107939, dur=0.1175366877605934),
    dict(tic=214311906, P=5.687252960485503, dur=0.0624972852800604),
    dict(tic=55829034, P=0.6157542375343312, dur=0.0553486955087039),
    dict(tic=61141412, P=6.52516823484286, dur=0.034708341674696),
]
NNULL = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def one(c):
    tic, P, dur = int(c["tic"]), float(c["P"]), float(c["dur"])
    secs = [s for s in stage2.load_all(tic) if s["sector"] != 104]
    per = []
    for s in secs:
        r = stage2.sector_test(s, P, dur, nnull=NNULL, rng=np.random.default_rng(tic + s["sector"]))
        if np.isfinite(r["snr"]):
            per.append(dict(sector=s["sector"], prov=s["prov"], snr=float(r["snr"]), fap=float(r["fap"])))
    fisher = -2 * np.sum(np.log([x["fap"] for x in per]))
    return dict(tic=tic, P=P, nnull=NNULL, per=per, p_comb=float(chi2.sf(fisher, 2 * len(per))))


if __name__ == "__main__":
    from multiprocessing import get_context
    c = CANDS
    with get_context("spawn").Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as pool:
        res = pool.map(one, c)
    json.dump(res, open(os.path.join(HERE, "recurrence_retest.json"), "w"), indent=1)
    for r in sorted(res, key=lambda r: r["p_comb"]):
        print(r["tic"], f"P={r['P']:.4f}", f"p_comb={r['p_comb']:.2g}", [(x["sector"], round(x["fap"], 4)) for x in r["per"]])
