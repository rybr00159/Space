"""r06b: is the recurring Delta-chi2 = 119.5 maximum of the analyst's inverted+rolled nulls a wrap-junction artefact?

Rolling a sector's flux array relative to its fixed time stamps (analyst s02 roll_within) makes the last and the first
points of the sector contiguous wherever the junction lands. If those edge points deviate (detrending edge effects),
the junction is a fixed 'event' present in every rolled realisation, regardless of the roll offset.
Here: the analyst's search input (search/data/m3v_clean.npz: det_clean, 5.69-d windows removed, 10-min bins, per-sector
median and robust weights, exactly as s01 search_data) is used; for each sector the inverted junction segment
(last n + first n bins, n = 1..30) is scored with the single-event box statistic S^2/W (S = -sum w y over the
segment, after inversion), and the best duration and depth are reported. Also the edge points' mean residuals.
Output: out/r06b_junction.json
RESULT: refuted (junction segments give Delta-chi2 <= 7). With argument "gap": gap_events() finds the real cause
(out/r06b_gap_events.json).
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *


def main():
    z = np.load(f"{BASE}/search/data/m3v_clean.npz")
    m = ~z["intr"]
    tb, yb, nb, sb = binlc(z["t"][m], z["det_clean"][m] - 1, z["sector"][m], z["cad"][m].astype(float), 10.0)
    w = np.zeros(len(yb)); y = yb.copy()
    for q in np.unique(sb):
        k = sb == q; y[k] -= np.median(y[k]); w[k] = 1 / rstd(y[k]) ** 2
    y -= np.sum(w * y) / np.sum(w)
    yi = -y
    out = {}
    for q in np.unique(sb):
        idx = np.where(sb == q)[0]
        best = (0, None)
        for n1 in range(0, 31):
            for n2 in range(0, 31):
                if n1 + n2 < 2: continue
                seg = np.concatenate([idx[len(idx) - n1:] if n1 else np.array([], int), idx[:n2]])
                S = -np.sum(w[seg] * yi[seg]); W = np.sum(w[seg])
                if S > 0 and S * S / W > best[0]:
                    best = (S * S / W, dict(n_end=n1, n_start=n2, depth_ppt=float(S / W * 1e3), dur_h=float((n1 + n2) * max(600, np.median(np.diff(tb[idx])) * 86400) / 3600)))
        # also single-segment (no junction) maximum within the sector for comparison: best inverted box of the same lengths
        out[int(q)] = dict(junction_dchi2=float(best[0]), **(best[1] or {}),
                           edge_end_mean_ppt=float(np.mean(y[idx[-12:]]) * 1e3), edge_start_mean_ppt=float(np.mean(y[idx[:12]]) * 1e3),
                           sector_rms_ppt=float(rstd(y[idx]) * 1e3))
        print(q, out[int(q)], flush=True)
    jdump(out, f"{ME}/out/r06b_junction.json")


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


def gap_events():
    """Index-contiguous (ignoring time gaps) inverted boxes of 2-30 bins: the largest S^2/W, and whether the box spans
    a time gap (> 30 min) in the real data. In rolled nulls, flux runs are re-assigned to contiguous time stamps, so
    an excursion split by a gap in the real data is fully exposed in most realisations."""
    z = np.load(f"{BASE}/search/data/m3v_clean.npz")
    m = ~z["intr"]
    tb, yb, nb, sb = binlc(z["t"][m], z["det_clean"][m] - 1, z["sector"][m], z["cad"][m].astype(float), 10.0)
    w = np.zeros(len(yb)); y = yb.copy()
    for q in np.unique(sb):
        k = sb == q; y[k] -= np.median(y[k]); w[k] = 1 / rstd(y[k]) ** 2
    y -= np.sum(w * y) / np.sum(w)
    res = []
    for sign, nm in ((-1, "inverted"), (1, "real")):
        yi = sign * y
        cw = np.concatenate([[0], np.cumsum(w)]); cy = np.concatenate([[0], np.cumsum(w * yi)])
        best = []
        for L in range(2, 31):
            S = -(cy[L:] - cy[:-L]); W = cw[L:] - cw[:-L]
            d2 = np.where(S > 0, S * S / W, 0)
            j = int(np.argmax(d2))
            idx = np.arange(j, j + L)
            same_sec = len(np.unique(sb[idx])) == 1
            span_h = (tb[idx[-1]] - tb[idx[0]]) * 24
            maxgap_h = float(np.max(np.diff(tb[idx])) * 24)
            best.append(dict(sign=nm, L=L, dchi2=float(d2[j]), t_start=float(tb[idx[0]]), sector=int(sb[j]), same_sector=bool(same_sec),
                             span_h=float(span_h), max_gap_h=maxgap_h, depth_ppt=float(S[j] / W[j] * 1e3),
                             near_transit_window=float(np.min(np.abs(phase_dist(tb[idx])))) ))
        bb = sorted(best, key=lambda r: -r["dchi2"])[:3]
        res += bb
        for r in bb: print(r)
    jdump(res, f"{ME}/out/r06b_gap_events.json")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "gap":
    gap_events()
