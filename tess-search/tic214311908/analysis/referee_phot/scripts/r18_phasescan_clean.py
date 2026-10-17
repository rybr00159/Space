"""Referee r18: clean re-aggregation of the r08 phase scan. r08 windows (half-width 0.35 d; 0.4 d S39; 0.6 d S12)
centred at phases 0.06-0.94 can overlap the real transit near phases 0.06 and 0.94, which biases the per-pixel
baseline and produces spurious 'dips'. Here every window whose span (centre +- half + T14/2) reaches a real transit is
dropped, the empirical single-window error per sector is recomputed from the clean windows, and the phase scan, the
phase-0.5 depth, and the scan maximum are recomputed. The look-elsewhere expectation for the maximum S/N is computed
by permuting the window depths among phases within each sector (500 permutations, keeping the per-phase window count).
Output: out/phase_scan_clean.csv, out/phase_scan_clean.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import half_for
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T0, P, T14 = G["T0"], G["P"], G["T14_h"] / 24
SC = pd.read_csv(f"{R}/out/phase_scan_windows.csv")
dist = np.abs(((SC.tc - T0) / P + 0.5) % 1 - 0.5) * P
lim = SC.sector.map(half_for) + T14 / 2 + 0.01
SC["clean"] = dist > lim
C = SC[SC.clean].copy()
sig = {s: rstd(C[C.sector == s].depth) for s in SECTORS}
C["e"] = C.sector.map(sig)


def agg(D):
    g = D.groupby("phase")
    w = 1 / D.e ** 2
    num = (D.depth * w).groupby(D.phase).sum(); den = w.groupby(D.phase).sum()
    a = pd.DataFrame(dict(depth=num / den, err=1 / np.sqrt(den), n=g.size())).reset_index()
    a["snr"] = a.depth / a.err; return a


A = agg(C); A.to_csv(f"{R}/out/phase_scan_clean.csv", index=False)
i = int(np.argmax(A.snr)); p5 = A.iloc[np.argmin(np.abs(A.phase - 0.5))]
# biased (dropped) windows: what did they show?
B = SC[~SC.clean].copy(); B["e"] = B.sector.map(sig); Ab = agg(B) if len(B) else None
rng = np.random.default_rng(7); mx = []
for k in range(500):
    D = C.copy()
    for s in SECTORS:
        m = (D.sector == s).values; D.loc[m, "depth"] = rng.permutation(D.depth.values[m])
    mx.append(agg(D).snr.max())
mx = np.array(mx)
out = dict(n_windows_total=len(SC), n_clean=len(C), phase_range_clean=[float(A.phase.min()), float(A.phase.max())],
           single_window_err=sig, phase05=dict(depth_e_s=float(p5.depth), err_e_s=float(p5.err), snr=float(p5.snr), depth_ppt_F979=float(p5.depth / 979 * 1e3), err_ppt_F979=float(p5.err / 979 * 1e3)),
           scan_max_snr=float(A.snr[i]), scan_max_phase=float(A.phase[i]), scan_max_depth=float(A.depth[i]), scan_min_snr=float(A.snr.min()),
           scan_snr_std=float(A.snr.std()), scan_snr_rstd=float(rstd(A.snr)), perm_max_snr_median=float(np.median(mx)), perm_p_value=float(np.mean(mx >= A.snr[i])),
           dropped_windows_max_snr=(float(Ab.snr.max()) if Ab is not None else None), dropped_windows_max_phase=(float(Ab.phase[np.argmax(Ab.snr)]) if Ab is not None else None))
save_json(out, f"{R}/out/phase_scan_clean.json")
print(json.dumps(out, indent=0))
