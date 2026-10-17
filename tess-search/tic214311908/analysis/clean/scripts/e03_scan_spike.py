"""e03 (exploratory): identify the window responsible for the single-step REF phase-scan spike at phase 0.5375."""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from c08_phasescan import *
for ph in (0.535, 0.5375, 0.54):
    for s in SECTORS:
        S = SectorC(s, "all"); t = S.ttdb[S.qok]
        E = np.arange(np.floor((t.min() - T0) / P) - 1, np.ceil((t.max() - T0) / P) + 1)
        for e in E:
            tc = T0 + P * (e + ph)
            if tc < t.min() or tc > t.max(): continue
            w = S.window(tc, half=HALF)
            if w is None: continue
            tt = w["t"]; ff = 1 + w["aM"] / F_M3V
            if not accept(tt, tc, HALF, S.exptime): continue
            d, de, n = common_depth([(tt, ff, tc, S.exptime)])
            if abs(d / de) > 2.5: print(f"phase {ph} S{s} E{int(e)} tc {tc:.3f}: depth {d:.1f} +- {de:.1f} ppt (S/N {d/de:.1f}); n {len(tt)}, n_oot {w['n_oot']}, B max {np.max(w['B']):.0f}")
