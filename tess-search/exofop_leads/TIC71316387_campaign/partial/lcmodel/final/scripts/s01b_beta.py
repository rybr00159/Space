"""Red-noise check: un-floored time-averaging beta per sector for the default set, computed per window (no concatenation
across window boundaries) at 30/45/60/90 min; also the same for fake windows at phases 0.25/0.75."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import lclib as L
for lc in L.default_set() + [L.load_qlp(95), L.load_qlp(105)]:
    lc = L.clip_outliers(lc)
    for phs in [(0.0, 0.5), (0.25, 0.75)]:
        ws = L.windows(lc, phases=phs)
        R = []
        for w in ws:
            i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
            if phs[0] == 0.0:
                c = np.polyfit(dt[oot], y[oot], 2); R.append((y[oot] / np.polyval(c, dt[oot]) - 1, dt[oot]))
            else:
                c = np.polyfit(dt, y, 2); R.append((y / np.polyval(c, dt) - 1, dt))
        s1 = np.std(np.concatenate([r for r, _ in R]))
        out = []
        for sc in (30, 45, 60, 90):
            n = int(round(sc / (lc['texp'] * 1440)))
            if n < 2: out.append(np.nan); continue
            bm = []
            for r, dt in R:
                # bin within contiguous segments only
                m = r.size // n
                if m < 2: continue
                bm.append(r[:m * n].reshape(m, n).mean(1))
            bm = np.concatenate(bm)
            out.append(np.std(bm) / (s1 / np.sqrt(n)))
        print(f"S{lc['sector']:3d} {lc['prod']:9s} phases {phs}: sigma {s1*1e6:6.0f} ppm  beta(30,45,60,90 min) = {np.round(out,2)}")
