"""Per-transit (single-window) fits: amplitude on the M3V / pair / Mira and the single-window free-fit position;
flags windows that disagree, and lists windows rejected for coverage."""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from l1img import window_sets
rows = []
for sec in SECTORS:
    S = Sector(sec); Dd = Data(S)
    tr, _ = window_sets(Dd)
    for meth in ('DIFF', 'AMP'):
        M = SectorModel(sec, meth, S=S, Dd=Dd)
        used = [tc for tc, _ in M.real]
        miss = [round(tc, 3) for tc in tr if all(abs(tc - u) > 1e-6 for u in used) and Dd.t.min() - 0.05 < tc < Dd.t.max() + 0.05]
        for (tc, _), v in zip(M.real, M.rv):
            d, s = M._single(v)
            # use the sector sigma map scaled to this window's noise
            c2, a, b, sa = M.fit(d, 'H')
            cg, *_ = M.fit(d, 'G'); pos, cmin = M.grid_refine(cg)
            h = M.hyp; ii = {n: np.where(h.name == n)[0][0] for n in ('M3V', 'pairA', 'Mira')}
            k = int(np.round((tc - T0) / P))
            scale = np.median(s * M.infl_pix / M.sig)
            rows.append(dict(sec=sec, meth=meth, epoch=k, tc=tc, a_M3V=a[ii['M3V']], sa=sa[ii['M3V']] * scale,
                             d_pair=(c2[ii['pairA']] - c2[ii['M3V']]) / scale ** 2, d_mira=(c2[ii['Mira']] - c2[ii['M3V']]) / scale ** 2,
                             free_E=pos[0], free_N=pos[1]))
            r = rows[-1]
            print(f"S{sec} {meth} E{k:4d} tc {tc:.3f}: a(M3V) {r['a_M3V']:5.1f}+-{r['sa']:.1f}  chi2(pair)-chi2(M3V) {r['d_pair']:+6.1f}  chi2(Mira)-chi2(M3V) {r['d_mira']:+6.1f}  free ({pos[0]:+.0f},{pos[1]:+.0f})")
        if miss: print(f'S{sec} {meth}: transits without enough coverage: {miss}')
pd.DataFrame(rows).to_csv(f'{BASE}/out/per_window.csv', index=False)
