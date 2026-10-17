"""Quick look at the real data with l2fit: per-sector chi2 of the main hypotheses, free-fit sky position (M-dwarf
colour), amplitudes; both methods."""
import sys, time
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
for sec in SECTORS:
    t0 = time.time()
    S = Sector(sec); Dd = Data(S)
    for meth in ('DIFF', 'AMP'):
        M = SectorModel(sec, meth, S=S, Dd=Dd)
        c2, a, b, sa = M.fit(M.d_real, 'H')
        cg, ag, _, _ = M.fit(M.d_real, 'G')
        pos, cmin = M.grid_refine(cg)
        nul = M.null_chi2(M.d_real)
        h = M.hyp
        def row(nm):
            i = np.where(h.name == nm)[0][0]; return f'{nm} dchi2 {c2[i]-cmin:6.1f} a {a[i]:5.1f}+-{sa[i]:.1f}'
        order = np.argsort(c2)[:3]
        print(f'S{sec} {meth}: n_real {len(M.rv)} n_fake {len(M.fv)} infl {M.infl_med:.2f} npix {M.npix}; null-chi2 - min {nul-cmin:.1f}; '
              f'free fit (dE,dN)=({pos[0]:+.1f},{pos[1]:+.1f})"; ' + '; '.join(row(n) for n in ('M3V', 'pairA', 'Mira')) +
              f'; best stars: ' + ', '.join(f"{h.name[i]}({h.sep[i]:.0f}\", T{h['T'][i]:.1f}) {c2[i]-cmin:.1f}" for i in order) + f'  [{time.time()-t0:.0f}s]')
