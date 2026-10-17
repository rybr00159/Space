"""Noise diagnostics for the DIFF and AMP per-window products: real-transit coverage, fake-window noise compared
with the formal (white) errors, pixel-to-pixel noise correlation near the Mira, and a quick look at the real images."""
import sys, time
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
from l1img import *
g = stars()
for sec in SECTORS:
    t0 = time.time()
    S = Sector(sec); Dd = Data(S)
    tr, fk = window_sets(Dd)
    # independent fakes: spaced by >= 0.7 d
    fk_ind = [fk[0]]
    for f in fk[1:]:
        if f - fk_ind[-1] >= 0.7: fk_ind.append(f)
    real_d = [(tc, Dd.diff_window(tc)) for tc in tr]
    real_d = [(tc, r) for tc, r in real_d if r is not None]
    fake_d = [Dd.diff_window(tc) for tc in fk_ind]; fake_d = [r for r in fake_d if r is not None]
    real_a = [(tc, Dd.amp_window(tc)) for tc in tr]; real_a = [(tc, r) for tc, r in real_a if r is not None]
    fake_a = [Dd.amp_window(tc) for tc in fk_ind]; fake_a = [r for r in fake_a if r is not None]
    xy = S.star_xy(g)
    im = int(np.where(g.source_id == IDS['Mira'])[0][0]); i3 = int(np.where(g.source_id == IDS['M3V'])[0][0])
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    cen = np.hypot(xx - xy[im, 0], yy - xy[im, 1]).ravel() <= 3.5
    # DIFF: normalized fake residuals
    Z = np.array([r['D'] / np.sqrt(r['V']) for r in fake_d])
    ratio = Z.std(0)
    # AMP: single-window alpha = num/den; formal sigma = sig_amp/sqrt(den)
    Za = np.array([r['num'] / r['den'] / (Dd.sig_amp / np.sqrt(r['den'])) for r in fake_a])
    ratio_a = Za.std(0)
    # pixel correlation (DIFF) among the central pixels
    Zc = Z[:, cen]; Cc = np.corrcoef(Zc.T); off = Cc[~np.eye(len(Cc), dtype=bool)]
    # real combined images
    w = np.array([1 / r['V'] for _, r in real_d]); Dm = (np.array([r['D'] for _, r in real_d]) * w).sum(0) / w.sum(0)
    De = 1 / np.sqrt(w.sum(0)) * np.median(ratio[cen])
    A = sum(r['num'] for _, r in real_a) / sum(r['den'] for _, r in real_a)
    Ae = Dd.sig_amp / np.sqrt(sum(r['den'] for _, r in real_a)) * np.median(ratio_a[cen])
    # crude totals: sum of the image within 2.5 px of the M3V (fraction of PRF ~0.9)
    ap = np.hypot(xx - xy[i3, 0], yy - xy[i3, 1]).ravel() <= 2.5
    print(f'S{sec}: real windows DIFF {len(real_d)} AMP {len(real_a)} (of {len(tr)}); fake indep {len(fake_d)}/{len(fake_a)}; '
          f'noise/formal DIFF median(central) {np.median(ratio[cen]):.2f} (all {np.median(ratio):.2f}), AMP {np.median(ratio_a[cen]):.2f} '
          f'(all {np.median(ratio_a):.2f}); DIFF central pixel corr: median {np.median(off):+.2f}, 90th pct {np.percentile(off, 90):+.2f}')
    print(f'      DIFF sum(r<=2.5 px of M3V) = {Dm[ap].sum():.1f} +- {np.sqrt((De[ap]**2).sum()):.1f} e/s; AMP = {A[ap].sum():.1f} +- {np.sqrt((Ae[ap]**2).sum()):.1f}; '
          f'max S/N pixel DIFF {np.nanmax(Dm / De):.1f} AMP {np.nanmax(A / Ae):.1f}; per-pixel sigma (central median) DIFF {np.median(De[cen]):.2f} AMP {np.median(Ae[cen]):.2f} e/s; {time.time()-t0:.0f}s')
