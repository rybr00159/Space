"""PCA of fake-window DIFF images (11x11 px around the Mira): what does the correlated noise look like?
Compares with (i) a uniform background, (ii) the Mira's PRF (flux fluctuation), (iii) PRF x/y derivatives (jitter)."""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
from l1img import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
g = stars()
fig, axs = plt.subplots(6, 5, figsize=(13, 15))
for row, sec in enumerate(SECTORS):
    S = Sector(sec); Dd = Data(S, bkg=False); psf = S.prf()   # deliberately WITHOUT background subtraction
    tr, fk = window_sets(Dd, step=0.1)
    fd = [Dd.diff_window(tc) for tc in fk]; fd = [r['D'] for r in fd if r is not None]
    F = np.array(fd)
    mira = np.array(json.load(open(f'{BASE}/out/mira_pos.json'))[str(sec)]['base']['dEN'])
    xy = S.star_xy(g); im = int(np.where(g.source_id == IDS['Mira'])[0][0])
    xm, ym = S.sky_to_xy(mira)[0]
    xc, yc = int(round(xm)), int(round(ym))
    sl = (slice(yc - 5, yc + 6), slice(xc - 5, xc + 6))
    Fc = F.reshape(len(F), S.ny, S.nx)[:, sl[0], sl[1]].reshape(len(F), -1)
    Fc0 = Fc - Fc.mean(0)
    U, s, Vt = np.linalg.svd(Fc0, full_matrices=False)
    frac = s ** 2 / (s ** 2).sum()
    # templates
    p0 = psf(xm, ym)[sl].ravel(); dx = (psf(xm + 0.02, ym) - psf(xm - 0.02, ym))[sl].ravel() / 0.04
    dy = (psf(xm, ym + 0.02) - psf(xm, ym - 0.02))[sl].ravel() / 0.04
    T = np.column_stack([np.ones_like(p0), p0, dx, dy])
    for k in range(3):
        v = Vt[k]
        c, res, *_ = np.linalg.lstsq(T, v, rcond=None)
        r2 = 1 - ((v - T @ c) ** 2).sum() / ((v - v.mean()) ** 2).sum() if True else 0
        cu = np.corrcoef(v, np.ones_like(v))[0, 1] if v.std() > 0 else 0
        print(f'S{sec} PC{k+1}: var frac {frac[k]:.2f}; mean {v.mean():+.3f} (uniform would be {1/np.sqrt(len(v)):.3f}); '
              f'R^2 of [const, Mira PRF, d/dx, d/dy] fit {r2:.2f}; coefs {np.round(c, 3)}')
        ax = axs[row, k]; vv = np.abs(v).max()
        ax.imshow(v.reshape(11, 11), origin='lower', cmap='RdBu_r', vmin=-vv, vmax=vv); ax.set_title(f'S{sec} PC{k+1} ({frac[k]*100:.0f}%)', fontsize=9)
    axs[row, 3].imshow(Fc0.std(0).reshape(11, 11), origin='lower'); axs[row, 3].set_title('fake DIFF std', fontsize=9)
    axs[row, 4].imshow(np.log10(np.median(S.flux[S.good], 0)[sl]), origin='lower'); axs[row, 4].set_title('median image', fontsize=9)
    print(f'S{sec}: n fake {len(F)}; std of fake DIFF per pixel: median {np.median(Fc0.std(0)):.2f} e/s; edge-pixel std {np.median(F.std(0)):.2f}')
fig.tight_layout(); fig.savefig(f'{BASE}/figs/noise_pca.png', dpi=70)
