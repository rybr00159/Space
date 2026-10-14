"""v4 real localisation per sector: matched-filter amplitude images, stacked per eclipse type, TESS-PRF fits.
offset = J @ (pos(amplitude image) - pos(direct image)); direct image fitted with Gaia neighbours tied (T<19).
Also per-event offsets and a few analysis variants (env: V4_TAG, V4_BOX, V4_BASEDEG, V4_SHAPE, V4_W, V4_PLANE)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

TAG = os.environ.get('V4_TAG', '')
PLANE = os.environ.get('V4_PLANE', '0') == '1'
RATIO = 0.47


def windows(S):
    fn = f'{V4}/results/win_s{S.sec}{TAG}.npz'
    if os.path.exists(fn) and os.environ.get('V4_REDO', '0') != '1':
        z = np.load(fn, allow_pickle=True)
        return {k: z[k] for k in z.files}
    W = S.all_windows()
    np.savez_compressed(fn, **W)
    return W


def main(sec):
    t0 = time.time()
    S = Sector(sec)
    W = windows(S)
    G = PRFGrid(S)
    kind = W['kind']
    out = dict(sector=sec, cam=S.cam, ccd=S.ccd, J=S.J, xt_wcs=S.xt, yt_wcs=S.yt, n_cut=S.ncut, texp_min=S.texp * 1440,
               epoch=S.epoch, box=BOX, basedeg=BASEDEG, shape=SHAPE, whalf=WHALF, plane=PLANE,
               n_windows={k: int((kind == k).sum()) for k in ('deep', 'shallow', 'fake')})
    sets = dict(deep=(np.where(kind == 'deep')[0], 1.0), shallow=(np.where(kind == 'shallow')[0], 1.0))
    ib = np.r_[sets['deep'][0], sets['shallow'][0]]
    sets['both'] = (ib, np.r_[np.ones(len(sets['deep'][0])), RATIO * np.ones(len(sets['shallow'][0]))])
    imgs = {}
    for name, (idx, rel) in sets.items():
        A, sA, B = stack(W, idx, rel)
        rB = G.fit(B, sA, direct=True, plane=PLANE)
        rA = G.fit(A, sA, plane=PLANE)
        rc = G.fit_cont(A, sA, rA['x'], rA['y'])
        off = S.J @ np.array([rA['x'] - rB['x'], rA['y'] - rB['y']])
        offc = S.J @ np.array([rc['x'] - rB['x'], rc['y'] - rB['y']])
        # chi2 with the source fixed at the direct-image position (target hypothesis)
        b = A[G.sl].ravel(); w = 1 / sA[G.sl].ravel()
        m = G.img(rB['x'], rB['y'])[G.sl].ravel()
        X = np.column_stack([m, np.ones_like(m)]) * w[:, None]
        cc, *_ = np.linalg.lstsq(X, b * w, rcond=None)
        chiT = float(((X @ cc - b * w) ** 2).sum())
        out[name] = dict(n=int(len(idx)), xB=rB['x'], yB=rB['y'], FB=rB['a'], xA=rA['x'], yA=rA['y'], amp=rA['a'],
                         dE=float(off[0]), dN=float(off[1]), dE_cont=float(offc[0]), dN_cont=float(offc[1]),
                         chi2=rA['chi2'], npix=rA['npix'], dchi2_target=chiT - rA['chi2'], edge=rA['edge'],
                         depth_peak_ppm=rA['a'] / rB['a'] * 1e6, snr_peak=float(cc[0] / np.sqrt(np.linalg.inv(X.T @ X)[0, 0])))
        imgs[name] = (A, sA, B, rA, rB)
        print(f"S{sec} {name:7s} n={len(idx):2d} depth_pk={out[name]['depth_peak_ppm']:6.0f}ppm S/N={out[name]['snr_peak']:5.1f} "
              f"off=({off[0]:+.2f},{off[1]:+.2f})\" cont=({offc[0]:+.2f},{offc[1]:+.2f}) dchi2T={chiT - rA['chi2']:.1f} "
              f"chi2/npix={rA['chi2']:.0f}/{rA['npix']}", flush=True)
    # per-event offsets (direct position from the deep+shallow direct stack)
    _, _, Ball = stack(W, ib)
    sAll = stack(W, ib)[1]
    rB = G.fit(Ball, sAll, direct=True)
    per = []
    for i in ib:
        rA = G.fit(W['A'][i], W['sA'][i])
        off = S.J @ np.array([rA['x'] - rB['x'], rA['y'] - rB['y']])
        per.append(dict(kind=str(kind[i]), tc=float(W['tc'][i]), dE=float(off[0]), dN=float(off[1]), amp=rA['a'], edge=rA['edge'],
                        w=float(1 / np.median(W['sA'][i] ** 2))))
    out['per_event'] = per
    to_json(out, f'{V4}/results/real_s{sec}{TAG}.json')
    fig, axs = plt.subplots(2, 4, figsize=(15, 7.5))
    for row, name in enumerate(('deep', 'shallow')):
        A, sA, B, rA, rB = imgs[name]
        ax = axs[row]
        ax[0].imshow(np.log10(np.clip(B - np.median(B), 1, None)), origin='lower', cmap='gray_r'); ax[0].set_title(f'S{sec} direct')
        v = np.nanmax(np.abs(A[G.sl])); ax[1].imshow(A, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v); ax[1].set_title(f'{name} amplitude')
        ax[2].imshow(A / sA, origin='lower', cmap='RdBu_r', vmin=-8, vmax=8); ax[2].set_title('A / sigma')
        mod = rA['a'] * G.img(rA['x'], rA['y'])
        r = (A - mod); r -= np.median(r[G.sl])
        ax[3].imshow(r / sA, origin='lower', cmap='RdBu_r', vmin=-4, vmax=4); ax[3].set_title('resid / sigma')
        for a in ax:
            a.plot(rB['x'], rB['y'], '+', color='lime', ms=14, mew=2); a.plot(rA['x'], rA['y'], 'x', color='m', ms=10, mew=2)
            for dx, dy, rr, _ in G.nb:
                if rr > 0.003:
                    a.plot(rB['x'] + dx, rB['y'] + dy, 'o', mfc='none', mec='orange')
            a.set_xlim(S.ix - 5.5, S.ix + 5.5); a.set_ylim(S.iy - 5.5, S.iy + 5.5)
    plt.tight_layout(); plt.savefig(f'{V4}/plots/real_s{sec}{TAG}.png', dpi=60); plt.close()
    print(f'S{sec} done [{time.time() - t0:.0f}s]', flush=True)


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        main(s)
