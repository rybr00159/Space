"""Real difference-image localisation per sector (deep, shallow, both), v3 pipeline.
Builds and caches the per-window difference images (real + fake windows), then:
  direct image B (stack) -> PRF fit with tied Gaia neighbours -> target position (xT, yT)
  difference image D (stack) -> free PRF + constant -> (xD, yD);  offset = J @ (xD - xT, yD - yT)
Also: per-eclipse offsets (each real window alone, direct position from the stack), a Gaussian-PSF variant
(widths fitted to the direct image) and a planar-background variant.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3lib import *
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RATIO = float(os.environ.get('V3_RATIO', 0.48))
TAG = os.environ.get('V3_TAG', '')


def get_windows(sec):
    fn = f'{V3}/results/win_s{sec}{TAG}.npz'
    if os.path.exists(fn):
        z = np.load(fn, allow_pickle=True)
        return {k: z[k] for k in z.files}
    C = Cube(sec)
    W = build_windows(C)
    np.savez_compressed(fn, **W)
    return W


def analyse(M, Ds, Ss, Bs, plane=False, rB=None):
    rB = M.fit(Bs, Ss, direct=True) if rB is None else rB
    rD = M.fit(Ds, Ss, x0=rB['x'], y0=rB['y'], plane=plane)
    off = M.C.jac @ np.array([rD['x'] - rB['x'], rD['y'] - rB['y']])
    cs = M.C.jac @ rD['cov'] @ M.C.jac.T
    c2T, aT = M.chi2_fixed(Ds, Ss, rB['x'], rB['y'])
    return dict(xT=rB['x'], yT=rB['y'], FT=rB['amp'], xD=rD['x'], yD=rD['y'], amp=rD['amp'], dE=float(off[0]), dN=float(off[1]),
                formal_cov=cs.tolist(), chi2=rD['chi2'], npix=rD['npix'], dchi2_target=c2T - rD['chi2'], bound=rD['bound'],
                depth_peak_ppm=rD['amp'] / rB['amp'] * 1e6), rB


def main(sec):
    C = Cube(sec, centroid_cut=False)   # geometry only (windows cached with the cut applied)
    W = get_windows(sec)
    M = Model(C)
    kind = W['kind']
    sets = {'deep': (np.where(kind == 'deep')[0], 1.0), 'shallow': (np.where(kind == 'shallow')[0], RATIO)}
    allidx = np.r_[sets['deep'][0], sets['shallow'][0]]
    allrel = np.r_[np.ones(len(sets['deep'][0])), RATIO * np.ones(len(sets['shallow'][0]))]
    out = dict(sector=sec, cam=C.cam, ccd=C.ccd, jac=C.jac, n_nb=len(M.nb), box=BOX)
    figs = []
    for name, (idx, rel) in list(sets.items()) + [('both', (allidx, allrel))]:
        Ds, Ss, Bs = stack(W, idx, rel if np.ndim(rel) else np.full(len(idx), rel))
        r, rB = analyse(M, Ds, Ss, Bs)
        rp, _ = analyse(M, Ds, Ss, Bs, plane=True, rB=rB)
        r['plane'] = dict(dE=rp['dE'], dN=rp['dN'])
        r['n'] = int(len(idx))
        out[name] = r
        figs.append((name, Ds, Ss, Bs, r))
        print(f"S{sec} {name:7s} n={len(idx):2d} depth_pk={r['depth_peak_ppm']:.0f}ppm off=({r['dE']:+.2f},{r['dN']:+.2f})\" "
              f"formal=({np.sqrt(r['formal_cov'][0][0]):.2f},{np.sqrt(r['formal_cov'][1][1]):.2f}) plane=({rp['dE']:+.2f},{rp['dN']:+.2f}) "
              f"dchi2T={r['dchi2_target']:.1f} chi2/npix={r['chi2']:.0f}/{r['npix']}", flush=True)
    # Gaussian PSF variant: widths fitted to the 'both' direct image
    Ds, Ss, Bs = stack(W, allidx, allrel)

    def obj(p):
        if min(p) < 0.3 or max(p) > 1.5:
            return 1e30
        return Model(C, kind='gauss', gsig=p).fit(Bs, Ss, direct=True)['chi2']
    gs = minimize(obj, [0.65, 0.65], method='Nelder-Mead', options=dict(xatol=3e-3, fatol=1e-2, maxiter=80)).x
    MG = Model(C, kind='gauss', gsig=gs)
    out['gauss_sigma'] = gs.tolist()
    for name, (idx, rel) in sets.items():
        Ds, Ss, Bs = stack(W, idx, np.full(len(idx), rel))
        r, _ = analyse(MG, Ds, Ss, Bs)
        out[name]['gauss'] = dict(dE=r['dE'], dN=r['dN'])
        print(f"S{sec} {name:7s} Gaussian({gs[0]:.2f},{gs[1]:.2f}) off=({r['dE']:+.2f},{r['dN']:+.2f})", flush=True)
    # per-eclipse offsets (direct position from the 'both' stack)
    rB = M.fit(Bs, Ss, direct=True)
    per = []
    for i in allidx:
        rel = 1.0 if kind[i] == 'deep' else RATIO
        D1, S1, B1 = stack(W, [i], [rel])
        try:
            r1, _ = analyse(M, D1, S1, B1, rB=rB)
        except Exception:
            continue
        per.append(dict(kind=str(kind[i]), tc=float(W['tc'][i]), dE=r1['dE'], dN=r1['dN'], amp=r1['amp'], bound=r1['bound'],
                        w=float((rel * W['sb'][i]) ** 2 / np.median(W['S'][i] ** 2))))
    out['per_eclipse'] = per
    to_json(out, f'{V3}/results/real_s{sec}{TAG}.json')
    fig, axs = plt.subplots(3, 4, figsize=(15, 11))
    for row, (name, Ds, Ss, Bs, r) in enumerate(figs):
        ax = axs[row]
        ax[0].imshow(np.log10(np.clip(Bs - np.median(Bs), 1, None)), origin='lower', cmap='gray_r'); ax[0].set_title(f'S{sec} {name} direct')
        v = np.nanmax(np.abs(Ds[M.sl])); ax[1].imshow(Ds, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v); ax[1].set_title('D')
        ax[2].imshow(Ds / Ss, origin='lower', cmap='RdBu_r', vmin=-8, vmax=8); ax[2].set_title('D/sigma')
        mod = r['amp'] * M.img(r['xD'], r['yD'])
        ax[3].imshow((Ds - mod - np.median((Ds - mod)[M.sl])) / Ss, origin='lower', cmap='RdBu_r', vmin=-4, vmax=4); ax[3].set_title('resid/sigma')
        for a in ax:
            a.plot(r['xT'], r['yT'], '+', color='lime', ms=14, mew=2); a.plot(r['xD'], r['yD'], 'x', color='m', ms=10, mew=2)
            for dx, dy, _ in M.nb:
                a.plot(r['xT'] + dx, r['yT'] + dy, 'o', mfc='none', mec='orange')
            a.set_xlim(C.ix - 5.5, C.ix + 5.5); a.set_ylim(C.iy - 5.5, C.iy + 5.5)
    plt.tight_layout(); plt.savefig(f'{V3}/plots/real_s{sec}{TAG}.png', dpi=60); plt.close()


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        main(s)
