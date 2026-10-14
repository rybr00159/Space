"""Real difference-image localisation per sector: deep, shallow and combined stacks, plus per-eclipse fits."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ana import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RATIO = 0.49   # shallow/deep depth ratio (aperture fits, s00)
TAG = os.environ.get('LOC_TAG', '')
PSFK = os.environ.get('LOC_PSF', 'prf')


def sets(ctx):
    k = ctx.kind
    deep = np.where(k == 'deep')[0]; sh = np.where(k == 'shallow')[0]
    return {'deep': (deep, np.ones(len(deep))), 'shallow': (sh, np.ones(len(sh))),
            'all': (np.r_[deep, sh], np.r_[np.ones(len(deep)), RATIO * np.ones(len(sh))])}


def gauss_width(ctx, B, sig):
    from scipy.optimize import minimize
    def obj(p):
        if min(p) < 0.3 or max(p) > 2.0:
            return 1e30
        ctx.set_gauss(*p)
        return ctx.fit_direct(B, sig)['chi2']
    r = minimize(obj, [0.75, 0.75], method='Nelder-Mead', options=dict(xatol=2e-3, fatol=1e-3, maxiter=60))
    ctx.set_gauss(*r.x)
    return r.x


def main(sec):
    ctx = Ctx(sec)
    out = dict(sector=sec, camera=int(ctx.z['camera']), ccd=int(ctx.z['ccd']), psf=PSFK, rfit=RFIT,
               jac=ctx.jac.tolist(), others=ctx.others, n_centroid_cut=int(ctx.z['nbad_centroid']))
    st = sets(ctx)
    if PSFK == 'gauss':
        Ds, ss, Bs = ctx.stack(*st['all'])
        out['gauss_sigma'] = gauss_width(ctx, Bs, ss).tolist()
    figs = []
    for name, (idx, an) in st.items():
        Ds, ss, Bs = ctx.stack(idx, an)
        r = ctx.analyse(Ds, ss, Bs)
        r['n_windows'] = int(len(idx))
        # residual rms of direct-image fit
        rd = ctx.fit_direct(Bs, ss)
        model_B = rd['amp'] * (ctx.psf(rd['x'], rd['y']) + ctx.others_img(rd['x'], rd['y'])) + rd['coef'][-1]
        m = ctx.mask
        r['direct_resid_rms_frac'] = float(np.sqrt(np.mean(((Bs - model_B) / Bs)[m] ** 2)))
        model_D = r['amp'] * ctx.psf(r['xD'], r['yD'])
        cD, _ = ctx.fit_fixed(Ds, ss, r['xD'], r['yD'])
        model_D = cD[0] * ctx.psf(r['xD'], r['yD']) + cD[1]
        r['diff_chi2_red'] = float(np.sum(((Ds - model_D) / ss)[m] ** 2) / (m.sum() - 4))
        out[name] = r
        print(f"S{sec} {name:7s} n={len(idx):2d} FT={r['FT']:.0f} amp={r['amp']:.1f} depth_pk={r['depth_peak']*1e6:.0f}ppm "
              f"off=({r['dE']:+.2f},{r['dN']:+.2f})\" formal=({np.sqrt(r['cov_sky'][0][0]):.2f},{np.sqrt(r['cov_sky'][1][1]):.2f}) "
              f"dchi2(target)={r['dchi2_target']:.1f} chi2red={r['diff_chi2_red']:.2f} Bresid={r['direct_resid_rms_frac']*100:.2f}%", flush=True)
        figs.append((name, Ds, ss, Bs, r, model_D))
    # per-eclipse fits (every real window alone), direct position from the 'all' stack
    Ds, ss, Bs = ctx.stack(*st['all'])
    rd_all = ctx.fit_direct(Bs, ss)
    per = []
    for i in np.where((ctx.kind == 'deep') | (ctx.kind == 'shallow'))[0]:
        an = 1.0 if ctx.kind[i] == 'deep' else RATIO
        D1, s1, B1 = ctx.stack([i], [an])
        try:
            r1 = ctx.analyse(D1, s1, B1, direct=rd_all)
        except Exception as e:
            continue
        per.append(dict(kind=str(ctx.kind[i]), tc=float(ctx.z['tc'][i]), epoch=int(ctx.z['epoch'][i]),
                        dE=r1['dE'], dN=r1['dN'], amp=r1['amp'], dchi2_target=r1['dchi2_target'], at_bound=r1['at_bound'],
                        w=float((an * ctx.z['sbar'][i]) ** 2 / np.median(ctx.z['sig'][i][ctx.mask] ** 2))))
    out['per_eclipse'] = per
    to_json(out, f'{L}/results/real_s{sec}_{PSFK}{TAG}.json')
    # figure
    fig, axs = plt.subplots(3, 4, figsize=(16, 12))
    for row, (name, Ds, ss, Bs, r, mD) in enumerate(figs):
        ax = axs[row]
        ax[0].imshow(np.log10(np.clip(Bs, 1, None)), origin='lower', cmap='gray_r'); ax[0].set_title(f'S{sec} {name}: direct (log)')
        v = np.nanmax(np.abs(Ds))
        ax[1].imshow(Ds, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v); ax[1].set_title('difference (deep-amp units)')
        ax[2].imshow(Ds / ss, origin='lower', cmap='RdBu_r', vmin=-8, vmax=8); ax[2].set_title('S/N')
        ax[3].imshow((Ds - mD) / ss, origin='lower', cmap='RdBu_r', vmin=-4, vmax=4); ax[3].set_title('resid/sigma (free PSF)')
        for a in ax:
            a.plot(r['xT'], r['yT'], '+', color='lime', ms=14, mew=2)
            a.plot(r['xD'], r['yD'], 'x', color='magenta', ms=10, mew=2)
            for dx, dy, rr, sid in ctx.others:
                a.plot(r['xT'] + dx, r['yT'] + dy, 'o', mfc='none', mec='orange', ms=6)
            a.set_xlim(r['xT'] - 5, r['xT'] + 5); a.set_ylim(r['yT'] - 5, r['yT'] + 5)
    plt.tight_layout(); plt.savefig(f'{L}/plots/real_s{sec}_{PSFK}{TAG}.png', dpi=70); plt.close()


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        main(s)
