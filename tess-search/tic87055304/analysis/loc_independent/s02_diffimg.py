"""Step 1 + 2(b): real difference images, per sector and per transit, with
 (i) one PSF at a free position (offset from target and neighbour, arcsec), and
 (ii) amplitudes at fixed positions (target+neighbour; all T<17.4 stars within 30"; all G<20 within 30"),
plus single-host chi2 for every Gaia star within 60" able to produce the depth.
Uncertainties: bootstrap over transits (B=2000) and per-transit scatter.

usage: python s02_diffimg.py [prf|gauss] [rfit] [reg_variant] [tag] [pointing]
Outputs: diff_<tag>.json, diff_<tag>.npz
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pipeline import *

psfkind = sys.argv[1] if len(sys.argv) > 1 else 'prf'
rfit = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
regv = sys.argv[3] if len(sys.argv) > 3 else 'local_r4'
tag = sys.argv[4] if len(sys.argv) > 4 else f'{psfkind}_r{rfit:g}_{regv}'
use_pointing = (len(sys.argv) > 5 and sys.argv[5] == 'pointing')
NBOOT = 2000
rng = np.random.default_rng(42)


def summarize(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    return dict(median=float(np.median(v)), mean=float(np.mean(v)), std=float(np.std(v, ddof=1)),
                p16=float(np.percentile(v, 16)), p84=float(np.percentile(v, 84)),
                p2=float(np.percentile(v, 2.5)), p97=float(np.percentile(v, 97.5)))


def main():
    res = {}
    arrays = {}
    for sec in SECTORS:
        t0 = time.time()
        ctx = SectorContext(sec, psfkind=psfkind, reg_variant=regv, rfit=rfit)
        S = ctx.S
        regs = None
        if use_pointing:
            dxs, dys, sc = pointing_series(ctx)
            regs = np.column_stack([dxs, dys])
        wins = ctx.real_windows()
        if os.environ.get('LOC_CLEAN', '0') == '1':
            # drop windows with a scattered-light change > 2 e/s/px (median of D in a 3.5-6 px annulus)
            ann = circle_mask((S.ny, S.nx), ctx.pT[0], ctx.pT[1], 6.0) & ~circle_mask((S.ny, S.nx), ctx.pT[0], ctx.pT[1], 3.5)
            wins = [w for w in wins if abs(np.median(ctx.diff(w, regressors=regs)[0][ann])) < 2.0]
        Ds, sigs, per = [], [], []
        for w in wins:
            D, sg, rms = ctx.diff(w, regressors=regs)
            Ds.append(D); sigs.append(sg)
            a = ctx.analyse(D, sg, free=True, extra=False)
            a.update(epoch=w['epoch'], tc=w['tc'], n_in=w['n_in'], n_l=w['n_l'], n_r=w['n_r'])
            if use_pointing or True:
                pass
            per.append(a)
        Ds = np.array(Ds); sigs = np.array(sigs)
        Dst, sgst = stack(Ds, sigs, mask=ctx.mask)
        full = ctx.analyse(Dst, sgst, free=True, extra=True)
        # chi2 / dof of the free fit as a noise-model check
        full['free_chi2_per_dof'] = full['free_chi2'] / (ctx.mask.sum() - 4)
        # bootstrap over transits
        nT = len(wins)
        boot = []
        for b in range(NBOOT):
            k = rng.integers(0, nT, nT)
            Db, sb = stack(Ds[k], sigs[k], mask=ctx.mask)
            a = ctx.analyse(Db, sb, free=True, extra=False)
            boot.append([a['offT_E'], a['offT_N'], a['s_par'], a['perp_arcsec'], a['f_N'], a['a_T'], a['a_N'],
                         a['dchi2_TminusN'], a['free_amp']])
        boot = np.array(boot)
        names = ['offT_E', 'offT_N', 's_par', 'perp', 'f_N', 'a_T', 'a_N', 'dchi2', 'free_amp']
        bsum = {n: summarize(boot[:, k]) for k, n in enumerate(names)}
        bcov = np.cov(boot[:, :2].T)
        # per-transit scatter
        pt = np.array([[p['offT_E'], p['offT_N'], p['s_par'], p['f_N']] for p in per])
        pts = dict(offT_E_mean=float(pt[:, 0].mean()), offT_N_mean=float(pt[:, 1].mean()),
                   offT_E_sem=float(pt[:, 0].std(ddof=1) / np.sqrt(nT)), offT_N_sem=float(pt[:, 1].std(ddof=1) / np.sqrt(nT)),
                   s_par_median=float(np.median(pt[:, 2])), s_par_mean=float(pt[:, 2].mean()),
                   s_par_sem=float(pt[:, 2].std(ddof=1) / np.sqrt(nT)),
                   f_N_mean=float(pt[:, 3].mean()), f_N_sem=float(pt[:, 3].std(ddof=1) / np.sqrt(nT)))
        res[str(sec)] = dict(n_transits=nT, epochs=[w['epoch'] for w in wins], stacked=full, bootstrap=bsum,
                             boot_cov_EN=bcov.tolist(), per_transit=per, per_transit_summary=pts,
                             target_xy=ctx.pT.tolist(), neigh_xy=ctx.pN.tolist(), jac=ctx.jac.tolist(),
                             neigh_EN=(ctx.jac @ ctx.vecTN).tolist(), reg_dxdy=[ctx.dx, ctx.dy],
                             psf=ctx.psf.__class__.__name__ if psfkind == 'prf' else ctx.psf.name)
        arrays[f's{sec}_D'] = Dst; arrays[f's{sec}_sig'] = sgst; arrays[f's{sec}_Dk'] = Ds; arrays[f's{sec}_sigk'] = sigs
        arrays[f's{sec}_boot'] = boot
        print(f'S{sec} {psfkind} nT={nT}: free offset from target (E,N)=({full["offT_E"]:+.2f},{full["offT_N"]:+.2f})" '
              f'boot sd=({bsum["offT_E"]["std"]:.2f},{bsum["offT_N"]["std"]:.2f}); s_par={full["s_par"]:.2f}+-{bsum["s_par"]["std"]:.2f}; '
              f'f_N={full["f_N"]:.2f}+-{bsum["f_N"]["std"]:.2f}; dchi2(T-N)={full["dchi2_TminusN"]:.1f}; '
              f'amp={full["free_amp"]:.1f} e/s; chi2/dof={full["free_chi2_per_dof"]:.2f}; {time.time()-t0:.0f}s', flush=True)
    to_json(res, f'{HERE}/diff_{tag}.json')
    np.savez_compressed(f'{HERE}/diff_{tag}.npz', **arrays)


if __name__ == '__main__':
    main()
