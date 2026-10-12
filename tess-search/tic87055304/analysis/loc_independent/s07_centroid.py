"""Model-light check: flux-weighted centroid shift between difference image and direct image.

For a circular aperture (radius 1.5, 2.0, 2.5 px) centred on the registered target position:
   c_D = sum(x D)/sum(D)         (difference image, stacked)
   c_M = sum(x M)/sum(M)         (mean out-of-transit image)
   shift = J (c_D - c_M)  [arcsec], projected on the target->neighbour axis in units of the separation.
No PSF fit is involved in the measurement; a constant (median of the 3.5-6 px annulus) is subtracted from D.  Expectations: target host -> shift ~ -(blend fraction)
(i.e. slightly AWAY from the neighbour, about -0.1..-0.2); neighbour host -> ~ +0.8.
Predictions are also computed with the TESS PRF and the Gaussian PSF (model D for a host at the
target or the neighbour, plus the real M), and uncertainties by bootstrap over transits.
Outputs: centroid_shift.json
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pipeline import *

rng = np.random.default_rng(3)
arr = np.load(f'{HERE}/diff_prf_main.npz')
out = {}
for sec in SECTORS:
    S = Sector(sec)
    res = {}
    for psfkind, regv in [('prf', 'local_r4'), ('gauss', 'local_r4_ratiofix')]:
        ctx = SectorContext(sec, psfkind=psfkind, reg_variant=regv, S=S)
        if psfkind == 'prf':
            M = S.flux[S.good].mean(0)
            Dk = arr[f's{sec}_Dk']; Sk = arr[f's{sec}_sigk']
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]
        vT = ctx.jac @ ctx.vecTN
        for rap in (1.5, 2.0, 2.5):
            ap = circle_mask(M.shape, ctx.pT[0], ctx.pT[1], rap)
            def cen(img):
                w = img[ap]
                return np.array([(xx[ap] * w).sum() / w.sum(), (yy[ap] * w).sum() / w.sum()])
            cM = cen(M)
            ann = circle_mask(M.shape, ctx.pT[0], ctx.pT[1], 6.0) & ~circle_mask(M.shape, ctx.pT[0], ctx.pT[1], 3.5)
            def s_of(D, bsub=True):
                if bsub:
                    D = D - np.median(D[ann])   # remove constant (scattered-light) offset of the difference image
                sh = ctx.jac @ (cen(D) - cM)
                return float(sh @ vT / (vT @ vT)), sh
            key = f'r{rap}'
            if psfkind == 'prf':
                D, _ = stack(Dk, Sk, mask=ctx.mask)
                s_obs, sh = s_of(D)
                bs = []
                for b in range(2000):
                    k = rng.integers(0, len(Dk), len(Dk))
                    bs.append(s_of(stack(Dk[k], Sk[k], mask=ctx.mask)[0])[0])
                bs = np.array(bs)
                res[key] = dict(s_obs=s_obs, shift_EN=sh.tolist(), s_boot_sd=float(np.std(bs)),
                                s_boot_p16_p84=np.percentile(bs, [16, 84]).tolist())
            sT, _ = s_of(ctx.psf(*ctx.pT), bsub=False)
            sN, _ = s_of(ctx.psf(*ctx.pN), bsub=False)
            res[key][f'pred_T_{psfkind}'] = sT
            res[key][f'pred_N_{psfkind}'] = sN
    out[str(sec)] = res
    for k, v in res.items():
        print(f"S{sec} {k}: observed s={v['s_obs']:+.3f}+-{v['s_boot_sd']:.3f}   predicted target-host {v['pred_T_prf']:+.3f} (PRF) {v['pred_T_gauss']:+.3f} (Gauss);"
              f"  neighbour-host {v['pred_N_prf']:+.3f} (PRF) {v['pred_N_gauss']:+.3f} (Gauss)")
to_json(out, f'{HERE}/centroid_shift.json')
