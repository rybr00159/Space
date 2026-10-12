"""Check 5 (decisive, PSF-model-free core): relative transit depth in single pixels on the neighbour side vs the
target's own pixel and the far side.

Definitions (per sector, pixel frame; orientation of u verified in r1_wcs_check.py):
  u  = unit vector target -> neighbour; s(p) = (p - x_t).u
  a  = pixel nearest x_t + 1 px*u   (neighbour side)       b = pixel nearest x_t - 1 px*u  (far side)
  t  = pixel containing the target  n = pixel containing the neighbour (only S93, S104 differ from t)
  Ga = pixels with r <= 2 px and s >= +0.5 px            Gb = pixels with r <= 2 px and s <= -0.5 px
Per transit: per-pixel linear baseline from both OOT sides (0.75-2.0 T14 from mid-transit), in-transit = central
75 % of T14; D = baseline - in-transit (e/s); stacked with one scalar weight per transit (1/var of D in pixel t);
the far-field D model (plane a + b x + c y fitted to pixels > 5 px from the target; scattered-light residual)
is subtracted. Symmetric cleaning: windows (real and fake alike) whose far-field D median deviates from the
sector's fake-window median by > 10 robust sigma are dropped (env R3_CLEAN=none keeps all; R3_FLOOR=median).
O = weighted mean OOT image of the same windows minus background B = 10th percentile of O (variants 5th/25th).

Relative depth delta(p) = D(p)/O(p). If the TARGET hosts a transit of depth delta_T:
    delta(p) = delta_T * fT(p),  fT(p) = fraction of the light in p that comes from the target, 0 <= fT <= 1,
so for any pair Q(x,y) = delta(x)/delta(y) = fT(x)/fT(y) <= 1/fT(y). No PSF model enters this inequality.
If the NEIGHBOUR hosts: Q = fN(x)/fN(y), large when x is nearer the neighbour.

To turn the inequality into a number one needs a LOWER bound on fT(y). Conservative envelope: every non-target
Gaia star (T < 18) contributes the MAXIMUM, over the SPOC PRF and Gaussians of sigma 0.5/0.65/0.8 px (fitted TESS widths are 0.47-0.64 px), of its pixel
fraction at its true position, x1.25, with fluxes from T mag scaled to the observed stamp flux, plus 20 e/s of
background error. kappa = 1/fT_min(y); H_target => L = delta(x) - kappa*delta(y) <= 0.
Noise: 4000 stacks of fake transits at phases 0.20-0.80 (same estimator, same number of windows, same weights);
z = L / std(L_null); also z with a MAD-based sigma and the empirical fraction of null stacks >= observed. Bootstrap over real transits gives a ratio CI.

PSF-model predictions of Q (context only): Q_T = [PSF_T(x)/PSF_T(y)] / [O(x)/O(y)], Q_N likewise with PSF_N.
"""
import numpy as np
from scipy.special import erf
from rcommon import Cube, SECTORS, OUT, P, T0, T14, diff_image, dump, TID

RNG = np.random.default_rng(12345)
PHASES = np.round(np.arange(0.20, 0.801, 0.04), 3)
NBOOT = 4000
ZP_SCALE = 0.86      # observed/predicted stamp flux with a 15000 e/s (T=10) zero point (see main(), printed)


def gauss_psf(shape, x, y, s):
    xe = np.arange(shape[1] + 1) - 0.5; ye = np.arange(shape[0] + 1) - 0.5
    cx = 0.5 * np.diff(erf((xe - x) / (np.sqrt(2) * s))); cy = 0.5 * np.diff(erf((ye - y) / (np.sqrt(2) * s)))
    return np.outer(cy, cx)


def prf_model(c):
    import PRF
    p = PRF.TESS_PRF(c.cam, c.ccd, c.s, c.col0 + c.nx / 2, c.row0 + c.ny / 2, localdatadir=f'{OUT}/prf_files/')
    shp = (c.ny, c.nx)
    def f(x, y):
        if x < -5 or y < -5 or x > c.nx + 4 or y > c.ny + 4:
            return np.zeros(shp)
        return p.locate(x, y, shp)
    return f


import os
CLEAN = os.environ.get('R3_CLEAN', 'mad10')     # 'mad10': drop windows (real AND fake) whose far-field D floor deviates
                                               # from the sector's fake-window median by > 10 robust sigma; 'none'
FLOOR = os.environ.get('R3_FLOOR', 'plane')     # far-field D model subtracted: 'plane' (a + b x + c y) or 'median'


def collect(c, shift=0.0):
    out = []
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    far = np.hypot(xx - c.xt, yy - c.yt) > 5
    A = np.column_stack([np.ones(far.sum()), (xx - c.xt)[far], (yy - c.yt)[far]])
    for n in c.epochs():
        tc = T0 + (n + shift) * P
        w = c.window(tc)
        if not w['ok']:
            continue
        D, sD, O = diff_image(c.t, c.f, w)
        if FLOOR == 'plane':
            coef, *_ = np.linalg.lstsq(A, D[far], rcond=None)
            fl_img = coef[0] + coef[1] * (xx - c.xt) + coef[2] * (yy - c.yt)
        else:
            fl_img = np.full(D.shape, np.median(D[far]))
        out.append(dict(n=int(n), tc=tc, D=D, sD=sD, O=O, floor=fl_img, floor_val=float(np.median(D[far]))))
    return out


def clean(real, fakes):
    """Symmetric scattered-light cut (same rule for real and fake windows)."""
    ff = np.array([w['floor_val'] for w in fakes])
    med = np.median(ff); mad = 1.4826 * np.median(np.abs(ff - med))
    if CLEAN == 'none':
        return real, fakes, []
    keep = lambda w: abs(w['floor_val'] - med) <= 10 * mad
    return [w for w in real if keep(w)], [w for w in fakes if keep(w)], [w['n'] for w in real if not keep(w)]


def stack(wins, idx, bkg_pct=10):
    ws = np.array([wins[k]['wt'] for k in idx])
    D = np.tensordot(ws, np.array([wins[k]['D'] - wins[k]['floor'] for k in idx]), 1) / ws.sum()
    O = np.tensordot(ws, np.array([wins[k]['O'] for k in idx]), 1) / ws.sum()
    B = np.percentile(O, bkg_pct)
    return D, O - B, B


def group_vals(img, grp):
    return float(sum(img[y, x] for x, y in grp))


def main():
    res = {}
    nulls = {}
    for s in SECTORS:
        c = Cube(s)
        shp = (c.ny, c.nx)
        yy, xx = np.mgrid[0:c.ny, 0:c.nx]
        rr = np.hypot(xx - c.xt, yy - c.yt)
        far = rr > 5
        tp = (int(round(c.xt)), int(round(c.yt)))
        real = collect(c)
        allfake = [w for ph in PHASES for w in collect(c, shift=ph)]
        real, allfake, dropped_real = clean(real, allfake)
        for w in real + allfake:
            w['wt'] = 1.0 / w['sD'][tp[1], tp[0]] ** 2
        nr = len(real)
        # pixel groups
        a = (int(round(c.xt + c.u[0])), int(round(c.yt + c.u[1])))
        b = (int(round(c.xt - c.u[0])), int(round(c.yt - c.u[1])))
        npx = (int(round(c.xn)), int(round(c.yn)))
        sproj = (xx - c.xt) * c.u[0] + (yy - c.yt) * c.u[1]
        Ga = [(int(x), int(y)) for y, x in zip(*np.where((rr <= 2.0) & (sproj >= 0.5)))]
        Gb = [(int(x), int(y)) for y, x in zip(*np.where((rr <= 2.0) & (sproj <= -0.5)))]
        groups = {'a': [a], 'b': [b], 't': [tp], 'Ga': Ga, 'Gb': Gb}
        if npx != tp:
            groups['n'] = [npx]
        pairs = [('a', 'b'), ('a', 't'), ('t', 'b'), ('Ga', 'Gb'), ('Ga', 't')] + ([('n', 't')] if 'n' in groups else [])
        # models
        prf = prf_model(c)
        g = c.g
        fl = ZP_SCALE * 15000.0 * 10 ** (-0.4 * (g.tmag.values - 10.0))
        ok = np.isfinite(g.x.values) & (np.hypot(g.x.values - c.xt, g.y.values - c.yt) < 8) & (g.tmag.values < 18)
        models = {'prf': prf, 'g0.5': lambda x, y: gauss_psf(shp, x, y, 0.5),
                  'g0.65': lambda x, y: gauss_psf(shp, x, y, 0.65), 'g0.8': lambda x, y: gauss_psf(shp, x, y, 0.8)}
        env_contam = np.zeros(shp)       # conservative envelope of non-target light per pixel
        for j in np.where(ok)[0]:
            if g.source_id.values[j] == TID:
                continue
            env = np.max([m(g.x.values[j], g.y.values[j]) for m in models.values()], axis=0)
            env_contam += 1.25 * fl[j] * env
        env_contam += 20.0
        PT = {k: m(c.xt, c.yt) for k, m in models.items()}
        PN = {k: m(c.xn, c.yn) for k, m in models.items()}
        # self-consistent scene models for the 'most target-friendly' prediction: all Gaia stars (T<18, <8 px) with
        # each PSF model, global shifts of 0 / +-0.1 px along u and perpendicular, contaminant fluxes x0.7/1/1.3,
        # residual background -20/0/+20 e/s. fT(p) = F_T M_T(p) / (scene(p) + dB); Q_T = fT(x)/fT(y), Q_N likewise.
        scenes = []
        jj = [j for j in np.where(ok)[0]]
        shifts = [(0, 0)] + [tuple(sg * 0.1 * v) for v in (c.u, c.uperp) for sg in (1, -1)]
        for k, m in models.items():
            for sh in shifts:
                comps = {j: m(g.x.values[j] + sh[0], g.y.values[j] + sh[1]) * fl[j] for j in jj}
                T_ = comps[c.it]; N_ = comps[c.ineb]
                rest = sum(v for j, v in comps.items() if j not in (c.it, c.ineb))
                for cs in (0.7, 1.0, 1.3):
                    for dB in (-20.0, 0.0, 20.0):
                        tot = T_ + cs * N_ + cs * rest + dB   # neighbour and other stars scaled together
                        scenes.append(dict(model=k, shift=sh, cscale=cs, dB=dB, fT=T_ / tot, fN=N_ * cs / tot))
        # real stack
        D, O, B = stack(real, range(nr))
        # null stacks
        nullD, nullO = [], []
        for _ in range(NBOOT):
            idx = RNG.choice(len(allfake), nr, replace=False)
            Df, Of, _ = stack(allfake, idx)
            nullD.append(Df); nullO.append(Of)
        nullD = np.array(nullD); nullO = np.array(nullO)
        bootD, bootO = [], []
        for _ in range(NBOOT):
            idx = RNG.integers(0, nr, nr)
            Db, Ob, _ = stack(real, idx)
            bootD.append(Db); bootO.append(Ob)
        bootD = np.array(bootD); bootO = np.array(bootO)
        # per-transit stacks
        ptD = [stack(real, [k])[0] for k in range(nr)]
        ptO = [stack(real, [k])[1] for k in range(nr)]
        # D/O map near the target (for the figure)
        r = dict(n_transits=nr, epochs=[w['n'] for w in real], dropped_real_epochs=dropped_real,
                 n_fake_windows=len(allfake), B=B,
                 target_px=[c.xt, c.yt], neighbour_px=[c.xn, c.yn], u=c.u.tolist(),
                 groups={k: v for k, v in groups.items()}, total_dF_r3=float(D[rr <= 3].sum()),
                 pairs={})
        for x, y in pairs:
            gx, gy = groups[x], groups[y]
            dFx, dFy = group_vals(D, gx), group_vals(D, gy)
            Ox, Oy = group_vals(O, gx), group_vals(O, gy)
            dx_, dy_ = dFx / Ox, dFy / Oy
            Q = dx_ / dy_
            Cy = group_vals(env_contam, gy)
            fT_min = 1.0 - Cy / Oy
            kappa = 1.0 / fT_min if fT_min > 0.02 else np.inf
            L = dx_ - kappa * dy_ if np.isfinite(kappa) else np.nan
            L1 = dx_ - dy_
            # null
            ndx = np.array([group_vals(nullD[k], gx) / Ox for k in range(NBOOT)])
            ndy = np.array([group_vals(nullD[k], gy) / Oy for k in range(NBOOT)])
            Ln = ndx - (kappa if np.isfinite(kappa) else 0) * ndy
            L1n = ndx - ndy
            # Q interval from null noise added to the observed values
            Qn = (dx_ + ndx) / (dy_ + ndy)
            bq = np.array([(group_vals(bootD[k], gx) / group_vals(bootO[k], gx)) /
                           (group_vals(bootD[k], gy) / group_vals(bootO[k], gy)) for k in range(NBOOT)])
            pt = [(group_vals(ptD[k], gx) / group_vals(ptO[k], gx)) / (group_vals(ptD[k], gy) / group_vals(ptO[k], gy))
                  for k in range(nr)]
            # model predictions (context)
            pred = {}
            for k in models:
                RT = group_vals(PT[k], gx) / group_vals(PT[k], gy)
                RN = group_vals(PN[k], gx) / group_vals(PN[k], gy)
                pred[k] = dict(Q_T=RT / (Ox / Oy), Q_N=RN / (Ox / Oy))   # simple: observed O as denominator
            # self-consistent scene predictions (ranges over PSF model, shift, flux scale, background error)
            QTs, QNs = [], []
            for sc in scenes:
                # group fraction: target light estimated as fT_scene(p) * O(p), summed over the group, / sum O
                fTx = sum(sc['fT'][yy_, xx_] * O[yy_, xx_] for xx_, yy_ in gx) / Ox
                fTy = sum(sc['fT'][yy_, xx_] * O[yy_, xx_] for xx_, yy_ in gy) / Oy
                fNx = sum(sc['fN'][yy_, xx_] * O[yy_, xx_] for xx_, yy_ in gx) / Ox
                fNy = sum(sc['fN'][yy_, xx_] * O[yy_, xx_] for xx_, yy_ in gy) / Oy
                QTs.append(fTx / fTy); QNs.append(fNx / fNy)
            QTs, QNs = np.array(QTs), np.array(QNs)
            pred['scene'] = dict(Q_T_range=[float(QTs.min()), float(np.median(QTs)), float(QTs.max())],
                                 Q_N_range=[float(QNs.min()), float(np.median(QNs)), float(QNs.max())])
            QTmax = float(QTs.max()); QNmin = float(QNs.min()); QNmax = float(QNs.max())
            LT = dx_ - QTmax * dy_; LTn = ndx - QTmax * ndy
            mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
            pr_extra = dict(Q_T_max_models=QTmax, Q_N_range_models=[QNmin, QNmax],
                            LT_ppt=1e3 * LT, sigma_LT_null_ppt=1e3 * LTn.std(), z_T_models=float(LT / LTn.std()),
                            z_T_models_mad=float(LT / mad(LTn)), p_emp_T_models=float(np.mean(LTn >= LT)),
                            z_L_mad=float(L / mad(Ln)) if np.isfinite(L) else None,
                            p_emp_L=float(np.mean(Ln >= L)) if np.isfinite(L) else None,
                            # consistency with H_N: (delta_x - Q_N*delta_y)/sigma for the min and max model Q_N
                            z_N_models=[float((dx_ - q * dy_) / (ndx - q * ndy).std()) for q in (QNmin, QNmax)])
            pr = dict(**pr_extra, x=gx, y=gy, dF_x=dFx, dF_y=dFy, O_x=Ox, O_y=Oy, delta_x_ppt=1e3 * dx_, delta_y_ppt=1e3 * dy_,
                      Q=Q, Q_null_noise_16_50_84=np.percentile(Qn, [16, 50, 84]),
                      Q_null_noise_2p5_97p5=np.percentile(Qn, [2.5, 97.5]),
                      Q_boot_16_50_84=np.percentile(bq, [16, 50, 84]),
                      sigma_delta_x_null_ppt=1e3 * ndx.std(), sigma_delta_y_null_ppt=1e3 * ndy.std(),
                      corr_null=float(np.corrcoef(ndx, ndy)[0, 1]),
                      contam_env_y=Cy, fT_min_y=fT_min, kappa=kappa,
                      L_ppt=1e3 * L if np.isfinite(L) else None, z_L=float(L / Ln.std()) if np.isfinite(L) else None,
                      L1_ppt=1e3 * L1, z_L1=float(L1 / L1n.std()),
                      per_transit_Q=pt, n_transits_Q_gt_1=int(np.sum(np.array(pt) > 1)),
                      n_transits_Q_gt_kappa=int(np.sum(np.array(pt) > kappa)) if np.isfinite(kappa) else None,
                      model_context=pred)
            r['pairs'][f'{x}/{y}'] = pr
            nulls[f'S{s}_{x}_{y}_zT_null'] = LTn / LTn.std()
            print(f"S{s} {x+'/'+y:6s} x={gx if len(gx)==1 else len(gx)} y={gy if len(gy)==1 else len(gy)}"
                  f" dF={dFx:.2f}/{dFy:.2f} O={Ox:.0f}/{Oy:.0f} delta={1e3*dx_:.2f}/{1e3*dy_:.2f}ppt"
                  f" Q={Q:.2f} [{pr['Q_null_noise_2p5_97p5'][0]:.2f},{pr['Q_null_noise_2p5_97p5'][1]:.2f}]"
                  f" fTmin(y)={fT_min:.2f} kappa={kappa:.2f} z_L={pr['z_L']} z_L1={pr['z_L1']:.1f}"
                  f" | scene Q_T {np.round(pred['scene']['Q_T_range'],2)} Q_N {np.round(pred['scene']['Q_N_range'],2)}"
                  f" | z_T(models)={pr['z_T_models']:.1f} (mad {pr['z_T_models_mad']:.1f}) z_N={np.round(pr['z_N_models'],1)}"
                  f" | per-transit Q>1: {pr['n_transits_Q_gt_1']}/{nr}", flush=True)
        # keep maps for the figure
        r['maps'] = dict(D=D, O=O, delta=np.where(O > 30, D / np.where(O > 30, O, 1), np.nan),
                         sigD=nullD.std(0), env_contam=env_contam)
        res[str(s)] = r
    dump(res, f'{OUT}/r3_pixel_ratio_{CLEAN}_{FLOOR}.json')
    np.savez_compressed(f'{OUT}/r3_nulls_{CLEAN}_{FLOOR}.npz', **nulls)


if __name__ == '__main__':
    main()
