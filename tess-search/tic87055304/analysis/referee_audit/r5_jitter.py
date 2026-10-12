"""Check 3: is the pointing-jitter leverage argument right, and is jitter a viable shared systematic?

(a) Leverage. A free-position PSF fit to D = A*M(x - x_src) + F_pair * s . grad M  (an in-transit image shift s,
    D = baseline - in-transit) returns x_src - (F_pair/A) s. Verified numerically below with the SPOC PRF: inject
    s = 1 mpx into a target-hosted D, fit the position (variable-projection least squares, r <= 3 px), read off the shift.
(b) In-transit pointing shift measured from the difference image itself (no PSF model): in pixels > 3 px from the
    target, D_far = s_x dO/dx + s_y dO/dy + plane (O = mean OOT image, gradients by central differences; the
    gradient-regression scale is calibrated by a known spline shift of O). One (s_x, s_y) per real transit and per fake
    window (phases 0.20-0.80). Projected on u (target -> neighbour).
(c) Spurious s_par predicted per transit = -(F_pair/A) s_u / sep; stacked with the transit weights.
"""
import numpy as np
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares
from rcommon import Cube, SECTORS, OUT, P, T0, T14, diff_image, dump
from r3_pixel_ratio import prf_model, PHASES

ZP = 0.86 * 15000.0


def grad_fit(D, O, mask):
    gy, gx = np.gradient(O)
    yy, xx = np.mgrid[0:O.shape[0], 0:O.shape[1]]
    A = np.column_stack([gx[mask], gy[mask], np.ones(mask.sum()), xx[mask], yy[mask]])
    coef, *_ = np.linalg.lstsq(A, D[mask], rcond=None)
    return coef[:2]


def main():
    res = {}
    for s in SECTORS:
        c = Cube(s)
        yy, xx = np.mgrid[0:c.ny, 0:c.nx]
        rr = np.hypot(xx - c.xt, yy - c.yt)
        farm = (rr > 3) & (xx > 0) & (yy > 0) & (xx < c.nx - 1) & (yy < c.ny - 1)
        tp = (int(round(c.xt)), int(round(c.yt)))
        prf = prf_model(c)
        g = c.g
        F_T = ZP * 10 ** (-0.4 * (g.tmag.values[c.it] - 10)); F_N = ZP * 10 ** (-0.4 * (g.tmag.values[c.ineb] - 10))
        F_pair = F_T + F_N
        wins = []
        for sh in [0.0] + list(PHASES):
            for n in c.epochs():
                w = c.window(T0 + (n + sh) * P)
                if not w['ok']:
                    continue
                D, sD, O = diff_image(c.t, c.f, w)
                Ob = O - np.percentile(O, 10)
                wins.append(dict(real=(sh == 0.0), n=int(n), D=D, O=Ob, wt=1 / sD[tp[1], tp[0]] ** 2,
                                 floor=float(np.median(D[rr > 5]))))
        Omean = np.mean([w['O'] for w in wins if w['real']], 0)
        # calibrate the gradient regression: shift O by a known 0.01 px in x and y
        cal = []
        for sx, sy in ((0.01, 0.0), (0.0, 0.01)):
            Osh = ndshift(Omean, (sy, sx), order=3, mode='nearest')
            cal.append(grad_fit(Omean - Osh, Omean, farm) / 0.01)   # = D for an in-transit shift +s
        cal = np.array(cal)   # rows: injected x / y; columns: measured x / y. D = O(x) - O(x - s) => +s
        Ainv = np.linalg.inv(cal.T)
        A_amp = float(np.average([w['D'][rr <= 3].sum() - w['floor'] * (rr <= 3).sum() for w in wins if w['real']],
                                 weights=[w['wt'] for w in wins if w['real']]))
        # (a) numerical leverage with the PRF
        MT = prf(c.xt, c.yt); mask = rr <= 3
        gyM, gxM = np.gradient(MT)
        def fit_pos(img):
            def resid(p):
                M = prf(p[0], p[1])[mask]
                Am = np.column_stack([M, np.ones(mask.sum())])
                co, *_ = np.linalg.lstsq(Am, img[mask], rcond=None)
                return Am @ co - img[mask]
            return least_squares(resid, [c.xt, c.yt], diff_step=1e-4, x_scale=0.05).x
        s_test = 0.001 * c.u     # 1 mpx toward the neighbour
        img0 = A_amp * MT
        Mshift = prf(c.xt + s_test[0], c.yt + s_test[1])
        # in-transit image I(x) = O(x - s) - A M_T  =>  D = O - I = A M_T - F_pair (M_T(shifted by +s) - M_T)
        img1 = img0 - F_pair * (Mshift - MT)
        p0, p1 = fit_pos(img0), fit_pos(img1)
        lev_num = float(((p1 - p0) @ c.u) / 0.001)    # expected -(F_pair/A): source appears displaced opposite to s
        # (b) per-window in-transit pointing shift
        out = []
        for w in wins:
            m = grad_fit(w['D'] - w['floor'], Omean, farm)
            sxy = Ainv @ m
            su = float(sxy @ c.u); sp = float(sxy @ c.uperp)
            out.append(dict(real=w['real'], n=w['n'], s_u_mpx=1e3 * su, s_perp_mpx=1e3 * sp, wt=w['wt'],
                            spurious_s_par=float(-(F_pair / A_amp) * su / c.sep_px)))
        R = [o for o in out if o['real']]; Fk = [o for o in out if not o['real']]
        wts = np.array([o['wt'] for o in R])
        su_stack = float(np.average([o['s_u_mpx'] for o in R], weights=wts))
        mad = lambda v: 1.4826 * np.median(np.abs(np.asarray(v) - np.median(v)))
        r = dict(F_pair_e_s=F_pair, A_e_s=A_amp, leverage_analytic=F_pair / A_amp, leverage_numeric_PRF=lev_num,
                 sep_px=c.sep_px, s_u_needed_mpx=float(-1e3 * c.sep_px * A_amp / F_pair),
                 gradient_calibration_matrix=cal, real=R,
                 real_s_u_stack_mpx=su_stack, real_spurious_s_par_stack=float(-(F_pair / A_amp) * su_stack / 1e3 / c.sep_px),
                 fake_s_u_std_mpx=float(np.std([o['s_u_mpx'] for o in Fk])), fake_s_u_mad_mpx=float(mad([o['s_u_mpx'] for o in Fk])),
                 n_fake=len(Fk))
        res[str(s)] = r
        print(f"S{s}: F_pair {F_pair:.0f} e/s, A {A_amp:.1f} e/s, leverage analytic {F_pair/A_amp:.0f} numeric(PRF fit) {lev_num:.0f}; "
              f"needed s_u {r['s_u_needed_mpx']:.2f} mpx; cal diag {np.round(np.diag(cal),2)}; real s_u per transit "
              f"{np.round([o['s_u_mpx'] for o in R],2)} mpx (stack {su_stack:.2f}); predicted spurious s_par per transit "
              f"{np.round([o['spurious_s_par'] for o in R],2)} (stack {r['real_spurious_s_par_stack']:.2f}); "
              f"fake s_u std {r['fake_s_u_std_mpx']:.2f} mad {r['fake_s_u_mad_mpx']:.2f} mpx (n={len(Fk)})", flush=True)
    dump(res, f'{OUT}/r5_jitter.json')


if __name__ == '__main__':
    main()
