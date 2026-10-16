"""TIC 61332759: per-pixel harmonic amplitude maps at the 1.4806-d variability
(fundamental and first harmonic = 0.7403 d), PRF-fit localisation of the
variable source vs Gaia positions, and a small-aperture / decorrelated search
for 0.7403-d dips on the target itself."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from astropy.timeseries import LombScargle, BoxLeastSquares
from reflib import *
from diffimg import fit_prf, pix_to_EN

TIC = 61332759
PV = 1.480598
NBR = 6711181657530577280
OUT = {}


def seg_design(t, nk=2):
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.8)])
    cols = []
    for s in np.unique(seg):
        w = seg == s
        x = (t - t[w].mean()) / max(1e-3, np.ptp(t[w]))
        for k in range(nk + 1):
            cols.append(np.where(w, x ** k, 0.0))
    return np.vstack(cols).T


for sec in (13, 27, 67, 94, 104):
    cut = load_cut(TIC, sec)
    t = cut['time']
    F = cut['flux']
    n, ny, nx = F.shape
    B = seg_design(t, 3)
    bkg = np.median(cut['bkg'].reshape(n, -1), axis=1)
    reg = [bkg - bkg.mean()]
    for pc in (cut['pc1'], cut['pc2']):
        if np.all(np.isfinite(pc)) and np.std(pc) > 0:
            reg.append(pc - pc.mean())
    B = np.hstack([B, np.vstack(reg).T])
    H = np.vstack([np.sin(2 * np.pi * t / PV), np.cos(2 * np.pi * t / PV), np.sin(4 * np.pi * t / PV), np.cos(4 * np.pi * t / PV)]).T
    X = np.hstack([B, H])
    Y = F.reshape(n, -1)
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    res = Y - X @ coef
    sig = np.std(res, axis=0)
    cov_diag = np.diag(np.linalg.pinv(X.T @ X))
    s1, c1, s2, c2 = [coef[B.shape[1] + i].reshape(ny, nx) for i in range(4)]
    e1 = (sig * np.sqrt(cov_diag[B.shape[1]])).reshape(ny, nx)
    A1 = np.hypot(s1, c1); A2 = np.hypot(s2, c2)
    g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
    tgt = g.iloc[0]
    nb = g[g.source_id == NBR].iloc[0]
    xt, yt = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    xn, yn = [float(v) for v in cut['wcs'].world_to_pixel_values(nb['ra_ep'], nb['dec_ep'])]
    yy_, xx_ = np.mgrid[:ny, :nx]
    central = np.hypot(xx_ - xt, yy_ - yt) < 5
    # signed projection onto the phase of the strongest central pixel (for each harmonic)
    j1 = np.unravel_index(np.argmax(np.where(central, A1, -1)), A1.shape); ph1 = np.arctan2(c1[j1], s1[j1])
    j2 = np.unravel_index(np.argmax(np.where(central, A2, -1)), A2.shape); ph2 = np.arctan2(c2[j2], s2[j2])
    S1 = s1 * np.cos(ph1) + c1 * np.sin(ph1)
    S2 = s2 * np.cos(ph2) + c2 * np.sin(ph2)
    prf = PRFModel(cut)
    direct = np.median(F, axis=0)
    derr = np.sqrt((0.002 * np.abs(direct)) ** 2 + 5 ** 2)
    rec = dict(peak_pix_A1=[int(j1[1]), int(j1[0])], peak_pix_A2=[int(j2[1]), int(j2[0])], target_pix=[xt, yt], nbr_pix=[xn, yn])
    for name, S in [('f1_1.48d', S1), ('f2_0.74d', S2)]:
        ee = e1 if name.startswith('f1') else e1  # same design -> similar per-pixel errors
        Sm = np.where(central, S, -np.inf)
        j = np.unravel_index(np.argmax(Sm), S.shape)
        fr = fit_prf(S, ee, prf, j[1], j[0], box=3)
        # positions relative to: fitted target (direct image) and fitted neighbour (direct image)
        ft = fit_prf(direct, derr, prf, xt, yt, box=2)
        fn = fit_prf(direct - ft['A'] * prf(ft['x'], ft['y']), derr, prf, xn, yn, box=1)
        oT = pix_to_EN(cut, fr['x'], fr['y'], ft['x'], ft['y'])
        oN = pix_to_EN(cut, fr['x'], fr['y'], fn['x'], fn['y'])
        # Gaia-based expected offset of neighbour from target
        gN = (nb['dE'], nb['dN'])
        rec[name] = dict(fit_pix=[fr['x'], fr['y']], amp_eps=fr['A'], offset_from_target_arcsec=list(oT), offset_from_neighbour_arcsec=list(oN),
                         gaia_nbr_offset_from_target=list(gN), redchi2=fr['chi2'] / fr['dof'])
        print(f"S{sec} {name}: PRF-fit source is ({oT[0]:+.1f},{oT[1]:+.1f})\" from the target's direct-image position and ({oN[0]:+.1f},{oN[1]:+.1f})\" from the neighbour's; "
              f"Gaia neighbour offset ({gN[0]:+.1f},{gN[1]:+.1f})\"; peak pixel (x,y)={j[1]},{j[0]}; target pix ({xt:.1f},{yt:.1f}) nbr pix ({xn:.1f},{yn:.1f})")
    # fractional amplitude at target vs neighbour pixel
    it = (int(round(yt)), int(round(xt))); inb = (int(round(yn)), int(round(xn)))
    rec['A2_frac_target_pix'] = float(A2[it] / direct[it]); rec['A2_frac_nbr_pix'] = float(A2[inb] / direct[inb])
    rec['A1_frac_target_pix'] = float(A1[it] / direct[it]); rec['A1_frac_nbr_pix'] = float(A1[inb] / direct[inb])
    print(f"   fractional semi-amp at target pixel: f1 {rec['A1_frac_target_pix']*1e6:.0f} ppm, f2 {rec['A2_frac_target_pix']*1e6:.0f} ppm; "
          f"at neighbour pixel: f1 {rec['A1_frac_nbr_pix']*1e6:.0f} ppm, f2 {rec['A2_frac_nbr_pix']*1e6:.0f} ppm")

    # --- target-only photometry: r=1 px aperture; also decorrelate against neighbour aperture LC
    tt, ft_, et_, ext = aperture_lc(TIC, sec, r=1.0)
    tn, fn_, en_, exn = aperture_lc(TIC, sec, r=1.0, center=(nb['ra_ep'], nb['dec_ep']))
    Bt = seg_design(tt, 3)
    # neighbour LC amplitude
    Hn = np.vstack([np.sin(2 * np.pi * tn / PV), np.cos(2 * np.pi * tn / PV), np.sin(4 * np.pi * tn / PV), np.cos(4 * np.pi * tn / PV)]).T
    cn, *_ = np.linalg.lstsq(np.hstack([Bt, Hn]), fn_, rcond=None)
    An1, An2 = np.hypot(*cn[-4:-2]), np.hypot(*cn[-2:])
    # target LC: remove baseline + neighbour leakage (regress on neighbour LC residuals) + own harmonics, then BLS at 0.7403 d
    nres = fn_ - Bt @ np.linalg.lstsq(Bt, fn_, rcond=None)[0]
    Xd = np.hstack([Bt, nres[:, None]])
    cd, *_ = np.linalg.lstsq(Xd, ft_, rcond=None)
    flat = ft_ - Xd @ cd + 1
    Ht = np.vstack([np.sin(2 * np.pi * tt / PV), np.cos(2 * np.pi * tt / PV), np.sin(4 * np.pi * tt / PV), np.cos(4 * np.pi * tt / PV)]).T
    ch, *_ = np.linalg.lstsq(np.hstack([Bt, Ht]), ft_, rcond=None)
    At1, At2 = np.hypot(*ch[-4:-2]), np.hypot(*ch[-2:])
    def dip_depth(tt_, ff_):
        P2, T02 = 0.7402990, 3148.78573
        sh = 1 - trapezoid(((tt_ - T02 + P2 / 2) % P2) - P2 / 2, 0.0, 1.0, 2.8 / 24, 0.7 / 24)
        X = np.hstack([seg_design(tt_, 3), sh[:, None]])
        b, *_ = np.linalg.lstsq(X, ff_, rcond=None)
        rr = ff_ - X @ b
        cov = np.linalg.pinv(X.T @ X) * np.var(rr, ddof=X.shape[1])
        return -b[-1] * 1e6, np.sqrt(cov[-1, -1]) * 1e6
    d_raw = dip_depth(tt, ft_)
    d_dec = dip_depth(tt, flat)
    d_nbr = dip_depth(tn, fn_)
    rec['target_r1'] = dict(A1_ppm=At1 * 1e6, A2_ppm=At2 * 1e6, nbr_A1_ppm=An1 * 1e6, nbr_A2_ppm=An2 * 1e6,
                            dip_depth_raw_ppm=list(d_raw), dip_depth_decorrelated_ppm=list(d_dec), nbr_dip_depth_ppm=list(d_nbr))
    qlp = load_lc(TIC, sec, 'QLP', 'SAP_FLUX')
    if qlp is not None:
        tq, fq = qlp[0], qlp[1]
        rec['QLP_SAP_dip_ppm'] = list(dip_depth(tq, fq))
    print(f"   r=1 target LC harmonics: f1 {At1*1e6:.0f} ppm, f2 {At2*1e6:.0f} ppm | neighbour r=1 LC: f1 {An1*1e6:.0f} ppm, f2 {An2*1e6:.0f} ppm")
    print(f"   0.7403-d dip depth at ephemeris: target r=1 raw {d_raw[0]:.0f}+-{d_raw[1]:.0f} ppm; decorrelated vs neighbour {d_dec[0]:.0f}+-{d_dec[1]:.0f} ppm; "
          f"neighbour r=1 {d_nbr[0]:.0f}+-{d_nbr[1]:.0f} ppm; QLP SAP {rec.get('QLP_SAP_dip_ppm',[np.nan,np.nan])[0]:.0f}+-{rec.get('QLP_SAP_dip_ppm',[np.nan,np.nan])[1]:.0f} ppm")
    OUT[f'S{sec}'] = rec
save_json(OUT, f'{BASE}/results/t2_var.json')
