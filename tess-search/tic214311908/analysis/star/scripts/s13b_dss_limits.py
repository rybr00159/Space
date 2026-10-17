"""Archival-plate test for a background star now hidden behind TIC 214311908 (proper-motion 'look-back' test).
The M3V moves 142 mas/yr; on the SERC-J (1975.3), SERC-I (1979.4), SERC-QV (1987.7) and AAO-SES R (1992.6) DSS scans
it lay 7.3/6.7/5.5/4.8" from its 2026.3 (TESS S104) position and 6.3/5.7/4.6/3.8" from its 2019.4 (S12) position, so a
background star at those positions would appear as a separate source.
Per plate (8' x 8' DSS2 scans; Gaia DR3 sources within 330" propagated to the plate epoch, data/gaia_330.csv):
 1. Astrometric alignment: affine correction from centroids of isolated G 13.5-17 stars.
 2. PSF: Moffat fitted to isolated G 16-19.5 stars (median FWHM, beta).
 3. 'Scene' photometry: in a box of radius 14" around a star, all Gaia sources (G < 21.5) are fitted simultaneously
    with the PSF at their (aligned) positions + a plane; global shift free. Magnitude calibration
    mag = poly2(log10 amp) against Gaia BP (J, QV plates) or RP (I, R plates) for ~300 field stars.
 4. Target box (radius 18"): the same scene model, the M3V with its own free Moffat (saturated photographic image),
    plus a TEST source at the 2019.4 or 2026.3 position; its amplitude -> magnitude or upper limit.
 5. Empirical noise of the test amplitude: same fit with the test source at 100 random positions in the field that are
    > 4" from any Gaia source; detection limit from injection-recovery of fake stars at the 2026.3 position
    (8 per magnitude step, 3-sigma detection in >= 90%).
Output: out/dss_limits.json, figs/dss_resid.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.optimize import least_squares
B = 'tic214311908/analysis/star'
rng = np.random.default_rng(7)
g = pd.read_csv(f'{B}/data/gaia_330.csv'); g = g[g.phot_g_mean_mag < 21.5].reset_index(drop=True)
M3V = 5925373272850784768
im3v = int(np.where(g.source_id == M3V)[0][0]); r = g.iloc[im3v]
imira = int(np.where(g.source_id == 5925373268528684288)[0][0])
def pos(ra, dec, pmra, pmdec, ep):
    dt = ep - 2016.0
    return ra + np.nan_to_num(pmra) * dt / 3.6e6 / np.cos(np.radians(dec)), dec + np.nan_to_num(pmdec) * dt / 3.6e6
def moffat(xx, yy, x0, y0, fwhm, beta):
    al = fwhm / (2 * np.sqrt(2 ** (1 / beta) - 1))
    return (1 + ((xx - x0) ** 2 + (yy - y0) ** 2) / al ** 2) ** (-beta) * (beta - 1) / (np.pi * al ** 2)
plates = {'SERC-J 1975': ('big_poss2ukstu_blue.fits', 'BP'), 'SERC-I 1979': ('big_poss2ukstu_ir.fits', 'RP'),
          'SERC-QV 1987': ('big_phase2_gsc1.fits', 'BP'), 'AAO-SES R 1992': ('big_poss2ukstu_red.fits', 'RP')}
out = {}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 4, figsize=(16, 8))
G = g.phot_g_mean_mag.values
for col, (lab, (fn, band)) in enumerate(plates.items()):
    h = fits.open(f'{B}/data/img/dss/{fn}')[0]; w = WCS(h.header); im = h.data.astype(float)
    ep = Time(h.header['DATE-OBS'][:10]).jyear
    scale = abs(float(h.header.get('PLTSCALE', 67.2))) * float(h.header.get('XPIXELSZ', 25.0)) / 1000.0
    gra, gde = pos(g.ra.values, g.dec.values, g.pmra.values, g.pmdec.values, ep)
    gx, gy = w.all_world2pix(gra, gde, 0)
    mref = (g.phot_bp_mean_mag if band == 'BP' else g.phot_rp_mean_mag).values
    ny, nx = im.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    inside = (gx > 12) & (gx < nx - 13) & (gy > 12) & (gy < ny - 13)
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([gx, gy]))
    import time; T0 = time.time()
    def iso(i, rad_as, dmag):
        nb = tree.query_ball_point([gx[i], gy[i]], rad_as / scale)
        return not any((j != i) and (G[j] < G[i] + dmag) for j in nb)
    # 1. alignment from centroids
    sel = [i for i in np.where(inside & (G > 13.5) & (G < 17))[0] if iso(i, 10, 4)]
    dxy = []
    for i in sel:
        xs, ys = int(round(gx[i])), int(round(gy[i])); s = 3
        c = im[ys - s:ys + s + 1, xs - s:xs + s + 1] - np.median(im[ys - 8:ys + 9, xs - 8:xs + 9])
        c = np.clip(c, 0, None)
        if c.sum() <= 0: continue
        cx = (c * xx[ys - s:ys + s + 1, xs - s:xs + s + 1]).sum() / c.sum(); cy = (c * yy[ys - s:ys + s + 1, xs - s:xs + s + 1]).sum() / c.sum()
        dxy.append((gx[i], gy[i], cx - gx[i], cy - gy[i]))
    dxy = np.array(dxy)
    A = np.column_stack([np.ones(len(dxy)), dxy[:, 0], dxy[:, 1]])
    keep = np.ones(len(dxy), bool)
    for it in range(5):
        cxf, *_ = np.linalg.lstsq(A[keep], dxy[keep, 2], rcond=None); cyf, *_ = np.linalg.lstsq(A[keep], dxy[keep, 3], rcond=None)
        rr = np.hypot(dxy[:, 2] - A @ cxf, dxy[:, 3] - A @ cyf); keep = rr < 3 * np.median(rr) + 0.1
    gx = gx + np.column_stack([np.ones(len(gx)), gx, gy]) @ cxf; gy = gy + np.column_stack([np.ones(len(gy)), gx, gy]) @ cyf
    align_rms = float(np.median(rr[keep]) * scale)
    # 2. PSF
    fw = []
    for i in [i for i in np.where(inside & (G > 16) & (G < 19.5))[0] if iso(i, 8, 3)][:80]:
        xs, ys = int(round(gx[i])), int(round(gy[i])); s = int(np.ceil(6 / scale))
        cut = im[ys - s:ys + s + 1, xs - s:xs + s + 1]; cx, cy = xx[ys - s:ys + s + 1, xs - s:xs + s + 1], yy[ys - s:ys + s + 1, xs - s:xs + s + 1]
        f = least_squares(lambda p: (cut - p[1] - p[0] * moffat(cx, cy, gx[i] + p[2], gy[i] + p[3], p[4], p[5])).ravel(),
                          [max(cut.max() - np.median(cut), 1) * 5, np.median(cut), 0, 0, 2.5 / scale, 2.5],
                          bounds=([0, -np.inf, -1, -1, 0.5, 1.05], [np.inf, np.inf, 1, 1, 8 / scale, 10]))
        fw.append(f.x[4:6])
    fw = np.array(fw); FW, BE = np.median(fw[:, 0]), np.median(fw[:, 1])
    # 3. scene photometry
    def scene(xc, yc, rad_as, test=None, inj=None, m3v_free=False, fixed=None):
        hb = int(np.ceil(rad_as / scale)); xs, ys = int(round(xc)), int(round(yc))
        sl = (slice(ys - hb, ys + hb + 1), slice(xs - hb, xs + hb + 1)); data = im[sl].copy(); cx, cy = xx[sl], yy[sl]
        if inj is not None: data = data + inj(cx, cy)
        nb = [j for j in tree.query_ball_point([xc, yc], hb * 1.42 + 4 / scale)]
        nb_fixed = [j for j in nb if not (m3v_free and j in (im3v, imira))]
        def design(p):
            cols = [np.ones(data.size), (cx - xc).ravel(), (cy - yc).ravel()]
            for j in nb_fixed: cols.append(moffat(cx, cy, gx[j] + p[0], gy[j] + p[1], FW, BE).ravel())
            if m3v_free:
                cols.append(moffat(cx, cy, gx[im3v] + p[0], gy[im3v] + p[1], p[2], p[3]).ravel())
                cols.append(moffat(cx, cy, gx[imira] + p[0] + p[6], gy[imira] + p[1] + p[7], p[4], p[5]).ravel())
            if test is not None: cols.append(moffat(cx, cy, test[0] + p[0], test[1] + p[1], FW, BE).ravel())
            return np.column_stack(cols)
        def res(p):
            Ad = design(p); c, *_ = np.linalg.lstsq(Ad, data.ravel(), rcond=None); return data.ravel() - Ad @ c
        if fixed is not None:
            px = np.array(fixed)
        else:
            p0 = [0, 0] + ([FW * 1.5, BE, FW * 2.5, BE, 0, 0] if m3v_free else [])
            lo = [-1, -1] + ([0.5, 1.05, 0.5, 1.05, -2, -2] if m3v_free else []); hi = [1, 1] + ([12 / scale, 10, 20 / scale, 10, 2, 2] if m3v_free else [])
            px = least_squares(res, p0, bounds=(lo, hi), max_nfev=200).x
        Ad = design(px); c, *_ = np.linalg.lstsq(Ad, data.ravel(), rcond=None)
        return dict(c=c, nb=nb_fixed, p=px, resid=(data.ravel() - Ad @ c).reshape(data.shape), ext=(xs - hb - .5, xs + hb + .5, ys - hb - .5, ys + hb + .5), data=data)
    calst = [i for i in np.where(inside & (G > 14.5) & (G < 21) & np.isfinite(mref))[0] if iso(i, 3.0, 0.0)]
    calst = list(rng.choice(calst, min(350, len(calst)), replace=False))
    amps = []
    for i in calst:
        try:
            s_ = scene(gx[i], gy[i], 10, fixed=[0.0, 0.0])
            amps.append((mref[i], s_['c'][3 + s_['nb'].index(i)], G[i], g.bp_rp.values[i]))
        except Exception:
            pass
    amps = np.array(amps); ok = amps[:, 1] > 0
    print(f'  {lab}: calibration fits done in {time.time()-T0:.0f}s', flush=True)
    np.save(f'{B}/out/dss_calib_{fn.replace(".fits","")}.npy', amps)
    # amplitude model: log10 amp = a0 + a1 (G-18) + a2 (G-18)^2 + a3 (BP-RP - 1.2), robust fit over G 15.5-21
    A_ = amps[ok]; la = np.log10(A_[:, 1]); Gc = A_[:, 2]; cc_ = np.nan_to_num(A_[:, 3], nan=1.2)
    sel_ = (Gc > 15.5) & (Gc < 21)
    X_ = np.column_stack([np.ones(len(la)), Gc - 18, (Gc - 18) ** 2, cc_ - 1.2])
    kk = sel_.copy()
    for it in range(8):
        cc, *_ = np.linalg.lstsq(X_[kk], la[kk], rcond=None); rr = la - X_ @ cc; sd = 1.4826 * np.median(np.abs(rr[kk])); kk = sel_ & (np.abs(rr) < 3 * sd)
    def G_of_T(T, c): return T + 0.00522555 * c ** 3 - 0.0891337 * c ** 2 + 0.633923 * c - 0.0324473
    def amp_pred(T, c):
        Gx = G_of_T(T, c); return 10 ** (cc[0] + cc[1] * (Gx - 18) + cc[2] * (Gx - 18) ** 2 + cc[3] * (c - 1.2))
    mono = bool(cc[1] + 2 * cc[2] * 3 < 0)   # decreasing up to G = 21
    # 4. target
    tx, ty = {}, {}
    for ept in (2019.4, 2026.3):
        a_, b_ = w.all_world2pix(*pos(r.ra, r.dec, r.pmra, r.pmdec, ept), 0)
        a_ = a_ + np.array([1, a_, b_]) @ cxf; b_ = b_ + np.array([1, a_, b_]) @ cyf; tx[ept], ty[ept] = float(a_), float(b_)
    xc, yc = 0.5 * (tx[2026.3] + gx[im3v]), 0.5 * (ty[2026.3] + gy[im3v])
    s0 = scene(xc, yc, 18, test=None, m3v_free=True); PF = list(s0['p'])
    res_t = {}
    for ept in (2019.4, 2026.3):
        s_ = scene(xc, yc, 18, test=(tx[ept], ty[ept]), m3v_free=True, fixed=PF); res_t[ept] = float(s_['c'][-1])
        if ept == 2026.3: S26 = s_
    # 5. empirical noise at random empty positions (in the target box: > 4" from Gaia, > 5" from the M3V) and in the field
    rnd = []
    s26 = np.hypot(tx[2026.3] - gx[im3v], ty[2026.3] - gy[im3v])
    for k in range(3000):
        rr_ = s26 * rng.uniform(0.75, 1.3); th = rng.uniform(0, 2 * np.pi)
        px, py = gx[im3v] + rr_ * np.cos(th), gy[im3v] + rr_ * np.sin(th)
        dmin = tree.query([px, py])[0] * scale
        if dmin > 3.5 and np.hypot(px - gx[imira], py - gy[imira]) * scale > 10:
            s_ = scene(xc, yc, 18, test=(px, py), m3v_free=True, fixed=PF); rnd.append(s_['c'][-1])
        if len(rnd) >= 60: break
    print(f'  random positions used: {len(rnd)}', flush=True)
    rnd = np.array(rnd); sig = 1.4826 * np.median(np.abs(rnd - np.median(rnd))); off = np.median(rnd)
    inj = {}
    for (cT, m) in [(c_, m_) for c_ in (1.0, 2.5) for m_ in np.arange(15.0, 21.01, 0.5)]:
        amp = amp_pred(m, cT); rec = []
        for k in range(8):
            jx, jy = tx[2026.3] + rng.normal(0, 0.3 / scale), ty[2026.3] + rng.normal(0, 0.3 / scale)
            s_ = scene(xc, yc, 18, test=(tx[2026.3], ty[2026.3]), m3v_free=True, fixed=PF, inj=lambda cx, cy: amp * moffat(cx, cy, jx, jy, FW, BE))
            rec.append(s_['c'][-1])
        rec = np.array(rec)
        inj[f'T{m:.1f}_c{cT}'] = dict(T=float(m), bp_rp=cT, amp_in=float(amp), rec_med=float(np.median(rec)), frac_3sig=float(np.mean((rec - off) > 3 * sig)))
    lim = {}
    for cT in (1.0, 2.5):
        det = [v['T'] for v in inj.values() if v['bp_rp'] == cT and v['frac_3sig'] >= 0.9]
        # faintest T such that all brighter T are also detected
        Ts = sorted(v['T'] for v in inj.values() if v['bp_rp'] == cT); L_ = None
        for Tt in Ts:
            if Tt in det: L_ = Tt
            else: break
        lim[f'bp_rp={cT}'] = L_
    snr = {e: float((v - off) / sig) for e, v in res_t.items()}
    # grid of test positions within 3" of the TESS-epoch positions: T limit from the local baseline + 3 sigma
    gridlim = []
    Tsc = np.linspace(12, 24, 1201)
    for ept in (2019.4, 2026.3):
        for rad in (0.0, 1.0, 2.0, 3.0):
            for th in (np.arange(8) * np.pi / 4 if rad > 0 else [0.0]):
                px, py = tx[ept] + rad / scale * np.cos(th), ty[ept] + rad / scale * np.sin(th)
                dmin = tree.query([px, py])[0] * scale
                if dmin < 1.5: continue
                s_ = scene(xc, yc, 18, test=(px, py), m3v_free=True, fixed=PF); base = s_['c'][-1]
                need = off + 3 * sig - base
                for cT in (1.0, 2.5):
                    ap = amp_pred(Tsc, cT); Tl = float(np.interp(-np.log10(max(need, 1e-9)), -np.log10(ap), Tsc))
                    gridlim.append(dict(epoch=ept, r=rad, th=float(th), bp_rp=cT, baseline_snr=float((base - off) / sig), T_lim=Tl,
                                        sep_from_star_then=float(np.hypot(px - gx[im3v], py - gy[im3v]) * scale)))
    def brentq_T(a):
        Ts = np.linspace(10, 24, 1401); ap = amp_pred(Ts, 1.0); return float(np.interp(-np.log10(a), -np.log10(ap), Ts))
    out[lab] = dict(file=fn, epoch=float(ep), scale=float(scale), band=band, align_rms_arcsec=align_rms, psf_fwhm_arcsec=float(FW * scale), psf_beta=float(BE),
                    calib=dict(model='log10 amp = a0 + a1 (G-18) + a2 (G-18)^2 + a3 (BP-RP-1.2)', coef=cc.tolist(), robust_sd_dex=float(sd), n=int(kk.sum()), monotonic=bool(mono)),
                    sep_2026_from_star_then=float(np.hypot(tx[2026.3] - gx[im3v], ty[2026.3] - gy[im3v]) * scale),
                    sep_2019_from_star_then=float(np.hypot(tx[2019.4] - gx[im3v], ty[2019.4] - gy[im3v]) * scale),
                    test_amp=res_t, test_snr=snr, noise_amp=float(sig), noise_offset=float(off), n_random=int(len(rnd)),
                    injection=inj, limit_T_3sig_90pct=lim, grid_limits=gridlim,
                    grid_T_lim_min={f'bp_rp={c_}': float(min(x['T_lim'] for x in gridlim if x['bp_rp'] == c_)) for c_ in (1.0, 2.5)} if gridlim else None,
                    grid_T_lim_median={f'bp_rp={c_}': float(np.median([x['T_lim'] for x in gridlim if x['bp_rp'] == c_])) for c_ in (1.0, 2.5)} if gridlim else None,
                    grid_T_lim_min_within2={f'bp_rp={c_}': float(min(x['T_lim'] for x in gridlim if x['bp_rp'] == c_ and x['r'] <= 2)) for c_ in (1.0, 2.5)} if gridlim else None, test_T_equiv_c1=(None if res_t[2026.3] <= 0 else float(brentq_T(res_t[2026.3]))), m3v_fit_fwhm_arcsec=float(s0['p'][2] * scale))
    print(f"{lab}: {scale:.2f}\"/px align rms {align_rms:.2f}\" PSF FWHM {FW*scale:.1f}\" beta {BE:.1f}; calib n={kk.sum()} sd {sd:.2f} mag (monotonic {mono}); "
          f"sep 2026/2019 from star then {out[lab]['sep_2026_from_star_then']:.1f}/{out[lab]['sep_2019_from_star_then']:.1f}\"; "
          f"test S/N 2019 {snr[2019.4]:+.1f}, 2026 {snr[2026.3]:+.1f}; T limit {lim}; calib sd {sd:.2f} dex; grid (<=3\" of 2019/2026 pos) T_lim min/median "
          f"{out[lab]['grid_T_lim_min']} / {out[lab]['grid_T_lim_median']}", flush=True)
    for row, img in ((0, S26['data']), (1, S26['resid'])):
        ax = axs[row, col]; v1, v2 = np.percentile(img, [2, 99.7]) if row == 0 else (-4 * np.std(img), 4 * np.std(img))
        ax.imshow(img, origin='lower', cmap='gray_r' if row == 0 else 'RdBu_r', vmin=v1, vmax=v2, extent=S26['ext'])
        ax.plot(tx[2026.3], ty[2026.3], '+', color='C3', ms=14, mew=1.5); ax.plot(tx[2019.4], ty[2019.4], '+', color='C1', ms=14, mew=1.5)
        ax.plot(gx[im3v], gy[im3v], 'x', color='C2', ms=10)
        ax.set_title(f'{lab}: ' + ('data' if row == 0 else f'residual; 90% limit T > {lim}'), fontsize=8)
fig.tight_layout(); fig.savefig(f'{B}/figs/dss_resid.png', dpi=100)
json.dump(out, open(f'{B}/out/dss_limits.json', 'w'), indent=1, default=float)
