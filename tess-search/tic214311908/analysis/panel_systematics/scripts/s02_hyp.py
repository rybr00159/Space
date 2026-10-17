"""Which source does the dip image look like?  Raw-pixel depth images (s01) against templates, own code.

Per sector: the mean real depth image within R px of the M3V is fitted with
   H_M3V : a * PRF(M3V apparent position)          + nuisance
   H_Mira: a * Mira empirical image (median image minus PRFs of all other Gaia stars) + nuisance
   H_pairA, H_G175 (3.1" star), H_G177 (5.5" star) : a * PRF(star)          + nuisance
   H_pt  : nuisance only (pointing dipole + constant)
nuisance = gx, gy (gradients of the median image, i.e. a pointing excursion) + constant.
Per-pixel noise: robust std of the fake-window images / sqrt(N_real).  Delta-chi2 calibrated with injections of each
template (14 e-/s total flux for point sources; the Mira template scaled to the same in-region flux as the M3V PRF)
into random sets of fake windows.  Also: free-position PRF grid, box-aperture depths, SPOC-background and pointing
excursions at transit vs fake epochs, and the Mira flux-scaling test on the box depths.
Output: out/s02_hyp.json
"""
import sys, numpy as np
sys.path.insert(0, 'tic214311908/analysis/panel_systematics/scripts')
from ps_lib import *

TAG = sys.argv[1] if len(sys.argv) > 1 else 'pt'
R = 3.5
NINJ = 400
rng = np.random.default_rng(1)
g = gaia()
names = {v: k for k, v in ID.items()}


def rob(x, axis=0):
    return 1.4826 * np.median(np.abs(x - np.median(x, axis)), axis)


res = dict(tag=TAG, sectors={})
chi_tot = {}
inj_tot = {}
grid_tot = None
for sec in SECTORS:
    c = Cube(sec)
    d = np.load(f'{ROOT}/out/s01_img_s{sec}.npz')
    pos = d['pos']; sid = d['source_id']
    ip = {n: int(np.where(sid == ID[n])[0][0]) for n in ID}
    xm3, ym3 = pos[ip['M3V']]
    Ar, Af = d[f'{TAG}_real_A'], d[f'{TAG}_fake_A']
    N = len(Ar)
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    reg = np.hypot(xx - xm3, yy - ym3) <= R
    med = d['mean']
    gy_, gx_ = np.gradient(med)
    # Mira template: median image minus PRF models of the other stars within 6 px
    model_others = np.zeros_like(med)
    for i, (x, y) in enumerate(pos):
        if sid[i] == ID['Mira'] or not (g['T'].values[i] <= 17.8):
            continue
        if np.hypot(x - xm3, y - ym3) > 7:
            continue
        model_others += g['F'].values[i] * c.prf(x, y)
    mira_t = med - model_others
    tm = {n: c.prf(*pos[ip[n]]) for n in ('M3V', 'pairA', 'pairB', 'G175', 'G177', 'G170')}
    # scale Mira template to the same in-region sum as the M3V PRF
    tm['Mira'] = mira_t / mira_t[reg].sum() * tm['M3V'][reg].sum()
    sig_pix = rob(Af, 0) / np.sqrt(N)
    w = 1 / sig_pix[reg] ** 2
    nuis = np.column_stack([gx_[reg], gy_[reg], np.ones(reg.sum())])
    nuis = nuis / np.sqrt((nuis ** 2).mean(0))

    def fit(img, T=None):
        X = nuis if T is None else np.column_stack([T[reg], nuis])
        Xw = X * np.sqrt(w)[:, None]; yw = img[reg] * np.sqrt(w)
        b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ b
        return float(r @ r), b

    def allfits(img):
        out = {}
        for n in ('M3V', 'Mira', 'pairA', 'pairB', 'G175', 'G177', 'G170'):
            out[n] = fit(img, tm[n])
        out['pt'] = fit(img)
        return out

    mimg = Ar.mean(0)
    F = allfits(mimg)
    chi = {k: v[0] for k, v in F.items()}
    amp = {k: float(v[1][0]) for k, v in F.items() if k != 'pt'}
    for k, v in chi.items():
        chi_tot[k] = chi_tot.get(k, 0) + v
    # injections: M3V and Mira as truth
    inj = {}
    for truth in ('M3V', 'Mira', 'pairA', 'G175', 'none'):
        L = []
        for _ in range(NINJ):
            ii = rng.choice(len(Af), N, replace=False)
            img = Af[ii].mean(0) + (0 if truth == 'none' else 14.0 * tm[truth])
            f = allfits(img)
            L.append([f[k][0] for k in ('M3V', 'Mira', 'pairA', 'pairB', 'G175', 'G177', 'G170', 'pt')])
        inj[truth] = np.array(L)
        inj_tot.setdefault(truth, []).append(inj[truth])
    # free-position grid in sky offsets from the M3V (arcsec)
    ge = np.arange(-24, 24.01, 1.0)
    G = np.zeros((len(ge), len(ge)))
    for i, dn in enumerate(ge):
        for j, de in enumerate(ge):
            dpx = c.Ji @ np.array([de, dn])
            G[i, j] = fit(mimg, c.prf(xm3 + dpx[0], ym3 + dpx[1]))[0]
    grid_tot = G if grid_tot is None else grid_tot + G
    ib = np.unravel_index(np.argmin(G), G.shape)
    # Mira empirical-image light centroid and M3V offsets (arcsec)
    mreg = np.hypot(xx - xm3, yy - ym3) <= 3
    mt = np.clip(mira_t, 0, None) * mreg
    cxm, cym = (mt * xx).sum() / mt.sum(), (mt * yy).sum() / mt.sum()
    mira_cen_off = c.J @ np.array([cxm - xm3, cym - ym3])
    mira_gaia_off = c.J @ (pos[ip['Mira']] - pos[ip['M3V']])
    # box apertures
    ix, iy = int(round(xm3)), int(round(ym3))
    box = np.zeros_like(med, bool); box[iy - 1:iy + 2, ix - 1:ix + 2] = True
    r25 = np.hypot(xx - xm3, yy - ym3) <= 2.5
    ee_box = tm['M3V'][box].sum(); ee_r25 = tm['M3V'][r25].sum()
    br = Ar[:, box].sum(1); bf = Af[:, box].sum(1)
    rr = Ar[:, r25].sum(1); rf = Af[:, r25].sum(1)
    Abr = d[f'{TAG}_real_Ab'][:, box].sum(1); Abf = d[f'{TAG}_fake_Ab'][:, box].sum(1)
    ptr = np.hypot(d[f'{TAG}_real_ptx'], d[f'{TAG}_real_pty']); ptf = np.hypot(d[f'{TAG}_fake_ptx'], d[f'{TAG}_fake_pty'])
    mira_flux_win = d[f'{TAG}_real_mean_img'][:, box].sum(1)
    res['sectors'][sec] = dict(
        N=N, nfake=len(Af), chi=chi, amp=amp, npix=int(reg.sum()),
        dchi_vs_M3V={k: chi[k] - chi['M3V'] for k in chi},
        inj_dchi_q={t: {k: np.percentile(inj[t][:, j] - inj[t][:, 0], [2.5, 50, 97.5]).tolist()
                        for j, k in enumerate(('M3V', 'Mira', 'pairA', 'pairB', 'G175', 'G177', 'G170', 'pt'))}
                    for t in inj},
        free_best_off_arcsec=[float(ge[ib[1]]), float(ge[ib[0]])], free_dchi_vs_M3V=float(G.min() - chi['M3V']),
        mira_template_centroid_off_arcsec=mira_cen_off.tolist(), mira_gaia_off_arcsec=mira_gaia_off.tolist(),
        box_depth_e=float(br.mean()), box_err_e=float(rob(bf) / np.sqrt(N)), box_fake_mean=float(bf.mean()),
        box_per_window=br.tolist(), E=d[f'{TAG}_real_E'].tolist(), tc=d[f'{TAG}_real_tc'].tolist(),
        r25_depth_e=float(rr.mean()), r25_err_e=float(rob(rf) / np.sqrt(N)),
        ee_box=float(ee_box), ee_r25=float(ee_r25),
        box_depth_ppt=float(br.mean() / ee_box / FM3V * 1e3), box_err_ppt=float(rob(bf) / np.sqrt(N) / ee_box / FM3V * 1e3),
        spocbkg_box_real=float(Abr.mean()), spocbkg_box_err=float(rob(Abf) / np.sqrt(N)),
        spocbkg_box_per_window=Abr.tolist(), spocbkg_fake_rob=float(rob(Abf)),
        pointing_real_px=ptr.tolist(), pointing_fake_p50_p95_max=[float(np.median(ptf)), float(np.percentile(ptf, 95)), float(ptf.max())],
        pointing_real_signed=np.column_stack([d[f'{TAG}_real_ptx'], d[f'{TAG}_real_pty']]).tolist(),
        pointing_mean_real=[float(d[f'{TAG}_real_ptx'].mean()), float(d[f'{TAG}_real_pty'].mean())],
        pointing_err=[float(rob(d[f'{TAG}_fake_ptx']) / np.sqrt(N)), float(rob(d[f'{TAG}_fake_pty']) / np.sqrt(N))],
        box_flux_win=mira_flux_win.tolist(),
        qbits_in=d[f'{TAG}_real_qbits_in'].tolist(), nflag_in=d[f'{TAG}_real_nflag_in'].tolist(),
        nin_all=d[f'{TAG}_real_nin_all'].tolist(), dump_dt_d=d[f'{TAG}_real_dump_dt'].tolist(),
        fake_dump_dt_frac_lt_0p1=float(np.mean(d[f'{TAG}_fake_dump_dt'] < 0.1)),
        bkg_in=d[f'{TAG}_real_bkg_in'].tolist(),
    )
    print(sec, N, {k: round(v, 1) for k, v in res['sectors'][sec]['dchi_vs_M3V'].items()},
          'box', round(br.mean(), 1), '+-', round(rob(bf) / np.sqrt(N), 1), 'spocbkg', round(Abr.mean(), 2),
          'free', res['sectors'][sec]['free_best_off_arcsec'], flush=True)

res['combined_dchi_vs_M3V'] = {k: chi_tot[k] - chi_tot['M3V'] for k in chi_tot}
keys = ('M3V', 'Mira', 'pairA', 'pairB', 'G175', 'G177', 'G170', 'pt')
comb = {}
for t, L in inj_tot.items():
    S = np.sum(np.array(L), 0)            # sum over sectors (independent draws)
    comb[t] = {k: np.percentile(S[:, j] - S[:, 0], [0.15, 2.5, 50, 97.5, 99.85]).tolist() for j, k in enumerate(keys)}
    comb[t]['frac_reaching_real'] = {k: float(np.mean((S[:, j] - S[:, 0]) >= res['combined_dchi_vs_M3V'][k]))
                                     for j, k in enumerate(keys)}
    comb[t]['frac_le_real'] = {k: float(np.mean((S[:, j] - S[:, 0]) <= res['combined_dchi_vs_M3V'][k]))
                               for j, k in enumerate(keys)}
res['combined_inj'] = comb
ib = np.unravel_index(np.argmin(grid_tot), grid_tot.shape)
ge = np.arange(-24, 24.01, 1.0)
res['combined_free_best_off_arcsec'] = [float(ge[ib[1]]), float(ge[ib[0]])]
res['combined_free_dchi_vs_M3V'] = float(grid_tot.min() - chi_tot['M3V'])
np.save(f'{ROOT}/out/s02_grid_{TAG}.npy', grid_tot)
# box depth summary
bd = np.array([res['sectors'][s]['box_depth_ppt'] for s in SECTORS]); be = np.array([res['sectors'][s]['box_err_ppt'] for s in SECTORS])
wm = np.sum(bd / be ** 2) / np.sum(1 / be ** 2)
res['box_ppt_wmean'] = [float(wm), float(1 / np.sqrt(np.sum(1 / be ** 2)))]
res['box_ppt_chi2_const'] = float(np.sum((bd - wm) ** 2 / be ** 2))
bde = np.array([res['sectors'][s]['box_depth_e'] for s in SECTORS]); bee = np.array([res['sectors'][s]['box_err_e'] for s in SECTORS])
wme = np.sum(bde / bee ** 2) / np.sum(1 / bee ** 2)
res['box_e_wmean'] = [float(wme), float(1 / np.sqrt(np.sum(1 / bee ** 2)))]
res['box_e_chi2_const'] = float(np.sum((bde - wme) ** 2 / bee ** 2))
sb = np.array([res['sectors'][s]['spocbkg_box_real'] for s in SECTORS]); se = np.array([res['sectors'][s]['spocbkg_box_err'] for s in SECTORS])
res['spocbkg_box_wmean'] = [float(np.sum(sb / se ** 2) / np.sum(1 / se ** 2)), float(1 / np.sqrt(np.sum(1 / se ** 2)))]
print('combined', {k: round(v, 1) for k, v in res['combined_dchi_vs_M3V'].items()})
print('free best', res['combined_free_best_off_arcsec'], res['combined_free_dchi_vs_M3V'])
print('box ppt', res['box_ppt_wmean'], res['box_ppt_chi2_const'], 'box e', res['box_e_wmean'], res['box_e_chi2_const'])
print('spoc bkg', res['spocbkg_box_wmean'])
for t in comb:
    print('inj', t, {k: [round(x, 1) for x in comb[t][k]] for k in ('Mira', 'pairA', 'G175', 'pt')})
jdump(res, f'{ROOT}/out/s02_hyp_{TAG}.json')
