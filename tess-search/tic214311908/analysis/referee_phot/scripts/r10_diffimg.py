"""Referee r10: in-transit difference images from the r05 local-window residual frames (flat part, |t-tc|<0.02 d;
per-pixel baselines fitted out of transit within each window, so no global Mira template is involved).
Per sector: dI = - weighted mean over real strict windows of the mean residual frame (positive = dimming), noise
per pixel from the same statistic over random subsets of fake windows (same number of windows), pixels within 3 px of
the M3V or Mira. Fits: constant + PRF(M3V), PRF(pair A+B), PRF(Mira), b-image(Mira), PRF(bg3), PRF at free position.
Combined: chi2 summed over sectors, and the free-position centroid (E, N arcsec from the M3V) per sector with errors
from refits of dI + (fake-mean-subtracted) fake stacks; weighted mean. Output: out/diffimg.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector
from scipy.optimize import least_squares
rng = np.random.default_rng(5)
out = {}; tot = {}; cen = []
for s in SECTORS:
    S = Sector(s); z = np.load(f"{R}/data/win_s{s}.npz", allow_pickle=True)
    ep = z["strict_epoch"]; keep = ~((s == 12) & (ep == -388))
    di = -z["strict_dimg"][keep]; nflat = z["strict_n_flat"][keep]
    fd = -z["fake_dimg"]; fn = z["fake_n_flat"]
    n = len(di)
    w = nflat / nflat.sum(); dI = np.nansum(di * w[:, None], 0)
    stacks = []
    for k in range(300):
        j = rng.choice(len(fd), n, replace=False); ww = fn[j] / fn[j].sum(); stacks.append(np.nansum(fd[j] * ww[:, None], 0))
    stacks = np.array(stacks); err = np.std(stacks, 0); mu0 = stacks.mean(0)
    C = S.C; yy, xx = S.yy.ravel(), S.xx.ravel()
    m3 = S.pos["M3V"]; mi = S.pos["Mira"]
    fm = ((np.hypot(xx - m3[0], yy - m3[1]) <= 3) | (np.hypot(xx - mi[0], yy - mi[1]) <= 3)) & np.isfinite(dI) & (err > 0)
    wt = 1 / err[fm]
    T = dict(M3V=S.P_M3V.ravel(), pair=S.P_pair.ravel(), Mira_prf=S.P_Mira.ravel(), Mira_b=S.b_Mira.ravel(), bg3=C.psf(*S.pos["bg3"]).ravel())
    fits_ = {}
    for nm, tp in T.items():
        A = np.column_stack([tp[fm], np.ones(fm.sum())]) * wt[:, None]
        c, *_ = np.linalg.lstsq(A, dI[fm] * wt, rcond=None)
        fits_[nm] = dict(amp=float(c[0]), chi2=float(np.sum((A @ c - dI[fm] * wt) ** 2)))
    def rr(p, img):
        A = np.column_stack([C.psf(p[0], p[1]).ravel()[fm], np.ones(fm.sum())]) * wt[:, None]
        c, *_ = np.linalg.lstsq(A, img[fm] * wt, rcond=None); return A @ c - img[fm] * wt
    best = min(((np.sum(rr([m3[0] + a, m3[1] + b], dI) ** 2), [m3[0] + a, m3[1] + b]) for a in np.linspace(-1.5, 1.5, 13) for b in np.linspace(-1.5, 1.5, 13)), key=lambda q: q[0])
    q = least_squares(lambda p: rr(p, dI), best[1], diff_step=1e-3, x_scale=0.05)
    off = C.en(m3[0], m3[1], q.x[0] - m3[0], q.x[1] - m3[1])
    pos = []
    for k in range(80):
        img = dI + (stacks[k] - mu0)
        q2 = least_squares(lambda p: rr(p, img), q.x, diff_step=1e-3, x_scale=0.05)
        pos.append(C.en(m3[0], m3[1], q2.x[0] - m3[0], q2.x[1] - m3[1]))
    pos = np.array(pos)
    fits_["free"] = dict(chi2=float(np.sum(q.fun ** 2)), off_EN=off.tolist(), off_err=pos.std(0).tolist(), off_err_robust=[rstd(pos[:, 0]), rstd(pos[:, 1])])
    rel = {k: C.en(m3[0], m3[1], S.pos[k][0] - m3[0], S.pos[k][1] - m3[1]).tolist() for k in ("Mira", "pA", "bg3")}
    rel["pair"] = C.en(m3[0], m3[1], (S.pos["pA"][0] + 0.367 * S.pos["pB"][0]) / 1.367 - m3[0], (S.pos["pA"][1] + 0.367 * S.pos["pB"][1]) / 1.367 - m3[1]).tolist()
    out[s] = dict(n_transits=int(n), npix=int(fm.sum()), fits=fits_, rel_from_M3V=rel)
    cen.append((s, off[0], off[1], pos.std(0)[0], pos.std(0)[1]))
    for k_, v in fits_.items(): tot[k_] = tot.get(k_, 0) + v["chi2"]
    print(f"S{s}: n {n}, {fm.sum()} px; chi2 " + ", ".join(f"{k_} {v['chi2']:.1f}" + (f" (amp {v['amp']:.1f})" if 'amp' in v else "") for k_, v in fits_.items()) +
          f"; free offset (E,N) {np.round(off,1).tolist()} +- {np.round(pos.std(0),1).tolist()}\"; bg3 at {np.round(rel['bg3'],1).tolist()}, pair {np.round(rel['pair'],1).tolist()}", flush=True)
dof = sum(v["npix"] - 2 for v in out.values())
cen = np.array(cen); wE = 1 / cen[:, 3] ** 2; wN = 1 / cen[:, 4] ** 2
mE = float(np.sum(wE * cen[:, 1]) / wE.sum()); mN = float(np.sum(wN * cen[:, 2]) / wN.sum())
sE = float(1 / np.sqrt(wE.sum())); sN = float(1 / np.sqrt(wN.sum()))
chiE = float(np.sum(wE * (cen[:, 1] - mE) ** 2)); chiN = float(np.sum(wN * (cen[:, 2] - mN) ** 2))
bgE, bgN = np.mean([out[s]["rel_from_M3V"]["bg3"] for s in SECTORS], 0)
pE, pN = np.mean([out[s]["rel_from_M3V"]["pair"] for s in SECTORS], 0)
res = dict(per_sector=out, total_chi2=tot, total_dof=dof, centroid_E=mE, centroid_E_err=sE, centroid_N=mN, centroid_N_err=sN, chi2_sector_scatter=(chiE, chiN),
           bg3_EN=[float(bgE), float(bgN)], chi2_bg3_stat=float((mE - bgE) ** 2 / sE ** 2 + (mN - bgN) ** 2 / sN ** 2),
           chi2_bg3_sys15=float((mE - bgE) ** 2 / (sE ** 2 + 1.5 ** 2) + (mN - bgN) ** 2 / (sN ** 2 + 1.5 ** 2)),
           chi2_M3V_stat=float(mE ** 2 / sE ** 2 + mN ** 2 / sN ** 2), pair_EN=[float(pE), float(pN)],
           chi2_pair_stat=float((mE - pE) ** 2 / sE ** 2 + (mN - pN) ** 2 / sN ** 2))
save_json({str(k): v for k, v in res.items()}, f"{R}/out/diffimg.json")
print("total chi2:", {k: round(v, 1) for k, v in tot.items()}, "dof", dof)
print(f"centroid from M3V: E {mE:+.2f}+-{sE:.2f}, N {mN:+.2f}+-{sN:.2f} (sector scatter chi2 {chiE:.1f}, {chiN:.1f} / {len(cen)-1}); bg3 at ({bgE:.1f},{bgN:.1f}): chi2 {res['chi2_bg3_stat']:.1f} (stat) / {res['chi2_bg3_sys15']:.1f} (+1.5\" sys); M3V chi2 {res['chi2_M3V_stat']:.1f}; pair chi2 {res['chi2_pair_stat']:.0f}")
