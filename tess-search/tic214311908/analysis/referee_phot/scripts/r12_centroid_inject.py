"""Referee r12: pixel-level injection test of the difference-image centroid (calibrates r10).
For each trial and sector, n_real fake windows (n_real = number of real strict transits used in r10) receive a
15 e-/s transit (r06 global shape, at the window centre) injected with the PRF of one source: M3V, bg3 (Gaia DR3
5925373272827260544, ~3" NE), the pair, or the M3V with a mismatched PRF (displaced 0.1 px, blurred 0.15 px).
The stacked flat-part residual frame is fitted at free position exactly as in r10 (same pixel errors from fake
stacks); the six sector offsets are combined with the same weights as r10. 40 trials per hypothesis.
Output: out/centroid_inject.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, fake_centres
from scipy.optimize import least_squares
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T14, FL = G["T14_h"] / 24, G["flat"]
D = json.load(open(f"{R}/out/diffimg.json"))
rng = np.random.default_rng(17)
NT = 40
HYP = ("M3V", "bg3", "pair", "M3V_mis")
res = {h: [] for h in HYP}
per = {}
for s in SECTORS:
    S = Sector(s); C = S.C
    z = np.load(f"{R}/data/win_s{s}.npz", allow_pickle=True)
    nreal = D["per_sector"][str(s)]["n_transits"]
    fd = -z["fake_dimg"]; fn = z["fake_n_flat"]
    st = []
    for k in range(300):
        j = rng.choice(len(fd), nreal, replace=False); ww = fn[j] / fn[j].sum(); st.append(np.nansum(fd[j] * ww[:, None], 0))
    err = np.std(st, 0)
    m3 = S.pos["M3V"]; mi = S.pos["Mira"]; yy, xx = S.yy.ravel(), S.xx.ravel()
    fm = ((np.hypot(xx - m3[0], yy - m3[1]) <= 3) | (np.hypot(xx - mi[0], yy - mi[1]) <= 3)) & (err > 0)
    wt = 1 / err[fm]
    n_ = np.array([-2, -1, 0, 1, 2]) * 0.12; wg = np.exp(-0.5 * (n_ / 0.15) ** 2); wg /= wg.sum()
    imgs = dict(M3V=S.P_M3V, bg3=C.psf(*S.pos["bg3"]), pair=S.P_pair,
                M3V_mis=sum(wg[i] * wg[j] * C.psf(m3[0] + 0.1 + n_[i], m3[1] - 0.1 + n_[j]) for i in range(5) for j in range(5)))
    cents = fake_centres(S, 120, np.random.default_rng(1000 + s))
    per[s] = {}
    for h in HYP:
        offs = []
        for trial in range(NT):
            sel = rng.choice(len(cents), nreal, replace=False)
            dims = []; nfl = []
            for tc in cents[sel]:
                r = S.window(tc, inject=[(imgs[h], 15.0, tc, T14, FL)])
                if r is None: continue
                dims.append(-r["dimg"]); nfl.append(r["n_flat"])
            if len(dims) < max(1, nreal - 1): continue
            ww = np.array(nfl) / np.sum(nfl); dI = np.nansum(np.array(dims) * ww[:, None], 0)
            ok = fm & np.isfinite(dI)
            def rr(p):
                A = np.column_stack([C.psf(p[0], p[1]).ravel()[ok], np.ones(ok.sum())]) / err[ok][:, None]
                c, *_ = np.linalg.lstsq(A, dI[ok] / err[ok], rcond=None); return A @ c - dI[ok] / err[ok]
            best = min(((np.sum(rr([m3[0] + a, m3[1] + b]) ** 2), [m3[0] + a, m3[1] + b]) for a in np.linspace(-1.2, 1.2, 9) for b in np.linspace(-1.2, 1.2, 9)), key=lambda q: q[0])
            q = least_squares(rr, best[1], diff_step=1e-3, x_scale=0.05)
            offs.append(C.en(m3[0], m3[1], q.x[0] - m3[0], q.x[1] - m3[1]))
        offs = np.array(offs); per[s][h] = offs
        print(f"S{s} {h:8s}: n {len(offs)} mean offset (E,N) {np.round(offs.mean(0),2).tolist()} std {np.round(offs.std(0),2).tolist()} (robust {rstd(offs[:,0]):.2f}, {rstd(offs[:,1]):.2f})", flush=True)
# combine sectors trial-by-trial with the r10 weights (from the r10 per-sector errors)
wE = np.array([1 / D["per_sector"][str(s)]["fits"]["free"]["off_err"][0] ** 2 for s in SECTORS]); wN = np.array([1 / D["per_sector"][str(s)]["fits"]["free"]["off_err"][1] ** 2 for s in SECTORS])
summ = {}
for h in HYP:
    n = min(len(per[s][h]) for s in SECTORS)
    E = np.array([per[s][h][:n, 0] for s in SECTORS]); N = np.array([per[s][h][:n, 1] for s in SECTORS])
    cE = (wE[:, None] * E).sum(0) / wE.sum(); cN = (wN[:, None] * N).sum(0) / wN.sum()
    summ[h] = dict(n=n, mean=[float(cE.mean()), float(cN.mean())], std=[float(cE.std()), float(cN.std())], robust=[float(rstd(cE)), float(rstd(cN))],
                   combined=np.column_stack([cE, cN]).tolist())
    print(f"combined {h:8s}: mean (E,N) {np.round([cE.mean(), cN.mean()],2).tolist()} std {np.round([cE.std(), cN.std()],2).tolist()}")
obs = np.array([D["centroid_E"], D["centroid_N"]])
for h in HYP:
    c = np.array(summ[h]["combined"]); mu = c.mean(0); cv = np.cov(c.T)
    d2 = float((obs - mu) @ np.linalg.inv(cv) @ (obs - mu)); summ[h]["mahal2_obs"] = d2
    print(f"observed centroid {np.round(obs,2).tolist()} vs injected-{h} distribution: Mahalanobis^2 {d2:.1f} (2 dof)")
save_json(dict(summary=summ, observed=obs.tolist(), bg3_EN=D["bg3_EN"]), f"{R}/out/centroid_inject.json")
