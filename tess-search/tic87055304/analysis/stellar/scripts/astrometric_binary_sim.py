"""Can an unresolved binary orbit that produces the observed Gaia DR3 astrometric excess noise (2.32 mas) bias the
parallax enough for the true distance to be the GSP-Phot 998 pc (true parallax ~1.0 mas) while DR3 reports
2.807 +- 0.357 mas (and DR2 3.512 +- 0.735 mas)?

Method (forward simulation, per true parallax plx_t on a grid 0.6-4.5 mas):
  * Real DR3 scan times and scan position angles at the target position from the Gaia DR3 nominal scanning law
    (gaiaunlimited, 'dr3_nominal' with DR3 astrometry gaps); DR2 = the subset of scans before TCB day 2334
    (2016 May 23). Each transit = one along-scan (AL) measurement with sigma = 0.12/sqrt(8) mas.
  * Parallax factors from the barycentric Earth position (Gaia at L2 ignored, ~1 %).
  * Binary photocentre orbits: P log-normal (Raghavan+10: mu log10 P[d] = 5.03, sigma 2.28) truncated to 30 d-100 yr,
    e ~ U(0, 0.8) (e < 0.1 for P < 12 d), isotropic orientation, random phase, primary mass M1 (1.0 Msun for
    plx_t > 2, 1.2 Msun otherwise), q ~ U(0.05, 1.0) with a luminous MS companion (G-band flux ratio f = q^4) or, in a
    second prior, a dark companion (f = 0, q ~ U(0.05, 1.5)). Photocentre semi-major axis a_ph = a_rel (q/(1+q) - f/(1+f)).
  * 5-parameter weighted least squares fit per release; excess noise eps from the residual rms
    (eps^2 = sum(res^2)/(N-5) - sigma^2); fitted parallax plx_fit.
  * A third, extreme prior allows a dark companion up to q = 3 (e.g. a black hole) to test robustness.
  * Likelihood L(plx_t) = fraction of simulated binaries with DR3 eps in [1.9, 2.8] mas multiplied by the density of
    plx_fit(DR3) at 2.832 mas (L21 zero point applied) -> joint likelihood of (eps, plx). A version also requiring
    DR2 plx_fit within 3.54 +- 0.735 mas is given.
Output: out/astrometric_sim.json, fig/astrometric_sim.png
"""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["GAIAUNLIMITED_DATADIR"] = "/home/user/work/tic87/stellar/cache/gaiaunlimited"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
import astropy.units as u
ROOT = "/home/user/work/tic87/stellar"
RA, DEC = np.radians(272.60647748), np.radians(-43.99546646)
cache = os.path.join(ROOT, "out/dr3_scans.npz")
if not os.path.exists(cache):
    from gaiaunlimited import scanninglaw
    sl = scanninglaw.GaiaScanningLaw(version="dr3_nominal", gaplist="dr3/Astrometry")
    t1, t2 = sl.query(272.60647748, -43.99546646, count_only=False)
    pdf = sl.pointingdf
    tt = pdf.tcb_at_gaia.values
    def ang(times, col):
        i = np.searchsorted(tt, times); i = np.clip(i, 0, len(tt) - 1)
        return pdf[col].values[i]
    times = np.concatenate([t1, t2]); psi = np.concatenate([ang(t1, "scan_angle_fov1"), ang(t2, "scan_angle_fov2")])
    o = np.argsort(times)
    np.savez(cache, t=times[o], psi=psi[o])
z = np.load(cache); t_tcb, psi_deg = z["t"], z["psi"]
jd = t_tcb + 2455197.5
print("DR3 scans:", len(jd), Time(jd[0], format="jd").iso, Time(jd[-1], format="jd").iso)
pos = get_body_barycentric("earth", Time(jd, format="jd", scale="tcb"))
X, Y, Z = pos.x.to(u.au).value, pos.y.to(u.au).value, pos.z.to(u.au).value
Pa = X * np.sin(RA) - Y * np.cos(RA)
Pd = X * np.cos(RA) * np.sin(DEC) + Y * np.sin(RA) * np.sin(DEC) - Z * np.cos(DEC)
psi = np.radians(psi_deg)
sa, ca = np.sin(psi), np.cos(psi)
ty = (jd - 2457389.0) / 365.25  # years from J2016.0
SIG = 0.12 / np.sqrt(8)


def design(m):
    return np.column_stack([sa[m], ca[m], ty[m] * sa[m], ty[m] * ca[m], Pa[m] * sa[m] + Pd[m] * ca[m]])

m3 = np.ones(len(jd), bool); m2 = t_tcb < 2334.0
A3, A2 = design(m3), design(m2)
H3 = np.linalg.solve(A3.T @ A3, A3.T); H2 = np.linalg.solve(A2.T @ A2, A2.T)
R3 = np.eye(m3.sum()) - A3 @ H3; R2 = np.eye(m2.sum()) - A2 @ H2
print("DR2 scans:", m2.sum())
rng = np.random.default_rng(7)


def kepler_E(M, e):
    E = M.copy()
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def simulate(plx_t, n, dark=False):
    M1 = 1.0 if plx_t > 2 else 1.2
    lp = rng.normal(5.03, 2.28, 4 * n); lp = lp[(lp > np.log10(30)) & (lp < np.log10(100 * 365.25))][:n]
    n = len(lp); P = 10 ** lp / 365.25
    e = rng.uniform(0, 0.8, n); e[P < 12 / 365.25] = rng.uniform(0, 0.1, (P < 12 / 365.25).sum())
    if dark == "extreme":
        q = rng.uniform(0.05, 3.0, n); f = np.zeros(n)
    elif dark:
        q = rng.uniform(0.05, 1.5, n); f = np.zeros(n)
    else:
        q = rng.uniform(0.05, 1.0, n); f = q ** 4
    B = q / (1 + q) - f / (1 + f)
    a_rel = (M1 * (1 + q) * P ** 2) ** (1 / 3)
    a_ph = a_rel * B * plx_t  # mas
    cosi = rng.uniform(-1, 1, n); Om = rng.uniform(0, 2 * np.pi, n); w = rng.uniform(0, 2 * np.pi, n); tp = rng.uniform(0, 1, n) * P
    # Thiele-Innes
    sinI = np.sqrt(1 - cosi ** 2)
    TA = a_ph * (np.cos(w) * np.cos(Om) - np.sin(w) * np.sin(Om) * cosi)
    TB = a_ph * (np.cos(w) * np.sin(Om) + np.sin(w) * np.cos(Om) * cosi)
    TF = a_ph * (-np.sin(w) * np.cos(Om) - np.cos(w) * np.sin(Om) * cosi)
    TG = a_ph * (-np.sin(w) * np.sin(Om) + np.cos(w) * np.cos(Om) * cosi)
    Mn = 2 * np.pi * ((ty[None, :] - tp[:, None]) / P[:, None])
    E = kepler_E(np.mod(Mn, 2 * np.pi), e[:, None])
    Xo = np.cos(E) - e[:, None]; Yo = np.sqrt(1 - e[:, None] ** 2) * np.sin(E)
    dde = TA[:, None] * Xo + TF[:, None] * Yo      # Delta delta (north), mas
    dra = TB[:, None] * Xo + TG[:, None] * Yo      # Delta alpha* (east), mas
    al = dra * sa[None, :] + dde * ca[None, :]
    noise = rng.normal(0, SIG, al.shape)
    y = al + noise   # parallax/pm/position of the barycentre add linearly and do not change residuals
    out = {}
    for name, m, H, R in [("dr3", m3, H3, R3), ("dr2", m2, H2, R2)]:
        yy = y[:, m]
        par = yy @ H.T
        res = yy @ R.T
        eps2 = np.sum(res ** 2, axis=1) / (m.sum() - 5) - SIG ** 2
        out[name] = dict(plx_bias=par[:, 4], pmra=par[:, 2], pmdec=par[:, 3], eps=np.sqrt(np.clip(eps2, 0, None)))
    return out, dict(P=P, q=q, a_ph=a_ph, e=e)


PLX_OBS3, E3 = 2.807 + 0.0248, 0.357
PLX_OBS2, E2 = 3.512 + 0.03, 0.735
EPS3 = 2.317
grid = np.array([0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.3, 2.6, 2.83, 3.1, 3.4, 3.8, 4.5])
res = {}
for dark in (False, True, "extreme"):
    key = {False: "luminous", True: "dark", "extreme": "dark_extreme_q3"}[dark]
    rows = []
    for pt in grid:
        o, prm = simulate(pt, 200000, dark)
        e3 = o["dr3"]["eps"]
        acc = (e3 > 0.8 * EPS3) & (e3 < 1.2 * EPS3)
        f_eps = acc.mean()
        pf3 = pt + o["dr3"]["plx_bias"][acc]; pf2 = pt + o["dr2"]["plx_bias"][acc]
        dpm = np.hypot(o["dr3"]["pmra"][acc] - o["dr2"]["pmra"][acc], o["dr3"]["pmdec"][acc] - o["dr2"]["pmdec"][acc])
        # Gaussian measurement kernel for DR3 parallax (width = quoted error)
        k3 = np.exp(-0.5 * ((pf3 - PLX_OBS3) / E3) ** 2) / (np.sqrt(2 * np.pi) * E3)
        k2 = np.exp(-0.5 * ((pf2 - PLX_OBS2) / E2) ** 2) / (np.sqrt(2 * np.pi) * E2)
        kpm = np.exp(-0.5 * ((dpm - 4.3) / 1.2) ** 2)
        L3 = f_eps * (k3.mean() if acc.sum() else 0)
        L32 = f_eps * ((k3 * k2).mean() if acc.sum() else 0)
        L32pm = f_eps * ((k3 * k2 * kpm).mean() if acc.sum() else 0)
        rows.append(dict(plx_true=float(pt), n_acc=int(acc.sum()), f_eps=float(f_eps), L_dr3=float(L3), L_dr3_dr2=float(L32), L_dr3_dr2_dpm=float(L32pm),
                         bias_med=float(np.median(pf3 - pt)) if acc.sum() else None,
                         bias_p2=float(np.percentile(pf3 - pt, 2.3)) if acc.sum() else None,
                         bias_p98=float(np.percentile(pf3 - pt, 97.7)) if acc.sum() else None,
                         frac_within_1sig=float(np.mean(np.abs(pf3 - PLX_OBS3) < E3)) if acc.sum() else None,
                         P_med_acc_yr=float(np.median(prm["P"][acc])) if acc.sum() else None,
                         dpm_med=float(np.median(dpm)) if acc.sum() else None))
        print(key, rows[-1], flush=True)
    res[key] = rows
# posterior over plx_true with a prior uniform in distance volume (p(d) ~ d^2 -> p(plx) ~ plx^-4) and flat in plx
for key in list(res.keys()):
    L = np.array([r["L_dr3_dr2"] for r in res[key]])
    for prior_name, pr in [("flat_plx", np.ones_like(grid)), ("volume", grid ** -4.0)]:
        post = L * pr
        if post.sum() > 0:
            post = post / np.trapezoid(post, grid)
            res[f"{key}_post_{prior_name}"] = dict(plx=grid.tolist(), post=post.tolist(),
                                                   P_plx_lt_1p3=float(np.trapezoid(post[grid <= 1.3], grid[grid <= 1.3]) if (grid <= 1.3).sum() > 1 else 0))
    i998 = np.argmin(np.abs(grid - 1.0)); i356 = np.argmin(np.abs(grid - 2.83)); i843 = np.argmin(np.abs(grid - 1.2))
    for lab in ("L_dr3", "L_dr3_dr2", "L_dr3_dr2_dpm"):
        Ls = np.array([r[lab] for r in res[key]])
        res[f"{key}_{lab}_ratio_998_over_356"] = float(Ls[i998] / Ls[i356]) if Ls[i356] > 0 else None
        res[f"{key}_{lab}_ratio_843_over_356"] = float(Ls[i843] / Ls[i356]) if Ls[i356] > 0 else None
json.dump(res, open(os.path.join(ROOT, "out/astrometric_sim.json"), "w"), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for key, c in [("luminous", "C0"), ("dark", "C1"), ("dark_extreme_q3", "C2")]:
    for lab, ls in [("L_dr3", ":"), ("L_dr3_dr2", "-"), ("L_dr3_dr2_dpm", "--")]:
        L = np.array([r[lab] for r in res[key]])
        if L.max() > 0: ax[0].plot(grid, L / L.max(), ls, color=c, label=f"{key} {lab}")
    ax[1].plot(grid, [r["f_eps"] for r in res[key]], "o-", color=c, label=f"{key}: P(eps_DR3 = 2.3 +- 20%)")
for a in ax: a.axvline(1.0, color="r", ls=":"); a.axvline(2.83, color="k", ls=":")
ax[0].set_xlabel("true parallax (mas)"); ax[0].set_ylabel("relative likelihood"); ax[0].legend(fontsize=6)
ax[1].set_xlabel("true parallax (mas)"); ax[1].set_yscale("log"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/astrometric_sim.png"), dpi=110)
print(json.dumps({k: v for k, v in res.items() if not isinstance(v, list)}, indent=1))
