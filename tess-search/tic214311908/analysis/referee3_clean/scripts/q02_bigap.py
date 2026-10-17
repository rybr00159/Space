"""q02: PRF-independent per-transit flux loss from LARGE apertures in the raw TESScut cubes (referee3, own code).

Idea: the M3V, the pair and the Mira lie within ~0.9 px of each other. A large aperture (r = 2.5-4 px) holds
essentially all of the M3V's light, so the flux loss it measures depends on the PSF only through the encircled
energy EE (0.93-0.99), not through a PRF fit in a scene dominated by the Mira. Field-wide estimators are not used
except a background level from faint pixels (variables, bright stars and the Mira region excluded).

Per sector:
  * time: TESScut TIME + (Roemer(target) - TIMECORR) from q01 (TESS ephemeris from JPL Horizons), i.e. BTJD_TDB.
  * cadences: QUALITY == 0, finite.
  * background b(t): per-cadence median of 'faint' pixels (> 6 px from the Mira, > 2.5 px from any star with
    T < 14, > 3 px from any Gaia variable, median flux in the lower half of the remaining pixels).
  * pointing (optional regressor): mean flux-weighted centroid offset of isolated bright stars (T < 11.8, > 5 px
    from the Mira, not variable, no T < T+2 neighbour within 2 px, > 2.5 px from the edges).
  * apertures (pixel centres within r of the centre): 'm25','m30','m35' centred on the M3V (r = 2.5, 3.0, 3.5 px);
    'g40' centred on the Mira/pair/M3V light centre with r = 4.0; 'box3' = 3x3 box at the M3V pixel.
  * A(t) = sum(flux in aperture) - Npix * b(t).
  * per window (|t - tc| < 0.30 d): A = c0 + c1 tau + c2 tau^2 [+ cx dx + cy dy] - L * s(t), s = batman shape
    normalised to 1 at mid-transit (k 0.117, b 0.52, a/R* 25.3, u 0.125/0.533; exposure-integrated, 5/9/21 sub-steps),
    L = flux loss at mid-transit in e-/s. Acceptance: >= 2 points at |tau| < 0.028 d and >= 25 % of the out-of-
    transit cadences on each side (as the analyst).
  * fake epochs: every window centre at phase 0.12-0.88 (step 0.02, ~38 per real transit slot) gives the empirical
    per-window scatter of L; per-transit errors = sector fake std; t-distribution not assumed.
  * depth_ppt = L / (979 * EE) * 1000, EE = sum of the photlib SPOC PRF (13x13-normalised, i.e. 'PRF units', the
    units of the 979 e-/s normalisation) over the aperture pixels, at the M3V position.
Ephemeris: final combined T0 = 3840.91488, P = 5.6915020 (clean/out/ephemeris_final.json).
Variants via environment: Q2_BKG (faint|none|wide), Q2_ORDER (baseline polynomial order, default 2), Q2_HALF
(window half-width, default 0.30 d), Q2_QUAL (q0 | q0_2048), Q2_TAG (output suffix), Q2_FAKESTEP (fake phase step),
Q2_PHRANGE ('start,stop,step' phase grid replacing the fake grid, e.g. for a secondary-eclipse scan).
Output: out/q02_bigap_per_transit{tag}.csv, out/q02_bigap_fake{tag}.csv, out/q02_bigap{tag}.json
"""
import sys, glob, json
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import batman

R3 = "tic214311908/analysis/referee3_clean"
PH = "tic214311908/analysis/phot"
SR = "tic214311908/analysis/search"
sys.path.insert(0, f"{PH}/scripts")
T0, P = 3840.91488, 5.6915020
F_M3V = 979.0
K, B_, AR, U = 0.117, 0.52, 25.3, (0.125, 0.533)
IDS = dict(Mira=5925373268528684288, M3V=5925373272850784768, pairA=5925373272850782208, pairB=5925373272827253632)
SECTORS = (12, 39, 66, 93, 103, 104)
HALF = 0.30
import os
OPT = dict(bkg=os.environ.get("Q2_BKG", "faint"), order=int(os.environ.get("Q2_ORDER", "2")), half=float(os.environ.get("Q2_HALF", "0.30")),
           qual=os.environ.get("Q2_QUAL", "q0"), tag=os.environ.get("Q2_TAG", ""), fakestep=float(os.environ.get("Q2_FAKESTEP", "0.02")))
HALF = OPT["half"]


def shape(t, tc, ex):
    ns = 5 if ex * 86400 < 300 else (9 if ex * 86400 < 900 else 21)
    off = ((np.arange(ns) + 0.5) / ns - 0.5) * ex
    tt = (t[:, None] + off[None]).ravel()
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = tc, P, K, AR, np.degrees(np.arccos(B_ / AR)), 0., 90., list(U), "quadratic"
    f = batman.TransitModel(p, tt).light_curve(p).reshape(-1, ns).mean(1)
    p0 = batman.TransitModel(p, np.array([tc])).light_curve(p)[0]
    return (1 - f) / (1 - p0)


def fit_window(t, A, X, tc, ex):
    tau = t - tc
    s = shape(t, tc, ex)
    cols = [-s] + [tau ** o for o in range(OPT["order"] + 1)] + ([X[:, j] for j in range(X.shape[1])] if X is not None else [])
    M = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(M, A, rcond=None)
    r = A - M @ c
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    cov = np.linalg.pinv(M.T @ M) * sig ** 2
    return c[0], np.sqrt(cov[0, 0]), sig


def accept(t, tc, ex):
    tau = t - tc
    oot = np.abs(tau) > 0.07
    nexp = 2 * (HALF - 0.07) / ex
    return (np.sum(oot & (tau < 0)) >= 0.25 * nexp / 2) and (np.sum(oot & (tau > 0)) >= 0.25 * nexp / 2) and np.sum(np.abs(tau) < 0.028) >= 2


def main():
    gaia = pd.read_csv(f"{PH}/data/gaia/gaia_360arcsec.csv")
    var = json.load(open(f"{SR}/out/gaia_variables.json"))["variables"]
    var_ids = set(int(v["source_id"]) for v in var) - {IDS["Mira"]}
    addtab = np.load(f"{R3}/out/q01_addtab.npz")
    from photlib import Cutout, PRFModel
    rows, frows, summ = [], [], {}
    for sec in SECTORS:
        fn = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{sec:04d}-*_astrocut.fits")[0]
        with fits.open(fn) as h:
            d = h[1].data; w = WCS(h[2].header)
            t = np.array(d["TIME"], float); F = np.array(d["FLUX"], np.float64); q = np.array(d["QUALITY"], int)
        qsel = (q == 0) if OPT["qual"] == "q0" else ((q == 0) | (q == 2048))
        ok = qsel & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
        t, F = t[ok], F[ok]
        t = t + np.interp(t, addtab[f"t_{sec}"], addtab[f"add_{sec}"]) / 86400.0
        ex = float(np.median(np.diff(t)))
        ny, nx = F.shape[1:]
        ep = 2000.0 + (np.median(t) + 2457000 - 2451545) / 365.25
        g = gaia.copy(); dt = ep - 2016.0
        g["ra_e"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
        g["dec_e"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
        g["x"], g["y"] = w.all_world2pix(g.ra_e.values, g.dec_e.values, 0)
        pos = {k: g[g.source_id == v][["x", "y"]].values[0] for k, v in IDS.items()}
        yy, xx = np.mgrid[0:ny, 0:nx]
        dmira = np.hypot(xx - pos["Mira"][0], yy - pos["Mira"][1])
        # ---------- background pixels
        bad = dmira < 6
        for _, s in g[(g["T"] < 14)].iterrows():
            bad |= np.hypot(xx - s.x, yy - s.y) < 2.5
        for _, s in g[g.source_id.isin(var_ids)].iterrows():
            bad |= np.hypot(xx - s.x, yy - s.y) < 3.0
        med = np.median(F, 0)
        cand = ~bad
        thr = np.median(med[cand])
        bpix = cand & (med <= thr)
        b = np.median(F[:, bpix], axis=1)
        if OPT["bkg"] == "none":
            b = np.zeros(len(t))
        elif OPT["bkg"] == "wide":       # all pixels > 6 px from the Mira not within 1.5 px of a T < 12 star or 3 px of a variable
            bad2 = dmira < 6
            for _, s_ in g[(g["T"] < 12)].iterrows():
                bad2 |= np.hypot(xx - s_.x, yy - s_.y) < 1.5
            for _, s_ in g[g.source_id.isin(var_ids)].iterrows():
                bad2 |= np.hypot(xx - s_.x, yy - s_.y) < 3.0
            b = np.median(F[:, ~bad2], axis=1); bpix = ~bad2
        # ---------- pointing from isolated bright stars
        offs = []
        for _, s in g[(g["T"] < 11.8)].iterrows():
            if s.source_id in var_ids or s.source_id == IDS["Mira"]: continue
            if not (2.5 < s.x < nx - 3.5 and 2.5 < s.y < ny - 3.5): continue
            if np.hypot(s.x - pos["Mira"][0], s.y - pos["Mira"][1]) < 5: continue
            nb = g[(np.hypot(g.x - s.x, g.y - s.y) < 2.0) & (g["T"] < s["T"] + 2) & (g.source_id != s.source_id)]
            if len(nb): continue
            ix, iy = int(round(s.x)), int(round(s.y))
            sub = F[:, iy - 2:iy + 3, ix - 2:ix + 3] - b[:, None, None]
            tot = sub.sum((1, 2))
            cx = (sub * (np.arange(-2, 3)[None, None, :])).sum((1, 2)) / tot
            cy = (sub * (np.arange(-2, 3)[None, :, None])).sum((1, 2)) / tot
            offs.append((cx - np.median(cx), cy - np.median(cy), s["T"]))
        if offs:
            wts = np.array([10 ** (-0.4 * (o[2] - 10)) for o in offs])
            dxp = np.sum([o[0] * wi for o, wi in zip(offs, wts)], 0) / wts.sum()
            dyp = np.sum([o[1] * wi for o, wi in zip(offs, wts)], 0) / wts.sum()
        else:
            dxp = dyp = np.zeros(len(t))
        XP = np.column_stack([dxp, dyp])
        # ---------- apertures and EE
        C = Cutout(sec)
        prf = PRFModel(C)
        m3 = pos["M3V"]
        img = prf.image(m3[0], m3[1])
        # light centre of the group from the PRF scene with TIC-like fluxes (Mira ~ median aperture)
        Fm = {"Mira": 15000.0, "pairA": 655.0, "pairB": 240.0, "M3V": 979.0}
        gc = sum(Fm[k] * pos[k] for k in Fm) / sum(Fm.values())
        aps = {}
        for r_ in (2.5, 3.0, 3.5):
            aps[f"m{int(r_ * 10)}"] = np.hypot(xx - m3[0], yy - m3[1]) <= r_
        aps["g40"] = np.hypot(xx - gc[0], yy - gc[1]) <= 4.0
        aps["g50"] = np.hypot(xx - gc[0], yy - gc[1]) <= 5.0
        ix, iy = int(round(m3[0])), int(round(m3[1]))
        bx = np.zeros((ny, nx), bool); bx[iy - 1:iy + 2, ix - 1:ix + 2] = True
        aps["box3"] = bx
        EE = {k: float(img[a].sum()) for k, a in aps.items()}
        # sensitivity of EE: +-0.1 px shifts and 0.15 px blur
        EEs = {k: [float(prf.image(m3[0] + dx_, m3[1] + dy_)[a].sum()) for dx_, dy_ in ((0.1, 0), (-0.1, 0), (0, 0.1), (0, -0.1))]
               + [float(prf.image(m3[0], m3[1], blur=0.15)[a].sum())] for k, a in aps.items()}
        Aser = {k: F[:, a].sum(1) - a.sum() * b for k, a in aps.items()}
        # ---------- real and fake windows
        n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
        for e in range(n0, n1 + 1):
            phr = os.environ.get("Q2_PHRANGE")
            grid = np.arange(*map(float, phr.split(","))) if phr else np.arange(0.12, 0.881, OPT["fakestep"])
            for fph in np.r_[0.0, grid]:
                tc = T0 + (e + fph) * P
                m = np.abs(t - tc) < HALF
                if m.sum() < 8 or not accept(t[m], tc, ex): continue
                rec = dict(sector=sec, epoch=e, phase=round(fph, 3), tc=tc, n=int(m.sum()))
                for k, A in Aser.items():
                    for reg, X in (("noreg", None), ("pt", XP[m])):
                        L, Le, sg = fit_window(t[m], A[m], X, tc, ex)
                        rec[f"L_{k}_{reg}"] = L; rec[f"Lerr_{k}_{reg}"] = Le; rec[f"sig_{k}_{reg}"] = sg
                (rows if fph == 0.0 else frows).append(rec)
        summ[sec] = dict(npix_bkg=int(bpix.sum()), n_pointing_stars=len(offs), exptime_s=ex * 86400, m3v_xy=list(map(float, m3)),
                         mira_xy=list(map(float, pos["Mira"])), group_centre=list(map(float, gc)), EE=EE, EE_variants=EEs,
                         npix_ap={k: int(a.sum()) for k, a in aps.items()}, pointing_rms_px=[float(np.std(dxp)), float(np.std(dyp))])
        print(sec, summ[sec]["EE"], summ[sec]["n_pointing_stars"], summ[sec]["npix_bkg"], flush=True)
    tg = OPT["tag"]
    summ["options"] = OPT
    pd.DataFrame(rows).to_csv(f"{R3}/out/q02_bigap_per_transit{tg}.csv", index=False)
    pd.DataFrame(frows).to_csv(f"{R3}/out/q02_bigap_fake{tg}.csv", index=False)
    json.dump(summ, open(f"{R3}/out/q02_bigap{tg}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
