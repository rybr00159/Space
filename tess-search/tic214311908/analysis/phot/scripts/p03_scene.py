"""Local static scene model around the Mira, per sector, on the median image (quality==0 cadences).

Positions: Gaia DR3 propagated to the sector epoch + TESScut WCS + global shift from p02 (normal-colour field stars).
The Mira gets a free position offset (its extreme colour shifts its TESS image) and optionally its own extra blur.
Free fluxes: Mira, M3V, pair (A+B combined, B/A fixed from Gaia G, at their flux-weighted position), and all other
stars with T < 15.5 within 6 px of the Mira; fainter stars enter with Gaia-predicted fluxes scaled by the p02 faint
scale (fixed); stars farther than 6 px are fixed at their p02 fluxes. Background: plane over the fit region.
Model variants compared by chi2 (weights 2% of flux + 5 e/s): (a) Mira at Gaia position; (b) Mira free position;
(c) Mira free position and free extra blur.
Outputs: out/scene_s{sector}.json (positions, fluxes, offsets, chi2 per variant), figs/scene_s{sector}.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RFIT = 4.5      # px around the Mira
TFREE_LOC = 15.5
RFREE = 6.0


def build_stars(C, W):
    """Star list with positions in the shifted WCS frame; pair A+B merged."""
    g = C.gaia(18.5)
    dx, dy = W["global_shift"]
    g["x"] = g.x + dx; g["y"] = g.y + dy
    iA = int(np.where(g.source_id == IDS["pairA"])[0][0]); iB = int(np.where(g.source_id == IDS["pairB"])[0][0])
    # pair B has a 2-parameter solution: propagate with A's proper motion
    dt = C.epoch - 2016.0
    raB = g.ra[iB] + g.pmra[iA] * dt / 3.6e6 / np.cos(np.radians(g.dec[iB])); deB = g.dec[iB] + g.pmdec[iA] * dt / 3.6e6
    xB, yB = C.wcs.all_world2pix([raB], [deB], 0); xB, yB = float(xB[0]) + dx, float(yB[0]) + dy
    rBA = 10 ** (-0.4 * (g.G[iB] - g.G[iA]))
    g.loc[iA, "x"] = (g.x[iA] + rBA * xB) / (1 + rBA); g.loc[iA, "y"] = (g.y[iA] + rBA * yB) / (1 + rBA)
    g.loc[iA, "T"] = g["T"][iA] - 2.5 * np.log10(1 + rBA)
    g = g.drop(index=iB).reset_index(drop=True)
    g["name"] = ""
    for k in ("Mira", "M3V", "pairA", "ref"):
        g.loc[g.source_id == IDS[k], "name"] = k
    g.loc[g.name == "pairA", "name"] = "pair"
    return g, rBA


def scene_setup(sector):
    C = Cutout(sector)
    W = json.load(open(f"{ROOT}/out/wcs_s{sector}.json"))
    g, rBA = build_stars(C, W)
    prf = PRFModel(C); prf.blur = W["blur"]
    im = g.index[g.name == "Mira"][0]
    xm, ym = g.x[im], g.y[im]
    dist = np.hypot(g.x - xm, g.y - ym)
    near = dist < RFIT + 7
    g = g[near].reset_index(drop=True); dist = np.hypot(g.x - xm, g.y - ym)
    free = [i for i in range(len(g)) if (g.name[i] in ("Mira", "M3V", "pair")) or (g["T"][i] < TFREE_LOC and dist[i] < RFREE)]
    ff = W["free_fluxes"]
    fixed_flux = np.array([ff.get(str(int(s)), tmag_flux(T) * W["faint_scale"]) for s, T in zip(g.source_id, g["T"])])
    return C, W, g, prf, free, fixed_flux, rBA


def fit_static(C, g, prf, free, fixed_flux, img, mask, variant, sig=None):
    """variant: 'gaia' (Mira at Gaia position), 'pos' (free Mira offset), 'posblur' (free offset and Mira blur)."""
    im = int(g.index[g.name == "Mira"][0])
    fixed_img = np.zeros(img.shape)
    for i in range(len(g)):
        if i not in free: fixed_img += fixed_flux[i] * prf.image(g.x[i], g.y[i])
    base_cols = {i: prf.image(g.x[i], g.y[i])[mask] for i in free if i != im}
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]; u = (xx[mask] - g.x[im]) / 5; v = (yy[mask] - g.y[im]) / 5
    b = (img - fixed_img)[mask]
    if sig is None: sig = 0.02 * np.abs(img[mask]) + 5.0
    def solve(p):
        mb = prf.blur if variant != "posblur" else np.hypot(prf.blur, p[2])
        cols = [prf.image(g.x[im] + p[0], g.y[im] + p[1], blur=mb)[mask]] + [base_cols[i] for i in free if i != im] + [np.ones(mask.sum()), u, v]
        A = np.column_stack(cols)
        c, *_ = np.linalg.lstsq(A / sig[:, None], b / sig, rcond=None)
        return A, c
    def res(p):
        A, c = solve(p); return (A @ c - b) / sig
    if variant == "gaia":
        p = np.array([0.0, 0.0]); A, c = solve(p); r = res(p); cov = None
    else:
        best = None
        for gx in np.linspace(-1.2, 1.2, 9):
            for gy in np.linspace(-1.2, 1.2, 9):
                p0 = [gx, gy] + ([0.2] if variant == "posblur" else [])
                c2 = np.sum(res(p0) ** 2)
                if best is None or c2 < best[0]: best = (c2, p0)
        q = least_squares(res, best[1], diff_step=1e-3, x_scale=0.05)
        p = q.x; A, c = solve(p); r = q.fun
        cov = np.linalg.pinv(q.jac.T @ q.jac) * (np.sum(r ** 2) / (mask.sum() - len(p) - A.shape[1]))
    order = [im] + [i for i in free if i != im]
    fluxes = {int(g.source_id[i]): float(c[k]) for k, i in enumerate(order)}
    model = np.full(img.shape, np.nan); model[mask] = A @ c + fixed_img[mask]
    return dict(p=np.array(p).tolist(), cov=None if cov is None else cov.tolist(), chi2=float(np.sum(r ** 2)), fluxes=fluxes,
                bkg=c[-3:].tolist(), model=model)


def run(sector):
    C, W, g, prf, free, fixed_flux, rBA = scene_setup(sector)
    good = C.fin & (C.q == 0)
    img = np.median(C.flux[good], 0)
    im = int(g.index[g.name == "Mira"][0])
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mask = np.hypot(xx - g.x[im], yy - g.y[im]) <= RFIT
    res = {}
    for v in ("gaia", "pos", "posblur"):
        res[v] = fit_static(C, g, prf, free, fixed_flux, img, mask, v)
    J = C.jacobian(g.x[im], g.y[im])
    out = dict(sector=sector, npix=int(mask.sum()), rfit=RFIT, wcs_shift=W["global_shift"], blur=W["blur"], pairB_over_A=float(rBA),
               stars=[dict(source_id=int(g.source_id[i]), name=g.name[i], T=float(g["T"][i]), x=float(g.x[i]), y=float(g.y[i]),
                           free=bool(i in free), fixed_flux=float(fixed_flux[i])) for i in range(len(g))])
    for v, r in res.items():
        d = dict(p=r["p"], cov=r["cov"], chi2=r["chi2"], bkg=r["bkg"],
                 flux={nm: r["fluxes"][int(g.source_id[g.name == nm].iloc[0])] for nm in ("Mira", "M3V", "pair")},
                 flux_all={str(k): vv for k, vv in r["fluxes"].items()})
        if v != "gaia":
            d["mira_offset_EN_arcsec"] = (J @ np.array(r["p"][:2])).tolist()
            d["mira_offset_err_px"] = np.sqrt(np.diag(r["cov"]))[:2].tolist()
        out[v] = d
    exp = {nm: tmag_flux(g["T"][g.name == nm].iloc[0]) for nm in ("M3V", "pair")}
    out["gaia_expected_flux"] = exp
    save_json(out, f"{ROOT}/out/scene_s{sector}.json")
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.5))
    ext = None
    ax[0].imshow(img, origin="lower", vmin=np.percentile(img, 5), vmax=np.percentile(img, 99.5)); ax[0].set_title(f"S{sector} median")
    for k, v in enumerate(("gaia", "pos", "posblur")):
        R = (img - res[v]["model"]) / np.sqrt(np.abs(img))
        vm = 3 * np.nanstd(R[mask])
        ax[k + 1].imshow(np.where(mask, R, np.nan), origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm)
        ax[k + 1].set_title(f"{v}: chi2 {res[v]['chi2']:.0f}; M3V {res[v]['fluxes'][IDS['M3V']]:.0f}, pair {res[v]['fluxes'][IDS['pairA']]:.0f}")
    for a in ax:
        for nm, c in (("Mira", "r"), ("M3V", "y"), ("pair", "m")):
            r_ = g[g.name == nm].iloc[0]; a.plot(r_.x, r_.y, "+", color=c, ms=10, mew=2)
        p = res["posblur"]["p"]; a.plot(g.x[im] + p[0], g.y[im] + p[1], "x", color="k", ms=10, mew=2)
        a.set_xlim(g.x[im] - 5.5, g.x[im] + 5.5); a.set_ylim(g.y[im] - 5.5, g.y[im] + 5.5)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/scene_s{sector}.png", dpi=75); plt.close(fig)
    msg = f"S{sector}: "
    for v in ("gaia", "pos", "posblur"):
        r = out[v]; msg += f"[{v}] chi2 {r['chi2']:.0f} Mira {r['flux']['Mira']:.0f} M3V {r['flux']['M3V']:.0f} pair {r['flux']['pair']:.0f}"
        if v != "gaia": msg += f" off(px) {np.round(r['p'], 3).tolist()} EN\" {np.round(r['mira_offset_EN_arcsec'], 1).tolist()}"
        msg += "; "
    print(msg + f"expected M3V {exp['M3V']:.0f} pair {exp['pair']:.0f}", flush=True)


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
