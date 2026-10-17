"""Referee r04b: static position of the Mira in the median image. Model within 3.5 px of the Mira: Mira = SPOC PRF with
free position, flux and extra Gaussian blur (5x5 shifted-PRF average); M3V, pair A, pair B and all other stars
(T < 18) at Gaia positions + r01 offset with catalogue fluxes x the r01b PRF flux scale (fixed); plane background.
Variant: M-dwarf fluxes free. Output: out/mirastatic.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from scipy.optimize import least_squares
FS = json.load(open(f"{R}/out/fluxscale.json"))["summary"]
tic = pd.read_csv(f"{R}/data/tic330.csv", dtype={"GAIA": str}); tic = tic[tic.GAIA.notna()]; TT = dict(zip(tic.GAIA.astype(np.int64), tic.Tmag))
out = {}
for s in SECTORS:
    C = Cut(s); g = C.stars(18.0); W = json.load(open(f"{R}/out/wcs_s{s}.json")); ox, oy = W["offset_px"]; g["x"] += ox; g["y"] += oy
    sc = FS[str(s)]["median"]
    a = np.load(f"{R}/data/aux_s{s}.npz")
    ok = C.fin & (C.q == 0) & (a["B"] < np.nanmedian(a["B"][C.q == 0]) + 60)
    img = np.median(C.f[ok] - a["B"][ok, None, None], 0)
    mira = g[g.source_id == GID["Mira"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]; m = np.hypot(xx - mira.x, yy - mira.y) <= 3.5
    fixed = np.zeros((C.ny, C.nx)); mdw = {}
    for r in g.itertuples():
        if r.source_id == GID["Mira"] or np.hypot(r.x - mira.x, r.y - mira.y) > 8: continue
        T_ = TT.get(r.source_id, r.T) if r.source_id != GID["pB"] else r.T
        Pp = C.psf(r.x, r.y)
        if r.source_id in (GID["M3V"], GID["pA"], GID["pB"]): mdw[r.source_id] = (Pp, tflux(T_) * sc)
        else: fixed += tflux(T_) * sc * Pp
    u = (xx[m] - mira.x) / 3; v = (yy[m] - mira.y) / 3
    w = 1 / np.sqrt(np.abs(img[m]) + 200 + (0.01 * img[m]) ** 2)
    def psfb(x, y, bl):
        bl = abs(bl)
        if bl < 0.01: return C.psf(x, y)
        n_ = np.array([-2, -1, 0, 1, 2]) * bl * 0.8; ww = np.exp(-0.5 * (n_ / bl) ** 2); ww /= ww.sum()
        return sum(ww[i] * ww[j] * C.psf(x + n_[i], y + n_[j]) for i in range(5) for j in range(5))
    for var in ("mdw_fixed", "mdw_free"):
        def model(p):
            cols = [psfb(p[0], p[1], p[2])[m], np.ones(m.sum()), u, v]
            dat = img[m] - fixed[m]
            if var == "mdw_fixed":
                dat = dat - sum(f_ * P_[m] for P_, f_ in mdw.values())
            else:
                cols += [mdw[GID["M3V"]][0][m], mdw[GID["pA"]][0][m] + 0.367 * mdw[GID["pB"]][0][m]]
            A = np.column_stack(cols); c, *_ = np.linalg.lstsq(A * w[:, None], dat * w, rcond=None)
            return A, c, dat
        def rr(p):
            A, c, dat = model(p); return (A @ c - dat) * w
        best = min(((np.sum(rr([mira.x + dx_, mira.y + dy_, 0.2]) ** 2), [mira.x + dx_, mira.y + dy_, 0.2]) for dx_ in np.linspace(-1, 1, 9) for dy_ in np.linspace(-1, 1, 9)), key=lambda z: z[0])
        q = least_squares(rr, best[1], diff_step=1e-3, x_scale=[0.05, 0.05, 0.05])
        A, c, dat = model(q.x)
        off = C.en(mira.x, mira.y, q.x[0] - mira.x, q.x[1] - mira.y)
        out.setdefault(s, {})[var] = dict(offset_EN=off.tolist(), blur=float(abs(q.x[2])), mira_flux=float(c[0]), chi2=float(2 * q.cost), npix=int(m.sum()),
                                          mdw_fluxes=(c[4:].tolist() if var == "mdw_free" else None),
                                          expected=dict(M3V=mdw[GID["M3V"]][1], pair=mdw[GID["pA"]][1] + mdw[GID["pB"]][1]))
    o = out[s]
    print(f"S{s}: Mira static PRF offset (E,N) fixed-Mdw {np.round(o['mdw_fixed']['offset_EN'],1).tolist()}\" blur {o['mdw_fixed']['blur']:.2f} chi2 {o['mdw_fixed']['chi2']:.0f}; free-Mdw {np.round(o['mdw_free']['offset_EN'],1).tolist()}\" "
          f"chi2 {o['mdw_free']['chi2']:.0f}, M3V {o['mdw_free']['mdw_fluxes'][0]:.0f} (exp {o['mdw_free']['expected']['M3V']:.0f}), pair {o['mdw_free']['mdw_fluxes'][1]:.0f} (exp {o['mdw_free']['expected']['pair']:.0f}); Mira flux {o['mdw_fixed']['mira_flux']:.0f}", flush=True)
save_json({str(k): v for k, v in out.items()}, f"{R}/out/mirastatic.json")
