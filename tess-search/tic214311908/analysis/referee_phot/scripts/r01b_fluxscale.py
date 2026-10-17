"""Referee r01b: PRF flux scale from a global fit of the field. All Gaia stars with T < 15 get free fluxes (positions
fixed: Gaia propagated + r01 offset), fainter stars (T < 18) enter with catalogue fluxes times one free scale, plus a
quadratic background; pixels within 4 px of the Mira masked. Ratio fitted flux / TIC-catalogue flux (15000 e/s at
T=10) for stars with T < 14 and no neighbour within 1 px brighter than T+1. Output: out/fluxscale.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
tic = pd.read_csv(f"{R}/data/tic330.csv", dtype={"GAIA": str, "ID": str}); tic = tic[tic.GAIA.notna()]
TT = dict(zip(tic.GAIA.astype(np.int64), tic.Tmag))
out = {}
allr = []
for s in SECTORS:
    C = Cut(s); g = C.stars(18.0); W = json.load(open(f"{R}/out/wcs_s{s}.json")); ox, oy = W["offset_px"]
    g["x"] += ox; g["y"] += oy
    good = C.fin & (C.q == 0); fl = C.f.reshape(C.n, -1); bk = np.nanpercentile(fl, 10, axis=1)
    good &= bk < np.nanmedian(bk[good]) + 3 * rstd(bk[good])
    img = np.median(C.f[good], 0)
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mira = g[g.source_id == GID["Mira"]].iloc[0]
    mask = (np.hypot(xx - mira.x, yy - mira.y) > 4.0)
    ins = g[(g.x > -4) & (g.x < C.nx + 3) & (g.y > -4) & (g.y < C.ny + 3)].reset_index(drop=True)
    free = ins[(ins["T"] < 15.0) & (np.hypot(ins.x - mira.x, ins.y - mira.y) > 2.5)]
    faint = ins[(ins["T"] >= 15.0) & (np.hypot(ins.x - mira.x, ins.y - mira.y) > 2.5)]
    cols = [C.psf(r.x, r.y)[mask] for r in free.itertuples()]
    fa = np.zeros(mask.sum())
    for r in faint.itertuples(): fa += tflux(r.T) * C.psf(r.x, r.y)[mask]
    u = (xx[mask] - 10) / 10; v = (yy[mask] - 10) / 10
    A = np.column_stack(cols + [fa, np.ones(mask.sum()), u, v, u * u, v * v, u * v])
    w = 1 / np.sqrt(np.abs(img[mask]) + 100)
    c, *_ = np.linalg.lstsq(A * w[:, None], img[mask] * w, rcond=None)
    res = img[mask] - A @ c
    f = c[:len(free)]
    free = free.assign(fit=f, Ttic=free.source_id.map(TT))
    rows = []
    for r in free.itertuples():
        d = np.hypot(ins.x - r.x, ins.y - r.y)
        iso = not np.any((d > 0) & (d < 1.0) & (ins["T"] < r.T + 1.0))
        inside = 1.5 < r.x < C.nx - 2.5 and 1.5 < r.y < C.ny - 2.5
        rows.append(dict(source_id=int(r.source_id), T=float(r.T), Ttic=float(r.Ttic), bp_rp=float(r.bp_rp), fit=float(r.fit),
                         ratio=float(r.fit / tflux(r.Ttic)), iso=bool(iso), inside=bool(inside)))
    d = pd.DataFrame(rows); d["sector"] = s; allr.append(d)
    sel = d[d.iso & d.inside & (d.Ttic < 14.0)]
    nrm = sel[sel.bp_rp < 2.0]; red = sel[sel.bp_rp >= 2.0]
    out[s] = dict(n=len(sel), ratio_all=(float(np.median(sel.ratio)), float(rstd(sel.ratio))), n_normal=len(nrm),
                  ratio_normal=(float(np.median(nrm.ratio)), float(rstd(nrm.ratio))), n_red=len(red),
                  ratio_red=(float(np.median(red.ratio)) if len(red) else np.nan), faint_scale=float(c[len(free)]),
                  resid_rel=float(rstd(res / np.clip(img[mask], 50, None))))
    print(f"S{s}: n {len(sel)}, PRF-fit/TIC flux median {out[s]['ratio_all'][0]:.3f} (rstd {out[s]['ratio_all'][1]:.3f}); normal {out[s]['ratio_normal'][0]:.3f} (n {len(nrm)}), "
          f"red {out[s]['ratio_red']:.3f} (n {len(red)}; {np.round(red[['bp_rp','ratio']].values,3).tolist()}); faint scale {c[len(free)]:.2f}; resid rel {out[s]['resid_rel']:.3f}", flush=True)
D = pd.concat(allr); D.to_csv(f"{R}/out/fluxscale_stars.csv", index=False)
sel = D[D.iso & D.inside & (D.Ttic < 14.0)]
# per-sector median ratio and its uncertainty
summ = {str(s): dict(median=float(np.median(sel[sel.sector == s].ratio)), err=float(1.2533 * rstd(sel[sel.sector == s].ratio) / np.sqrt((sel.sector == s).sum())))
        for s in SECTORS}
out["summary"] = summ
out["overall_median"] = float(np.median([v["median"] for v in summ.values()]))
print(summ, out["overall_median"])
save_json({str(k): v for k, v in out.items()}, f"{R}/out/fluxscale.json")
