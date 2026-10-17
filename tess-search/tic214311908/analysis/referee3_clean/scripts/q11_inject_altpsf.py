"""q11: does the analyst's INV (cleaned, mode 'all') and REF pipeline respond to a constant transit in a
SECTOR-DEPENDENT way when the injected PSF is NOT the fitting template?  (The analyst's matched-template injections
recover 1.000 by construction: both pipelines are linear least squares with the injected image in the design matrix.)

Injected images (unit flux, at the M3V position; the fitting templates are never changed):
  'mirab'  : the empirical Mira image b (d pixel / d F_Mira from the analyst's own no-injection INV regression, pixels
             within 4 px of the Mira), shifted to the M3V position with a cubic-spline sub-pixel shift and scaled to
             the same total as the SPOC PRF model over the same pixels -> the detector's real PSF ~1 px away
             (but for a much redder star);
  'altpos' : the SPOC PRF evaluated 400 detector pixels away (towards the CCD centre) at the same sub-pixel phase
             -> size of PRF-model spatial-interpolation error.
Amplitude 15.37 e-/s at mid-transit (15.7 ppt of 979 e-/s), analyst's batman shape, 10 fake ephemerides per sector
(phase offsets 0.13-0.87, as c05), paired no-injection subtraction.  INV: cleanlib.photometry + make_lc + c05.fit_window
on raw (offset+slope); REF: refclean.SectorC windows + c05.fit_window.
The analyst's code is imported read-only; all outputs go to referee3_clean/out.
Output: out/q11_inject_s{sector}.csv, out/q11_summary.json (via --summary)
"""
import sys, json
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
import numpy as np, pandas as pd
from scipy.ndimage import shift as ndshift
R3 = "tic214311908/analysis/referee3_clean"


def run(sector):
    from cleanlib import setup, photometry, make_lc, in_windows, t_tdb_of, F_M3V, P1, T01
    from refclean import SectorC
    import c05_inject as c5
    from photlib import PRFModel
    S = setup(sector, "all"); C = S["C"]
    t_tdb = t_tdb_of(sector, C.t, C.timecorr)
    g = S["g"]; m3 = g[g.name == "M3V"].iloc[0]; mi = g[g.name == "Mira"].iloc[0]
    prf = S["prf"]
    R0 = photometry(sector, S, fit_models=("B",))
    b = R0["b_img"].copy()
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    near = np.hypot(xx - mi.x, yy - mi.y) <= 4
    b[~near] = 0.0
    bs = ndshift(b, (m3.y - mi.y, m3.x - mi.x), order=3, mode="constant")
    nearm = np.hypot(xx - m3.x, yy - m3.y) <= 4
    bs[~nearm] = 0.0
    P0 = prf.image(m3.x, m3.y)
    img_b = bs * (P0[nearm].sum() / bs[nearm].sum())
    # alternative-position PRF
    class Cshift: pass
    C2 = Cshift(); C2.sector = C.sector; C2.camera = C.camera; C2.ccd = C.ccd; C2.ny, C2.nx = C.ny, C.nx
    dc = 400 if C.col0 < 1000 else -400; dr = 400 if C.row0 < 1000 else -400
    C2.col0 = C.col0 + dc; C2.row0 = C.row0 + dr
    prf2 = PRFModel(C2); prf2.blur = prf.blur
    img_alt = prf2.image(m3.x, m3.y)
    imgs = {"mirab": img_b, "altpos": img_alt}
    diag = {k: dict(corr_with_template=float(np.corrcoef(v[nearm], P0[nearm])[0, 1]),
                    l2_frac_diff=float(np.linalg.norm((v - P0)[nearm]) / np.linalg.norm(P0[nearm])),
                    projection_on_template=float((v[nearm] @ P0[nearm]) / (P0[nearm] @ P0[nearm]))) for k, v in imgs.items()}
    print(sector, diag, flush=True)
    R_ = SectorC(sector, "all"); flat0 = R_.flat.copy()
    rows = []
    tfin = np.isfinite(t_tdb)
    for f in c5.FAKE_F:
        T0f = T01 + f * P1
        k = np.arange(np.floor((np.nanmin(t_tdb) - T0f) / P1) - 1, np.ceil((np.nanmax(t_tdb) - T0f) / P1) + 2)
        fake_ep = T0f + k * P1
        fake_ep = fake_ep[(fake_ep > np.nanmin(t_tdb) - 0.2) & (fake_ep < np.nanmax(t_tdb) + 0.2)]
        amp = np.zeros(len(t_tdb)); amp[tfin] = c5.AMP * c5.shape_frac(t_tdb[tfin], T0f, C.exptime)
        ampR = np.zeros(len(R_.ttdb)); fr = np.isfinite(R_.ttdb); ampR[fr] = c5.AMP * c5.shape_frac(R_.ttdb[fr], T0f, R_.exptime)
        res = {}
        for var in ("none",) + tuple(imgs):
            inj = None if var == "none" else [dict(img=imgs[var], amp=amp)]
            R = photometry(sector, S, inject=inj, extra_mask_times=fake_ep, fit_models=("B",))
            extra = in_windows(t_tdb, fake_ep, 0.1)
            L = make_lc(sector, R, extra_mask=extra, models=(("B", "M3V"),))
            ok = L["quality_ok"] & np.isfinite(L["raw_M3V_B"])
            R_.flat = flat0 if var == "none" else flat0 - ampR[:, None] * imgs[var].ravel()[None]
            for j, tc in enumerate(fake_ep):
                out = {}
                m = ok & (np.abs(L["t_tdb"] - tc) < c5.HALF)
                if m.sum() >= 8 and c5.accept(L["t_tdb"][m], tc, c5.HALF, C.exptime):
                    r = c5.fit_window(L["t_tdb"][m], L["raw_M3V_B"][m], tc, C.exptime, free=False)
                    if r: out["inv"] = r["depth_fix"]
                w = R_.window(tc)
                if w is not None:
                    tt = w["t"]; ff = 1 + w["aM"] / F_M3V
                    if c5.accept(tt, tc, {12: 0.6, 39: 0.4}.get(sector, 0.35), R_.exptime):
                        r = c5.fit_window(tt, ff, tc, R_.exptime, free=False)
                        if r: out["ref"] = r["depth_fix"]
                res[(j, var)] = out
        R_.flat = flat0
        for j, tc in enumerate(fake_ep):
            base = res[(j, "none")]
            for var in imgs:
                o = res[(j, var)]
                row = dict(sector=sector, fake_phase=f, tc=tc, variant=var)
                for p in ("inv", "ref"):
                    if p in o and p in base: row[f"{p}_resp"] = (o[p] - base[p]) / c5.DEP_PPT
                rows.append(row)
    pd.DataFrame(rows).to_csv(f"{R3}/out/q11_inject_s{sector}.csv", index=False)
    json.dump(diag, open(f"{R3}/out/q11_diag_s{sector}.json", "w"), indent=1)


def summary():
    import glob
    out = {}
    for fn in sorted(glob.glob(f"{R3}/out/q11_inject_s*.csv")):
        d = pd.read_csv(fn); s = int(d.sector.iloc[0])
        dg = json.load(open(f"{R3}/out/q11_diag_s{s}.json"))
        out[s] = {}
        for var, g in d.groupby("variant"):
            out[s][var] = {p: dict(mean=float(g[f"{p}_resp"].mean()), sem=float(g[f"{p}_resp"].std() / np.sqrt(g[f"{p}_resp"].notna().sum())),
                                   n=int(g[f"{p}_resp"].notna().sum())) for p in ("inv", "ref") if f"{p}_resp" in g}
            out[s][var]["diag"] = dg[var]
    for var in ("mirab", "altpos"):
        for p in ("inv", "ref"):
            v = {s: out[s][var][p]["mean"] for s in out if var in out[s] and p in out[s][var]}
            if all(k in v for k in (39, 66, 103, 104)):
                ratio = (v[103] + v[104]) / (v[39] + v[66])
                out.setdefault("late_early_ratio", {})[f"{var}_{p}"] = float(ratio)
            print(var, p, {k: round(x, 3) for k, x in v.items()})
    print(out.get("late_early_ratio"))
    json.dump(out, open(f"{R3}/out/q11_summary.json", "w"), indent=1)


if __name__ == "__main__":
    if "--summary" in sys.argv: summary()
    else:
        for s in [int(a) for a in sys.argv[1:]]:
            run(s)
