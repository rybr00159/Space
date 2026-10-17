"""s03: control-star light curves processed with the phot difference-image pipeline (p05 photometry, read-only import,
option `controls`: per-cadence PRF + plane fit to the same difference images within 2.5 px of each control star).

Controls: Gaia stars in the 21x21 cutout with T 11.5-15.0, > 3 px from the Mira, > 2.5 px from the cutout edge, no
neighbour within 1.5 px brighter than T + 1.5 (the same 3-5 stars per sector; the field is sparse).
Light curves: same cadence mask as the M3V (phot quality_ok), flux normalised by 1.107 x TIC flux (PRF units), wotan
biweight 0.75 d with the 5.69-d windows masked (identical to the M3V treatment).
Also 'phantom' light curves: a PRF fitted at empty positions at the M3V-Mira distance (0.9 px from the Mira) and
1.5-2 px from the M3V, jointly with the Mira empirical image, the M3V PRF and a plane (model-B-like), on the cached
difference images of s04 (crop 11x11). These share the Mira-induced systematics of the M3V light curve.
Output: data/ctrl_lcs.npz (dict of light curves: t, det, sector, cad, name), out/controls.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
sys.path.insert(0, "tic214311908/analysis/search/scripts")
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import photlib
from p05_phot import setup, photometry
import slib
from slib import detrend, transit_mask, robust_std, save_json
import numpy as np, json

OUT = slib.ROOT


def pick(C, g):
    mira = g[g.name == "Mira"].iloc[0]
    cand = []
    for r in g.itertuples():
        if not (11.5 <= r.T <= 15.0): continue
        if not (2.5 < r.x < C.nx - 3.5 and 2.5 < r.y < C.ny - 3.5): continue
        if np.hypot(r.x - mira.x, r.y - mira.y) < 3.0: continue
        d = np.hypot(g.x - r.x, g.y - r.y)
        if np.any((d > 0) & (d < 1.5) & (g["T"] < r.T + 1.5)): continue
        cand.append((int(r.source_id), float(r.x), float(r.y), float(r.T)))
    return cand


def phantoms(s):
    z = np.load(f"{OUT}/data/pix_s{s}.npz")
    names = list(z["pos_names"]); xy = z["pos_xy"]
    mira = xy[names.index("Mira")]; m3 = xy[names.index("M3V")]; pair = xy[names.index("pair")]
    r0 = np.hypot(*(m3 - mira)); a0 = np.arctan2(*(m3 - mira)[::-1])
    pts = []
    for ang in np.radians([90, 135, 180, 225, 270]):
        p = mira + r0 * np.array([np.cos(a0 + ang), np.sin(a0 + ang)])
        if np.hypot(*(p - pair)) > 0.45 and np.hypot(*(p - m3)) > 0.6: pts.append(p)
    for ang in np.radians([60, 120, 180, 240, 300]):
        p = m3 + 1.8 * np.array([np.cos(a0 + ang), np.sin(a0 + ang)])
        if np.hypot(*(p - pair)) > 0.6 and np.hypot(*(p - mira)) > 0.8: pts.append(p)
    return pts


def main():
    lcs = {}; meta = {}
    for s in photlib.SECTORS:
        S = setup(s); C, W, J, g, prf, rBA = S
        ctr = pick(C, g)
        R = photometry(s, S=S, verbose=False, controls={str(c[0]): (c[1], c[2]) for c in ctr})
        L = np.load(f"{photlib.ROOT}/lc_allstars_s{s:03d}.npz")
        t = L["t_tdb"]; ok = L["quality_ok"].astype(bool) & np.isfinite(t)
        cad = float(L["exptime_s"])
        for c in ctr:
            F = R["CTRL"][str(c[0])]
            fref = 1.107 * photlib.tmag_flux(c[3])
            raw = 1 + F / fref
            m = ok & np.isfinite(raw)
            tr = detrend(t[m], raw[m], transit_mask(t[m]))
            det = raw[m] / tr
            key = f"ctrl_{c[0]}"
            lcs.setdefault(key, []).append((t[m], det, np.full(m.sum(), s), np.full(m.sum(), cad)))
            meta.setdefault(key, dict(T=c[3], sectors=[]))["sectors"].append(s)
        # phantoms (cached difference images; model: Mira b + M3V PRF + phantom PRF + plane)
        z = np.load(f"{OUT}/data/pix_s{s}.npz")
        D = z["D"].astype(float); names = list(z["prf_names"]); b = z["b"]
        T3 = z["prf_imgs"][names.index("M3V")]
        x0, y0 = int(z["x0"]), int(z["y0"])
        ny, nx = D.shape[1:]; yy, xx = np.mgrid[0:ny, 0:nx]
        xym = z["pos_xy"][list(z["pos_names"]).index("Mira")]; xy3 = z["pos_xy"][list(z["pos_names"]).index("M3V")]
        fm = ((np.hypot(xx - xym[0], yy - xym[1]) <= 3) | (np.hypot(xx - xy3[0], yy - xy3[1]) <= 3)).ravel()
        q0 = ok
        Dg = D[q0].reshape(q0.sum(), -1)[:, fm]
        sig = 1.4826 * np.nanmedian(np.abs(np.diff(Dg, axis=0)), axis=0) / np.sqrt(2) + 0.5
        u = (xx - xym[0]) / 3; v = (yy - xym[1]) / 3
        for k, p in enumerate(phantoms(s)):
            Tp = prf.image(p[0] + x0, p[1] + y0)[y0:y0 + ny, x0:x0 + nx]
            A = np.column_stack([b.ravel()[fm], T3.ravel()[fm], Tp.ravel()[fm], np.ones(fm.sum()), u.ravel()[fm], v.ravel()[fm]]) / sig[:, None]
            Y = (D.reshape(len(t), -1)[:, fm] / sig).T
            cov = np.linalg.pinv(A.T @ A); c_ = cov @ A.T @ np.nan_to_num(Y)
            F = c_[2]
            raw = 1 + F / 979.0
            m = ok & np.isfinite(raw) & np.all(np.isfinite(D.reshape(len(t), -1)[:, fm]), 1)
            tr = detrend(t[m], raw[m], transit_mask(t[m]))
            det = raw[m] / tr
            key = f"phantom_{k}"
            lcs.setdefault(key, []).append((t[m], det, np.full(m.sum(), s), np.full(m.sum(), cad)))
            corr = cov[2, 1] / np.sqrt(cov[2, 2] * cov[1, 1])
            meta.setdefault(key, dict(sectors=[], pos=[], corr_with_M3V=[]))
            meta[key]["sectors"].append(s); meta[key]["pos"].append([float(p[0]), float(p[1])]); meta[key]["corr_with_M3V"].append(float(corr))
        print(f"S{s}: controls " + ", ".join(f"{c[0]} T{c[3]:.2f}" for c in ctr) + f"; phantoms {len(phantoms(s))}", flush=True)
    store = {}
    for key, parts in lcs.items():
        t = np.concatenate([p[0] for p in parts]); det = np.concatenate([p[1] for p in parts])
        sec = np.concatenate([p[2] for p in parts]); cad = np.concatenate([p[3] for p in parts])
        store[f"{key}__t"] = t; store[f"{key}__det"] = det; store[f"{key}__sector"] = sec; store[f"{key}__cad"] = cad
        meta[key]["rstd_ppt_by_sector"] = {int(s): float(robust_std(det[sec == s] - 1) * 1e3) for s in np.unique(sec)}
    np.savez_compressed(f"{OUT}/data/ctrl_lcs.npz", **store)
    save_json(meta, f"{OUT}/out/controls.json")
    for k, v in meta.items(): print(k, v.get("T", ""), v["rstd_ppt_by_sector"])


if __name__ == "__main__":
    main()
