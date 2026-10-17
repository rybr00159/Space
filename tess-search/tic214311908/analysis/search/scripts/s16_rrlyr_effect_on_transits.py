"""s16: effect of the RRab leak (0.648528 d, 1.8-3.9 ppt first-harmonic semi-amplitude in the model-B M3V curve) on the
5.69-d transit measurements: global and per-transit depths (fixed trapezoid, T14 1.6 h, flat 70 %, referee ephemeris,
local linear baseline within +-0.3 d) on the phot detrended flux (det) and on the RRab-cleaned flux (det_clean), and
the RRab model value at each transit. Empirical per-transit errors from 100 fake epochs per sector.
Output: out/rrlyr_transit_effect.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s12_alias import trap
RNG = np.random.default_rng(16)


def local_depth(t, y, cad, tc, half=0.3):
    k = np.abs(t - tc) < half
    if k.sum() < 10: return np.nan
    tm = trap(t[k], P_B, tc, cad[k])
    if (tm > 0.5).sum() < 2: return np.nan
    X = np.column_stack([np.ones(k.sum()), t[k] - tc, -tm])
    c, *_ = np.linalg.lstsq(X, y[k], rcond=None)
    return float(c[2])


def main():
    z = np.load(f"{ROOT}/data/m3v_clean.npz")
    t, det, cl, s, cad = z["t"], z["det"], z["det_clean"], z["sector"], z["cad"].astype(float)
    rows = []
    for e in range(-400, 70):
        tc = T0_B + e * P_B
        if not np.any(np.abs(t - tc) < 0.05): continue
        d0 = local_depth(t, det - 1, cad, tc); d1 = local_depth(t, cl - 1, cad, tc)
        k = np.abs(t - tc) < 0.04
        rows.append(dict(epoch=e, sector=int(s[np.argmin(np.abs(t - tc))]), depth_det_ppt=d0 * 1e3, depth_clean_ppt=d1 * 1e3,
                         rr_model_in_transit_ppt=float(np.mean((det - cl)[k]) * 1e3)))
    D = pd.DataFrame(rows).dropna()
    # empirical single-transit error per sector (fake epochs on the cleaned flux)
    err = {}
    for sec in SECTORS:
        m = s == sec; tt = t[m]
        fk = []
        for _ in range(100):
            tc = RNG.uniform(tt.min() + 0.3, tt.max() - 0.3)
            if np.min(np.abs(((tc - T0_B) / P_B + 0.5) % 1 - 0.5) * P_B) < 0.3: continue
            fk.append(local_depth(t, cl - 1, cad, tc))
        err[sec] = float(np.nanstd(fk) * 1e3)
    D["err_emp_ppt"] = D.sector.map(err)
    def wstats(col):
        w = 1 / D.err_emp_ppt ** 2
        mu = float(np.sum(w * D[col]) / np.sum(w)); chi2 = float(np.sum(w * (D[col] - mu) ** 2))
        return mu, float(1 / np.sqrt(np.sum(w))), chi2
    odd = D.epoch % 2 == 1
    def oe(col):
        w = 1 / D.err_emp_ppt ** 2
        mo = np.sum(w[odd] * D[col][odd]) / np.sum(w[odd]); me = np.sum(w[~odd] * D[col][~odd]) / np.sum(w[~odd])
        return float(mo - me), float(np.sqrt(1 / np.sum(w[odd]) + 1 / np.sum(w[~odd])))
    out = dict(n=len(D), per_sector_err_ppt=err,
               det=dict(zip(("mean_ppt", "err_ppt", "chi2"), wstats("depth_det_ppt")), odd_minus_even=oe("depth_det_ppt")),
               clean=dict(zip(("mean_ppt", "err_ppt", "chi2"), wstats("depth_clean_ppt")), odd_minus_even=oe("depth_clean_ppt")),
               rms_change_ppt=float(np.std(D.depth_clean_ppt - D.depth_det_ppt)), max_abs_change_ppt=float(np.max(np.abs(D.depth_clean_ppt - D.depth_det_ppt))),
               table=D.to_dict(orient="records"))
    save_json(out, f"{ROOT}/out/rrlyr_transit_effect.json")
    print(json.dumps({k: v for k, v in out.items() if k != "table"}, indent=0))
    print(D.round(2).to_string())


if __name__ == "__main__":
    main()
