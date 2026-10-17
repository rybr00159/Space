"""s02: empirical false-alarm calibration of the BLS searches (same data, same grid, same statistic as s01).

Null realisations of the deblended M3V model-B light curve (5.69-d transits masked, 10-min bins):
  inv_roll : flux inverted (dips <-> bumps) and cyclically rolled within each sector by a random number of bins
             (time stamps, gaps and weights kept) -> same noise and systematics statistics, no real transits.
  roll     : not inverted, rolled within each sector -> real dips kept but inter-sector phase coherence destroyed.
Plus single realisations of the deblended pair (model A) and Mira (model B) light curves (same masking).
For every realisation: global max Delta-chi2 of the coherent search (0.3-30 d), max in period bands, and the
max of each per-sector search. Usage: s02_bls_null.py <kind> <n> <seed>
Outputs: out/null_<kind>_<seed>.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s01_bls_real import search_data

BANDS = [(0.3, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 30.0)]
GROUPS = {"S12": [12], "S39": [39], "S66": [66], "S93": [93], "S103": [103], "S104": [104], "S103_104": [103, 104]}


def star_data(star, model):
    d = load_star(star, model)
    m = ~d.intr.values & np.isfinite(d.det.values)
    tb, fb, nb_, sb, cb = bin_lc(d.t.values[m], d.det.values[m], d.sector.values[m], d.cad.values[m], 10.0)
    y, w = prep_weights(fb - 1, sb)
    return tb, y, w, sb


def summarize(t, y, w, sb, fr_all, fr_sec):
    r = run_bls(t, y, w, fr_all)
    P = 1 / r["freq"]; d2 = r["dchi2"]
    i = int(np.argmax(d2))
    res = dict(max=float(d2[i]), P=float(P[i]), depth_ppt=float(r["depth"][i] * 1e3), sde=float((d2[i] - d2.mean()) / d2.std()),
               bands={f"{a}-{b}": float(d2[(P >= a) & (P < b)].max()) for a, b in BANDS})
    sec = {}
    for g, ss in GROUPS.items():
        m = np.isin(sb, ss)
        yy = y[m] - np.sum(w[m] * y[m]) / w[m].sum()
        rs = run_bls(t[m], yy, w[m], fr_sec[g])
        j = int(np.argmax(rs["dchi2"]))
        sec[g] = dict(max=float(rs["dchi2"][j]), P=float(1 / rs["freq"][j]), sde=float((rs["dchi2"][j] - rs["dchi2"].mean()) / rs["dchi2"].std()))
    res["sectors"] = sec
    return res


def roll_within(y, w, sb, rng):
    y2 = y.copy(); w2 = w.copy()
    for s in np.unique(sb):
        m = np.where(sb == s)[0]
        k = rng.integers(len(m) // 10, len(m) - len(m) // 10)
        y2[m] = np.roll(y[m], k); w2[m] = np.roll(w[m], k)
    return y2, w2


def main():
    kind = sys.argv[1]; n = int(sys.argv[2]); seed = int(sys.argv[3])
    flux = sys.argv[4] if len(sys.argv) > 4 else "clean"
    rng = np.random.default_rng(seed)
    t, y, w, sb, cb = search_data("B", flux)
    B = t.max() - t.min()
    fr_all = freq_grid(0.3, 30.0, B)
    fr_sec = {}
    for g, ss in GROUPS.items():
        m = np.isin(sb, ss); span = t[m].max() - t[m].min()
        fr_sec[g] = freq_grid(0.3, min(30.0, span / 2), span)
    out = []
    if kind == "ctrl":
        # control stars and phantom positions (s03), each with >= 5 sectors, RRab leak removed per sector, and inverted
        from s07_periodic_model import harm_X, harm_fit
        Z = np.load(f"{ROOT}/data/ctrl_lcs.npz")
        names = sorted({k.split("__")[0] for k in Z.files})
        for nm in names:
            tt, dd, ss_, cc_ = Z[f"{nm}__t"], Z[f"{nm}__det"], Z[f"{nm}__sector"], Z[f"{nm}__cad"]
            if len(np.unique(ss_)) < 5: continue
            yy = dd - 1
            intr = transit_mask(tt)
            for sec in np.unique(ss_):
                m = (ss_ == sec) & ~intr & np.isfinite(yy)
                c, _, _ = harm_fit(tt[m], yy[m], np.ones(m.sum()), 0.648528)
                mf = ss_ == sec
                yy[mf] = yy[mf] - (harm_X(tt[mf], 0.648528) @ c - c[0])
            k = ~intr & np.isfinite(yy)
            tb, fb, nb_, sb_, cb_ = bin_lc(tt[k], 1 + yy[k], ss_[k], cc_[k], 10.0)
            y0, w0 = prep_weights(fb - 1, sb_)
            fr_sec2 = {}
            for g, ss in GROUPS.items():
                m = np.isin(sb_, ss)
                if m.sum() < 100: continue
                span = tb[m].max() - tb[m].min(); fr_sec2[g] = freq_grid(0.3, min(30.0, span / 2), span)
            for sign in (1, -1):
                GR = dict(GROUPS)
                res = summarize(tb, sign * y0, w0, sb_, fr_all, {g: fr_sec2[g] for g in fr_sec2}) if False else None
                r = run_bls(tb, sign * y0, w0, fr_all)
                P_ = 1 / r["freq"]; d2 = r["dchi2"]; i = int(np.argmax(d2))
                res = dict(name=nm, sign=sign, max=float(d2[i]), P=float(P_[i]), depth_ppt=float(r["depth"][i] * 1e3),
                           sde=float((d2[i] - d2.mean()) / d2.std()),
                           bands={f"{a}-{b}": float(d2[(P_ >= a) & (P_ < b)].max()) for a, b in BANDS})
                sec = {}
                for g in fr_sec2:
                    m = np.isin(sb_, GROUPS[g])
                    yy2 = sign * y0[m] - np.sum(w0[m] * sign * y0[m]) / w0[m].sum()
                    rs = run_bls(tb[m], yy2, w0[m], fr_sec2[g]); j = int(np.argmax(rs["dchi2"]))
                    sec[g] = dict(max=float(rs["dchi2"][j]), P=float(1 / rs["freq"][j]))
                res["sectors"] = sec
                out.append(res)
                print(nm, sign, f"{res['max']:.1f} P {res['P']:.3f}", flush=True)
                save_json(dict(kind=kind, seed=seed, flux=flux, n=len(out), res=out, partial=True), f"{ROOT}/out/null_{kind}_{flux}_{seed}.json")
    elif kind in ("pair_A", "Mira_B", "pair_A_inv", "Mira_B_inv"):
        star, model = kind.split("_")[:2]
        tt, yy, ww, ss = star_data(star, model)
        if kind.endswith("inv"): yy = -yy
        fr_all_s = freq_grid(0.3, 30.0, B)
        out.append(summarize(tt, yy, ww, ss, fr_all_s, fr_sec))
        print(kind, out[-1]["max"], out[-1]["P"], flush=True)
    else:
        for k in range(n):
            y2, w2 = roll_within(y, w, sb, rng)
            if kind == "inv_roll": y2 = -y2
            out.append(summarize(t, y2, w2, sb, fr_all, fr_sec))
            print(k, f"{out[-1]['max']:.1f} P {out[-1]['P']:.4f}", flush=True)
            save_json(dict(kind=kind, seed=seed, flux=flux, n=len(out), res=out, partial=True), f"{ROOT}/out/null_{kind}_{flux}_{seed}.json")
    save_json(dict(kind=kind, seed=seed, flux=flux, n=len(out), res=out), f"{ROOT}/out/null_{kind}_{flux}_{seed}.json")


if __name__ == "__main__":
    main()
