"""r05: injection-recovery spot check of the additional-planet search with referee-2 code.

Each injection: batman transit (quadratic LD u = 0.25, 0.35; circular; rho = 9.4 g/cm3; b ~ U(0, 0.85); random
epoch; exposure-integrated) multiplied into the undetrended M3V model-B flux (phot flux_raw_norm, rescaled), then:
wotan biweight 0.75 d per sector with the 5.69-d windows (+-0.1 d) masked -> per-sector 3-harmonic removal at the
Gaia RRab period (refit) -> 5.69-d windows removed -> 10-min bins -> per-sector weights -> bls2 on the r03 grid within
+-1.5 % of the injected frequency. Stored: best dchi2 at a matched period (|dP|/P < 1e-3, epoch within 0.6 T14).
Grid: P = 0.7, 1.5, 3, 6, 12, 22 d x Rp = 2.0, 2.5, 3.0, 3.5, 4.5, 5.5 R_earth, n per cell.
Recovery fractions are evaluated in r06 against the r03 null thresholds (global and per period band).
Usage: r05_inject.py <n_per_cell> <seed>  -> out/r05_inject_<seed>.csv
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
import bls2, batman
from wotan import flatten


def bat(t, cad, t0, Pp, k, b):
    pm = batman.TransitParams()
    pm.t0 = t0; pm.per = Pp; pm.rp = k; pm.a = a_R(Pp); pm.inc = np.degrees(np.arccos(b / pm.a))
    pm.ecc = 0.; pm.w = 90.; pm.u = [0.25, 0.35]; pm.limb_dark = "quadratic"
    out = np.ones(len(t))
    for c in np.unique(cad):
        m = cad == c
        out[m] = batman.TransitModel(pm, t[m], supersample_factor=max(1, int(c // 120)), exp_time=c / 86400).light_curve(pm)
    return out


def chain(t, s, cad, raw, intr):
    y = np.empty(len(t))
    for q in SECTORS:
        m = s == q
        mm = m & ~intr
        tr = flatten(t[mm], raw[mm], method="biweight", window_length=0.75, break_tolerance=0.5, edge_cutoff=0, return_trend=True)[1]
        ok = np.isfinite(tr)
        y[m] = raw[m] / np.interp(t[m], t[mm][ok], tr[ok]) - 1
    for q in SECTORS:
        m = (s == q) & ~intr
        X = harmX(t[m], P_RR_GAIA)
        c, *_ = np.linalg.lstsq(X, y[m], rcond=None)
        mf = s == q
        y[mf] -= harmX(t[mf], P_RR_GAIA) @ c - c[0]
    k = ~intr
    tb, yb, nb, sb = binlc(t[k], y[k], s[k], cad[k], 10.0)
    yy, ww = bls2.weights(yb, sb)
    return tb, yy, ww, sb


def main():
    n = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    d = load_m3v("B")
    t, s, cad, raw, intr = d.t.values, d.s.values, d.cad.values, d.raw.values, d.intr.values
    tb0, _, _, _ = chain(t, s, cad, raw, intr)
    FR = bls2.freq_grid(0.3, 30.0, tb0.max() - tb0.min())
    rows = []
    for Pp in (0.7, 1.5, 3.0, 6.0, 12.0, 22.0):
        for rp in (2.0, 2.5, 3.0, 3.5, 4.5, 5.5):
            for j in range(n):
                P_ = Pp * np.exp(rng.uniform(-0.05, 0.05))
                if abs(P_ / P - round(P_ / P)) < 0.02 or abs(P / P_ - round(P / P_)) < 0.02:
                    P_ *= 1.03                                    # avoid commensurability with the 5.69-d mask
                b = rng.uniform(0, 0.85); t0 = T0 + rng.uniform(0, 1) * P_
                k = rp * RE_RSUN / RSTAR
                mdl = bat(t, cad, t0, P_, k, b)
                tb, yy, ww, sb = chain(t, s, cad, raw * mdl, intr)
                f0 = 1 / P_
                sel = (FR > f0 * 0.985) & (FR < f0 * 1.015)
                r = bls2.run(tb, yy, ww, FR[sel])
                Pr = 1 / r["freq"]
                T14 = t14c(P_) * np.sqrt(max(1 - b ** 2, 0.05))
                tmid = 0.5 * (tb.min() + tb.max()); tref = t0 + np.round((tmid - t0) / P_) * P_
                dt = np.abs(((r["t0"] - tref) / Pr + 0.5) % 1 - 0.5) * Pr
                match = (np.abs(Pr - P_) / P_ < 1e-3) & (dt < 0.6 * max(T14, 0.03))
                ph = phase_dist(tb, P_, t0)
                ntr = len(np.unique(np.round((tb[np.abs(ph) < T14 / 2] - t0) / P_)))
                rows.append(dict(Pcell=Pp, Rcell=rp, P=P_, rp=rp, b=b, depth_ppt=(1 - mdl.min()) * 1e3, T14_h=T14 * 24, n_tr=ntr,
                                 d2_match=float(r["dchi2"][match].max()) if match.any() else 0.0, d2_best=float(r["dchi2"].max())))
            print(Pp, rp, np.round([x["d2_match"] for x in rows[-n:]], 1), flush=True)
            pd.DataFrame(rows).to_csv(f"{ME}/out/r05_inject_{seed}.csv", index=False)


if __name__ == "__main__":
    main()
