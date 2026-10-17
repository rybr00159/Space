"""s09: injection-recovery completeness of the coherent BLS search (s01 'clean' mode) for additional transiting planets.

Each injection: a limb-darkened (u = 0.25, 0.35), exposure-integrated batman transit (circular orbit, rho = 9.4 g/cm3,
impact parameter b ~ U(0, 0.9), random phase) is multiplied into the undetrended M3V model-B flux (phot
flux_raw_norm, rescaled to 979 e-/s), then the full chain is re-run: wotan biweight 0.75 d (5.69-d windows masked) ->
per-sector 3-harmonic fit of the RRab leak at 0.648528 d (refit on the injected data) and subtraction -> 5.69-d windows
removed -> 10-min bins -> per-sector white weights -> BLS.
Search: the frequencies of the full s01 grid within +-1.5 % of the injected frequency (the global threshold, from the
s02 nulls, already accounts for the look-elsewhere effect of the full grid).
Recovered (strict): Delta-chi2 above the threshold at a period within 0.1 % of P_inj and the folded mid-time within
0.6 T14 of the injected one. Also stored: the best Delta-chi2 and period in the window, to allow other thresholds.
Grid: P log-uniform 0.5-30 d, Rp log-uniform 0.8-6 R_earth (R* = 0.381 R_sun).
Usage: s09_inject.py <n> <seed>      Output: out/inject_<seed>.csv
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s07_periodic_model import harm_X, harm_fit

P_RR = 0.648528


def prepare():
    d = load_star("M3V", "B")
    return d


def chain(d, raw_inj):
    """raw -> detrend -> RRab-leak removal -> masked, binned, weighted arrays."""
    t = d.t.values; s = d.sector.values; intr = d.intr.values
    det = np.empty(len(t))
    for sec in SECTORS:
        m = s == sec
        tr = detrend(t[m], raw_inj[m], intr[m])
        det[m] = raw_inj[m] / tr
    y = det - 1
    for sec in SECTORS:
        m = (s == sec) & ~intr & np.isfinite(y)
        w = np.full(m.sum(), 1.0)
        c, _, _ = harm_fit(t[m], y[m], w, P_RR)
        mf = s == sec
        y[mf] = y[mf] - (harm_X(t[mf], P_RR) @ c - c[0])
    k = ~intr & np.isfinite(y)
    tb, fb, nb_, sb, cb = bin_lc(t[k], 1 + y[k], s[k], d.cad.values[k], 10.0)
    yy, ww = prep_weights(fb - 1, sb)
    return tb, yy, ww, sb


def main():
    n = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    d = prepare()
    t = d.t.values; raw = d.raw.values; cad = d.cad.values.astype(float)
    tb0, yy0, ww0, sb0 = chain(d, raw)
    B = tb0.max() - tb0.min()
    FR = freq_grid(0.3, 30.0, B)
    rows = []
    for i in range(n):
        P = float(np.exp(rng.uniform(np.log(0.5), np.log(30.0))))
        rp = float(np.exp(rng.uniform(np.log(0.8), np.log(6.0))))
        b = float(rng.uniform(0, 0.9))
        t0 = float(T0_B + rng.uniform(0, 1) * P)
        k = rp * REARTH_RSUN / RSTAR
        mdl = batman_model(t, t0, P, k, b, cad)
        tb, yy, ww, sb = chain(d, raw * mdl)
        f0 = 1 / P
        sel = (FR > f0 * (1 - 0.015)) & (FR < f0 * (1 + 0.015))
        r = run_bls(tb, yy, ww, FR[sel])
        j = int(np.argmax(r["dchi2"]))
        Pb = 1 / r["freq"][j]; T0b = r["t0"][j]
        # strict match
        good = np.abs(1 / r["freq"] - P) / P < 1e-3
        tmid = 0.5 * (tb.min() + tb.max())
        ep = np.round((tmid - t0) / P); tref = t0 + ep * P
        dt = np.abs(((r["t0"] - tref) / (1 / r["freq"]) + 0.5) % 1 - 0.5) * (1 / r["freq"])
        T14 = t14_central(P) * np.sqrt(max(1 - b ** 2, 0.05))
        match = good & (dt < 0.6 * max(T14, 0.03))
        d2m = float(r["dchi2"][match].max()) if match.any() else 0.0
        depth_inj = float(1 - mdl.min())
        # how many transits fall in the (masked, quality) data
        ph = ((tb - t0) / P + 0.5) % 1 - 0.5
        ntr = len(np.unique(np.round((tb[np.abs(ph * P) < T14 / 2] - t0) / P)))
        rows.append(dict(P=P, rp=rp, b=b, t0=t0, depth_ppt=depth_inj * 1e3, T14_h=T14 * 24, n_transits=ntr,
                         dchi2_best=float(r["dchi2"][j]), P_best=float(Pb), dchi2_match=d2m,
                         best_is_match=bool(match[j])))
        if i % 20 == 0: print(i, rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{ROOT}/out/inject_{seed}.csv", index=False)


if __name__ == "__main__":
    main()
