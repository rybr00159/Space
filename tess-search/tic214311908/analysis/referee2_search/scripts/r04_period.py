"""r04: is P = 5.6915 d unique?  Cycle count and aliases with referee-2 code.

A. Transit times: phot (out/timing_table.csv, model B strict, err_cal) and referee (referee_phot/out/per_transit.csv,
   strict, tc_err_cal; E-388 excluded as in phase 1). For every trial period on a fine grid 5.55-5.85 d (step 2e-6 d)
   the cycle numbers are re-assigned by rounding relative to the best-measured event, T0 and P refitted by weighted
   least squares; chi2 of each distinct assignment. Reported: best, next-best distinct assignment, and the cycle-count
   uncertainty of S93 from S103+S104 alone and of each earlier sector from the later ones.
   Also: with errors inflated x2 and with the 3 worst outliers removed (robustness).
B. Folded depths at phases 0.5, 1/3, 2/3 (P/2, P/3 and 2P/3 aliases) on the referee-2 cleaned light curve (r01) with a
   fixed trapezoid (1.6 h, flat 0.7) and empirical errors from 400 random phases.
Output: out/r04_period.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
rng = np.random.default_rng(404)


def wls(E, tc, err):
    w = 1 / err ** 2
    A = np.column_stack([np.ones(len(E)), E])
    cov = np.linalg.inv((A * w[:, None]).T @ A)
    c = cov @ ((A * w[:, None]).T @ tc)
    return c, cov, float(np.sum(w * (tc - A @ c) ** 2))


def scan(tc, err, Pmin=5.55, Pmax=5.85, step=2e-6):
    iref = int(np.argmin(err))
    seen = {}
    for Pt in np.arange(Pmin, Pmax, step):
        E = np.round((tc - tc[iref]) / Pt).astype(int)
        key = tuple(E)
        if key in seen: continue
        c, cov, ch = wls(E, tc, err)
        seen[key] = (float(c[1]), ch)
    lst = sorted(seen.values(), key=lambda x: x[1])
    return dict(n_assign=len(lst), best_P=lst[0][0], best_chi2=lst[0][1], next=[dict(P=a, dchi2=b - lst[0][1]) for a, b in lst[1:4]])


def bootstrap_count(tc, err, sec):
    order = [[103, 104], [93], [66], [39], [12]]
    used = np.isin(sec, order[0])
    E = np.round((tc[used] - tc[used][0]) / P).astype(int)
    c, cov, _ = wls(E, tc[used], err[used])
    steps = []
    T0r = tc[used][0]
    for g in order[1:]:
        m = np.isin(sec, g)
        if not m.any(): continue
        N = (tc[m] - c[0] - c[1] * 0) / c[1]   # cycles relative to the reference event
        sigN = np.sqrt(cov[0, 0] + N ** 2 * cov[1, 1] + 2 * N * cov[0, 1] + err[m] ** 2) / c[1]
        steps.append(dict(add=g, N_median=float(np.median(N)), frac=float(np.median(N - np.round(N))), sigma_N=float(np.max(sigN)),
                          P_before=float(c[1]), P_err_before=float(np.sqrt(cov[1, 1]))))
        used |= m
        E = np.round((tc[used] - c[0]) / c[1]).astype(int)
        c0, cov, _ = wls(E, tc[used], err[used])
        c = c0
    return steps


def main():
    out = {}
    T = pd.read_csv(f"{BASE}/phot/out/timing_table.csv")
    T = T[(T.model == "B") & (T.relaxed == False) & T.tc.notna()]
    R = pd.read_csv(f"{BASE}/referee_phot/out/per_transit.csv")
    R = R[(R.sel == "strict") & R.tc.notna() & (R.epoch != -388)]
    sets = dict(phot=(T.tc.values, T.err_cal.values, T.sector.values.astype(int)),
                referee=(R.tc.values, R.tc_err_cal.values, R.sector.values.astype(int)))
    for nm, (tc, err, sec) in sets.items():
        o = {}
        o["scan"] = scan(tc, err)
        o["scan_err_x2"] = scan(tc, 2 * err)
        # drop the 3 largest |O-C|/err
        E = np.round((tc - tc[np.argmin(err)]) / P).astype(int)
        c, cov, ch = wls(E, tc, err)
        pull = np.abs(tc - (c[0] + c[1] * E)) / err
        keep = np.argsort(pull)[:-3]
        o["scan_drop3"] = scan(tc[keep], err[keep])
        o["linear"] = dict(P=float(c[1]), P_err=float(np.sqrt(cov[1, 1])), chi2=ch, n=len(tc))
        o["cycle_bootstrap"] = bootstrap_count(tc, err, sec)
        out[nm] = o
        print(nm, json.dumps(o), flush=True)
    # ---- B alias phases on the cleaned light curve
    z = np.load(f"{ME}/data/m3v_clean_r2.npz")
    tb, yb, nb, sb = binlc(z["t"], z["clean"] - 1, z["s"], z["cad"], 10.0)
    w = np.zeros(len(yb))
    for q in np.unique(sb):
        m = sb == q; yb[m] -= np.median(yb[m]); w[m] = 1 / rstd(yb[m]) ** 2
    cadb = np.array([max(600.0, np.median(z["cad"][z["s"] == q])) for q in sb])

    def dfit(T0_):
        tm = trapz_periodic(tb, P, T0_, cadb)
        a = np.sum(w * tm * tm)
        return float(-np.sum(w * tm * yb) / a), float(1 / np.sqrt(a))
    fk = []
    for ph in rng.uniform(0.04, 0.96, 500):
        if min(abs(ph - 0.5), abs(ph - 1 / 3), abs(ph - 2 / 3)) < 0.03: continue
        fk.append(dfit(T0 + ph * P)[0])
        if len(fk) >= 400: break
    sd = float(np.std(fk))
    B = dict(err_emp_ppt=sd * 1e3, primary=dfit(T0)[0] * 1e3)
    for nm, ph in (("0.5", 0.5), ("1/3", 1 / 3), ("2/3", 2 / 3)):
        B[nm] = dfit(T0 + ph * P)[0] * 1e3
    B["P2_alias_sigma"] = (B["primary"] - B["0.5"]) / (np.sqrt(2) * B["err_emp_ppt"])
    out["alias_phases"] = B
    print(B)
    jdump(out, f"{ME}/out/r04_period.json")


if __name__ == "__main__":
    main()
