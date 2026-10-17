"""s12: is the 5.69-d period unique? Aliases (P/3, P/2, 2P/3, 3P/2, 2P, 3P) and the cycle count across the
2019-2026 gaps.

A. Light curve: coherent BLS (s01 machinery) of the cleaned M3V light curve WITHOUT the 5.69-d mask over the full grid
   (0.3-30 d); Delta-chi2 at the simple aliases relative to the peak; fine scan of the alias comb around P
   (teeth spaced by P^2 / baseline) with the trapezoid shape fixed.
B. Transit times (phot out/timing_table.csv, model B strict, injection-calibrated errors err_cal; and the referee's
   referee_phot/out/per_transit.csv, tc_err_cal): (i) per-sector-group local periods; (ii) progressive ephemeris
   bootstrapping from S103+S104 back to S12 with the cycle-count uncertainty sigma_N at each gap;
   (iii) timing chi2 as a function of trial period over +-5 alias teeth (T0 free, cycle numbers re-assigned).
C. Depth at phase 0.5 / 1/3 / 2/3 (for P/2, P/3 aliases) on the cleaned light curve with the same trapezoid shape and
   empirical errors (random phases), and odd/even for 2P.
Output: out/alias.json, figs/alias.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RNG = np.random.default_rng(12)
T14, FLAT = 1.6 / 24, 0.7


def trap(t, P, T0, cad_s):
    out = np.zeros(len(t))
    for c in np.unique(cad_s):
        m = cad_s == c
        ss = np.linspace(-0.5, 0.5, 7) * c / 86400
        acc = np.zeros(m.sum())
        for s_ in ss:
            x = np.abs((((t[m] + s_ - T0) / P + 0.5) % 1 - 0.5) * P)
            ing = (T14 - FLAT * T14) / 2
            acc += np.where(x <= FLAT * T14 / 2, 1.0, np.where(x < T14 / 2, (T14 / 2 - x) / ing, 0.0))
        out[m] = acc / len(ss)
    return out


def depth_fit(t, y, w, tm):
    a = np.sum(w * tm * tm)
    if a <= 0: return np.nan, np.nan
    return float(-np.sum(w * tm * y) / a), float(1 / np.sqrt(a))


def lin_fit(E, tc, err):
    w = 1 / err ** 2
    A = np.column_stack([np.ones(len(E)), E]) * np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(A, tc * np.sqrt(w), rcond=None)
    cov = np.linalg.inv(A.T @ A)
    chi2 = float(np.sum((A @ c - tc * np.sqrt(w)) ** 2))
    return c, cov, chi2


def main():
    out = {}
    z = np.load(f"{ROOT}/data/m3v_clean.npz")
    t, f, s, cad = z["t"], z["det_clean"], z["sector"], z["cad"].astype(float)
    tb, fb, nb_, sb, cb = bin_lc(t, f, s, cad, 10.0)
    y, w = prep_weights(fb - 1, sb)
    # ---- A: full-grid BLS without the mask
    fr = freq_grid(0.3, 30.0, tb.max() - tb.min())
    r = run_bls(tb, y, w, fr)
    i = int(np.argmax(r["dchi2"])); Pbest = 1 / fr[i]
    P = P_B
    al = {}
    for nm, fac in (("P/3", 1 / 3), ("P/2", 0.5), ("2P/3", 2 / 3), ("P", 1.0), ("3P/2", 1.5), ("2P", 2.0), ("3P", 3.0)):
        sel = np.abs(fr - 1 / (P * fac)) < 0.002 / (P * fac)
        j = np.argmax(r["dchi2"][sel])
        al[nm] = dict(P=float(1 / fr[sel][j]), dchi2=float(r["dchi2"][sel][j]))
    # second-highest distinct peak (outside +-0.2 % of P and its simple aliases)
    mask = np.ones(len(fr), bool)
    for fac in (1 / 3, 0.5, 2 / 3, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0):
        mask &= np.abs(fr - 1 / (P * fac)) > 0.004 / (P * fac)
    j2 = int(np.argmax(np.where(mask, r["dchi2"], -1)))
    out["bls_unmasked"] = dict(P_best=float(Pbest), dchi2_best=float(r["dchi2"][i]), depth_ppt=float(r["depth"][i] * 1e3),
                               aliases=al, next_unrelated=dict(P=float(1 / fr[j2]), dchi2=float(r["dchi2"][j2])))
    print("A:", json.dumps(out["bls_unmasked"]), flush=True)
    # fine comb scan with the trapezoid shape (T0 refit by a small grid per period)
    B = tb.max() - tb.min()
    dP = P ** 2 / B
    Ps = P + np.arange(-6, 6.001, 0.02) * dP
    chi_comb = []
    for Pt in Ps:
        best = -np.inf
        for dt in np.linspace(-0.03, 0.03, 31):
            tm = trap(tb, Pt, T0_B + dt, cb)
            d_, e_ = depth_fit(tb, y, w, tm)
            best = max(best, (d_ / e_) ** 2 if d_ > 0 else 0)
        chi_comb.append(best)
    chi_comb = np.array(chi_comb)
    k0 = int(np.argmax(chi_comb))
    mins = [k for k in range(1, len(Ps) - 1) if chi_comb[k] > chi_comb[k - 1] and chi_comb[k] >= chi_comb[k + 1]]
    teeth = sorted([(float(Ps[k]), float(chi_comb[k0] - chi_comb[k])) for k in mins], key=lambda x: x[1])[:6]
    out["comb"] = dict(tooth_spacing_d=float(dP), best_P=float(Ps[k0]), teeth=teeth)
    print("comb:", out["comb"], flush=True)
    # ---- B: transit times
    T = pd.read_csv(f"{PHOT}/out/timing_table.csv")
    T = T[(T.model == "B") & (T.relaxed == False) & T.tc.notna()].copy()
    R = pd.read_csv(f"{BASE}/referee_phot/out/per_transit.csv"); R = R[(R.sel == "strict") & R.tc.notna() & (R.epoch != -388)]
    sets = {"phot": (T.epoch.values.astype(int), T.tc.values, T.err_cal.values, T.sector.values.astype(int)),
            "referee": (R.epoch.values.astype(int), R.tc.values, R.tc_err_cal.values, R.sector.values.astype(int))}
    out["timing"] = {}
    for nm, (E, tc, err, sec) in sets.items():
        res = {}
        c, cov, chi2 = lin_fit(E, tc, err)
        res["global"] = dict(T0=float(c[0]), P=float(c[1]), P_err=float(np.sqrt(cov[1, 1])), chi2=chi2, n=len(E))
        groups = [("S103_104", [103, 104]), ("S93", [93]), ("S66", [66]), ("S39", [39]), ("S12", [12])]
        loc = {}
        for g, ss in groups:
            m = np.isin(sec, ss)
            if m.sum() >= 2:
                cc, cv, ch = lin_fit(E[m] - E[m][0], tc[m], err[m])
                loc[g] = dict(P=float(cc[1]), P_err=float(np.sqrt(cv[1, 1])), n=int(m.sum()), span_cycles=int(E[m].max() - E[m].min()))
        res["local_periods"] = loc
        # progressive bootstrapping: use only time differences, re-derive cycle numbers from scratch
        steps = []
        used = np.isin(sec, [103, 104])
        tref = tc[used][0]
        cc, cv, _ = lin_fit(np.round((tc[used] - tref) / P).astype(int), tc[used], err[used])
        for g, ss in groups[1:]:
            m = np.isin(sec, ss)
            if not m.any(): continue
            Nf = (tc[m] - cc[0]) / cc[1]
            sigN = np.sqrt(cv[0, 0] + (Nf ** 2) * cv[1, 1] + 2 * Nf * cv[0, 1] + err[m] ** 2) / cc[1]
            frac = Nf - np.round(Nf)
            steps.append(dict(add=g, N=float(np.median(Nf)), N_frac_offset=float(np.median(frac)), sigma_N=float(np.max(sigN)),
                              P_before=float(cc[1]), P_err_before=float(np.sqrt(cv[1, 1]))))
            used = used | m
            cc, cv, _ = lin_fit(np.round((tc[used] - cc[0]) / cc[1]).astype(int), tc[used], err[used])
        res["bootstrap"] = steps
        # timing chi2 vs trial period: each trial period defines a cycle assignment; distinct assignments are compared
        B_ = tc.max() - tc.min()
        dPt = P ** 2 / B_
        Pt = P + np.arange(-5, 5.0001, 0.005) * dPt
        seen = {}
        for p in Pt:
            En = np.round((tc - tc[np.argmin(err)]) / p).astype(int)
            key = tuple(En - En.min())
            if key in seen: continue
            cc_, cv_, ch_ = lin_fit(En, tc, err)
            seen[key] = (float(cc_[1]), ch_)
        lst = sorted(seen.values(), key=lambda x: x[1])
        res["chi2_scan"] = dict(tooth_spacing_d=float(dPt), n_assignments=len(lst), best=dict(P=lst[0][0], chi2=lst[0][1]),
                                next=[dict(P=a, dchi2=b - lst[0][1]) for a, b in lst[1:5]])
        out["timing"][nm] = res
        print(nm, json.dumps(res), flush=True)
    # ---- C: depths at alias phases (cleaned LC, all data, empirical errors)
    tm0 = trap(tb, P, T0_B, cb)
    d0, e0 = depth_fit(tb, y, w, tm0)
    fake = []
    for _ in range(300):
        ph = RNG.uniform(0.05, 0.95)
        if min(abs(ph - 0.5), abs(ph - 1 / 3), abs(ph - 2 / 3)) < 0.03: continue
        fake.append(depth_fit(tb, y, w, trap(tb, P, T0_B + ph * P, cb))[0])
    sd = float(np.std(fake))
    ph_res = {}
    for nm, ph in (("0.5", 0.5), ("1/3", 1 / 3), ("2/3", 2 / 3)):
        d_, e_ = depth_fit(tb, y, w, trap(tb, P, T0_B + ph * P, cb))
        ph_res[nm] = dict(depth_ppt=d_ * 1e3, err_emp_ppt=sd * 1e3, err_white_ppt=e_ * 1e3)
    # odd / even with empirical errors from random-phase fakes of the odd and even subsets
    ep = np.round((tb - T0_B) / P).astype(int)
    odd = (ep % 2 == 1)
    dodd = depth_fit(tb[odd], y[odd], w[odd], tm0[odd]); deven = depth_fit(tb[~odd], y[~odd], w[~odd], tm0[~odd])
    fo, fe = [], []
    for _ in range(300):
        ph = RNG.uniform(0.05, 0.95)
        if abs(ph - 0.5) < 0.03: continue
        tmf = trap(tb, P, T0_B + ph * P, cb)
        epf = np.round((tb - T0_B - ph * P) / P).astype(int); of = epf % 2 == 1
        fo.append(depth_fit(tb[of], y[of], w[of], tmf[of])[0]); fe.append(depth_fit(tb[~of], y[~of], w[~of], tmf[~of])[0])
    sdiff = float(np.std(np.array(fo) - np.array(fe)))
    out["phases"] = dict(primary_depth_ppt=d0 * 1e3, primary_err_emp_ppt=sd * 1e3, alias_phases=ph_res,
                         odd_ppt=dodd[0] * 1e3, even_ppt=deven[0] * 1e3, odd_minus_even_err_emp_ppt=sdiff * 1e3,
                         odd_minus_even_sigma_emp=float((dodd[0] - deven[0]) / sdiff),
                         note="model-B light curve with the RRab leak removed; the referee showed the model-B odd/even excess is a photometry artefact")
    print("C:", json.dumps(out["phases"]), flush=True)
    save_json(out, f"{ROOT}/out/alias.json")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    Pg = 1 / fr
    ax[0].plot(Pg, r["dchi2"], lw=0.3); ax[0].set_xscale("log"); ax[0].set_xlabel("P (d)"); ax[0].set_ylabel("BLS dchi2 (unmasked, cleaned)")
    for fac in (1 / 3, 0.5, 2 / 3, 1, 1.5, 2, 3): ax[0].axvline(P * fac, color="r", lw=0.5, ls=":")
    ax[1].plot((Ps - P) / dP, chi_comb, "-"); ax[1].set_xlabel("(P - 5.6915087) / (P^2/baseline)"); ax[1].set_ylabel("(depth/err)^2, trapezoid")
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/alias.png", dpi=80); plt.close(fig)


if __name__ == "__main__":
    main()
