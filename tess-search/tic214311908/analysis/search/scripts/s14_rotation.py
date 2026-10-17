"""s14: rotational modulation of the M3V in the deblended photometry, with pixel-level injections for sensitivity.

For each sector the phot photometry (p05, read-only import) is re-run on the TESScut cube (a) as is and (b) with a
sinusoid of semi-amplitude A (e-/s) and period P injected into the cube on the M3V PRF (p05 'inject' option), on a grid
P = 1, 2, 4, 6, 8.7, 12 d and A = 5, 10, 20, 40 e-/s (0.5-4 % of the M3V's 979 e-/s), random phase.
Light curve for rotation (as phot p12): the undetrended model-B (and model-A) M3V flux, quality_ok cadences outside the
5.69-d windows, regressed on [1, background excess, dx, dy] (option: + the Mira flux, as in p12), with the RRab leak
removed (3 harmonics at 0.648528 d), binned to 1 h.
Measured: Lomb-Scargle periodogram 0.5-14 d; injected signal 'detected' if the highest peak lies within 10 % of P_inj
and its power exceeds the 99th percentile of the highest peaks of the un-injected light curve's block-shuffled
(1-d blocks) realisations; recovered amplitude = sinusoid + linear-trend fit at P_inj, minus the same fit on the
un-injected data (vector difference).
The un-injected periodograms of M3V-B, M3V-A, pair-A and Mira-B are reported per sector (top peaks).
Output: out/rotation_inject.json, out/rotation_real.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/search/scripts")
sys.path.insert(1, "tic214311908/analysis/phot/scripts")
from slib import *
from s07_periodic_model import harm_X, harm_fit
from astropy.timeseries import LombScargle

PGRID = (1.0, 2.0, 4.0, 6.0, 8.7, 12.0)
AGRID = (5.0, 10.0, 20.0, 40.0)
RNG = np.random.default_rng(14)


def lc_from_R(R, L, key_model="B", star="M3V", mira_reg=False):
    names = R["FIT"][key_model]["names"]
    F = R["FIT"][key_model]["F"][:, names.index(star)]
    t = L["t_tdb"]; ok = L["quality_ok"].astype(bool) & ~L["intransit_window"].astype(bool) & np.isfinite(F) & np.isfinite(t)
    bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
    regs = [np.ones(ok.sum()), bexc[ok], R["dx"][ok] - np.nanmean(R["dx"][ok]), R["dy"][ok] - np.nanmean(R["dy"][ok])]
    if mira_reg:
        Fm = R["FIT"]["B"]["F"][:, R["FIT"]["B"]["names"].index("Mira")]
        regs.append(Fm[ok] - np.nanmean(Fm[ok]))
    X = np.column_stack(regs)
    c, *_ = np.linalg.lstsq(X, F[ok], rcond=None)
    r = F[ok] - X[:, 1:] @ c[1:]
    tt = t[ok]
    # RRab leak
    ch, _, _ = harm_fit(tt, r, np.ones(len(r)), 0.648528)
    r = r - (harm_X(tt, 0.648528) @ ch - ch[0])
    k = np.floor((tt - tt[0]) * 24).astype(int); u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    return np.bincount(inv, tt) / n, np.bincount(inv, r) / n


def ls_top(t, y, pmin=0.5, pmax=14.0, n=4):
    fr = np.linspace(1 / pmax, 1 / pmin, 8000); pw = LombScargle(t, y).power(fr)
    top = []
    for i in np.argsort(pw)[::-1]:
        if len(top) >= n: break
        if all(abs(fr[i] - fr[j]) > 0.05 * fr[j] for j in top): top.append(i)
    return [(float(1 / fr[i]), float(pw[i])) for i in top], fr, pw


def sin_vec(t, y, P):
    X = np.column_stack([np.ones(len(t)), np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P), (t - t.mean()) / 10, ((t - t.mean()) / 10) ** 2])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return c[1:3]


def null_maxpow(t, y, nsh=100):
    """Highest LS peak of 1-d block-shuffled versions of y (destroys periods > ~1 d, keeps the noise distribution)."""
    kb = np.floor((t - t.min())).astype(int); ub = np.unique(kb)
    fr = np.linspace(1 / 14.0, 1 / 0.5, 3000); mx = []
    for _ in range(nsh):
        perm = RNG.permutation(len(ub)); y3 = np.concatenate([y[kb == ub[i]] for i in perm])
        mx.append(LombScargle(t, y3).power(fr).max())
    return np.array(mx)


def main():
    from p05_phot import setup, photometry
    real, inj = {}, []
    for s in SECTORS:
        S = setup(s); C, W, J, g, prf, rBA = S
        L = np.load(f"{PHOT}/lc_allstars_s{s:03d}.npz")
        R0 = photometry(s, S=S, verbose=False)
        m3 = g[g.name == "M3V"].iloc[0]
        img = prf.image(m3.x, m3.y)
        series0 = {}
        rs = {}
        for key, (mod, star, mreg) in {"M3V_B": ("B", "M3V", False), "M3V_B_mirareg": ("B", "M3V", True), "M3V_A": ("A", "M3V", False),
                                      "pair_A": ("A", "pair", False), "Mira_B": ("B", "Mira", False)}.items():
            tb, yb = lc_from_R(R0, L, mod, star, mreg)
            series0[key] = (tb, yb)
            top, fr, pw = ls_top(tb, yb)
            rs[key] = dict(top=top, rms_1h_e_s=float(np.std(yb)))
        nm = null_maxpow(*series0["M3V_B"]); nmA = null_maxpow(*series0["M3V_A"])
        rs["null_maxpow_p99_M3V_B"] = float(np.percentile(nm, 99)); rs["null_maxpow_p99_M3V_A"] = float(np.percentile(nmA, 99))
        real[s] = rs
        print(f"S{s} real: " + "; ".join(f"{k} " + ", ".join(f"{p:.2f}({w:.2f})" for p, w in v['top'][:3]) for k, v in rs.items() if isinstance(v, dict)) +
              f"; null p99 B {rs['null_maxpow_p99_M3V_B']:.2f} A {rs['null_maxpow_p99_M3V_A']:.2f}", flush=True)
        t_all = R0["t"]
        for P in PGRID:
            for A in AGRID:
                ph = RNG.uniform(0, 2 * np.pi)
                amp = np.nan_to_num(A * np.sin(2 * np.pi * t_all / P + ph))
                R1 = photometry(s, S=S, verbose=False, inject=[dict(img=img, amp=-amp)])
                row = dict(sector=s, P=P, A_inj=A)
                for key, (mod, star) in {"M3V_B": ("B", "M3V"), "M3V_A": ("A", "M3V")}.items():
                    tb, yb = lc_from_R(R1, L, mod, star)
                    t0_, y0_ = series0[key]
                    v1 = sin_vec(tb, yb, P); v0 = sin_vec(t0_, y0_, P)
                    top, fr, pw = ls_top(tb, yb)
                    thr = rs["null_maxpow_p99_M3V_B"] if key == "M3V_B" else rs["null_maxpow_p99_M3V_A"]
                    det = (abs(top[0][0] - P) / P < 0.1) and top[0][1] > thr
                    row[f"{key}_A_rec"] = float(np.hypot(*(v1 - v0))); row[f"{key}_top"] = top[0]; row[f"{key}_det"] = bool(det)
                inj.append(row)
            print(f"  S{s} P {P}: " + "; ".join(f"A {r['A_inj']:.0f} rec B {r['M3V_B_A_rec']:.1f} det {int(r['M3V_B_det'])} / A {r['M3V_A_A_rec']:.1f} det {int(r['M3V_A_det'])}" for r in inj if r['sector'] == s and r['P'] == P), flush=True)
        save_json({str(k): v for k, v in real.items()}, f"{ROOT}/out/rotation_real.json")
        pd.DataFrame(inj).to_csv(f"{ROOT}/out/rotation_inject.csv", index=False)
    # summary: smallest A detected per (P) in >= 4 of 6 sectors
    D = pd.DataFrame(inj)
    summ = {}
    for P in PGRID:
        for key in ("M3V_B", "M3V_A"):
            fr_ = {A: float(D[(D.P == P) & (D.A_inj == A)][f"{key}_det"].mean()) for A in AGRID}
            rec = {A: float(np.median(D[(D.P == P) & (D.A_inj == A)][f"{key}_A_rec"] / A)) for A in AGRID}
            summ[f"{key}_P{P}"] = dict(det_frac=fr_, median_recovery_ratio=rec)
    save_json(summ, f"{ROOT}/out/rotation_inject_summary.json")
    print(json.dumps(summ, indent=0))


if __name__ == "__main__":
    main()
