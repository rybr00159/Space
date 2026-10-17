"""r08: (a) can an equal-eclipse EB at 2P = 11.383 d mimic the signal?  (b) X-ray activity vs rotation.

(a1) Odd - even timing offset (a secondary eclipse displaced from phase 0.5 would show up as an odd/even O-C offset):
     weighted mean O-C of odd and even epochs (phot model-B strict times with err_cal, and the referee's times with
     tc_err_cal) about a common linear ephemeris -> e cos(w) limit for a 2P binary: dt = (2 * 2P / pi) e cos w.
(a2) Main-sequence pairs: for a flat-bottomed (total/annular) eclipse the radius ratio must satisfy
     k <= (1 - f)/(1 + f) (b = 0; smaller for b > 0), f = T23/T14 >= 0.47 (phot) or ~0.6-0.73 (referee) ->
     k <= 0.36 (or <= 0.25). For primaries along an approximate dwarf sequence (Pecaut & Mamajek-type Teff-R table,
     hard-coded below) and secondaries of radius k R1 on the same sequence, the occultation/transit depth ratio
     ~ S2/S1 / 1.1 (S = blackbody photon surface brightness integrated over 600-1000 nm; 1.1 = limb-darkening
     enhancement of a central transit). The predicted odd-even difference for a 2P EB with the observed mean depth
     (15.7 ppt) is compared with the measured odd-even values.
(b)  L_X/L_bol = 4.2e-4 (eRASS1, 7 +- 3 counts): Rossby number and P_rot from Wright et al. 2011 (Ro_sat 0.13,
     log Rx_sat -3.13, beta -2.18, tau(M) of W11) and Wright et al. 2018 (Ro_sat 0.14, log Rx_sat -3.05, beta -2.3,
     log tau = 2.33 - 1.50 M + 0.31 M^2), M = 0.37 M_sun; with +-0.3 dex (scatter) and +-0.19 dex (counting) shifts.
Output: out/r08_2p_xray.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *

# approximate dwarf sequence (Teff [K], R [R_sun]); Pecaut & Mamajek 2013 style values (approximate)
SEQ = np.array([[6550, 1.33], [5930, 1.10], [5660, 0.92], [5270, 0.83], [4830, 0.76], [4440, 0.68], [4050, 0.63],
                [3850, 0.59], [3660, 0.50], [3560, 0.44], [3430, 0.37], [3210, 0.27], [3060, 0.20], [2810, 0.14],
                [2680, 0.12], [2570, 0.11], [2380, 0.10]])


def teff_of_R(R):
    Rs, Ts = SEQ[::-1, 1], SEQ[::-1, 0]
    return np.interp(R, Rs, Ts, left=2000.0)


def S_tess(T):
    lam = np.linspace(0.6, 1.0, 200)          # micron
    x = 14387.77 / (lam * T)
    return np.trapezoid(1 / lam ** 4 / np.expm1(x), lam)   # photon-number weighting ~ B_lambda * lambda


def main():
    out = {}
    # ---- a1 odd/even timing
    T = pd.read_csv(f"{BASE}/phot/out/timing_table.csv"); T = T[(T.model == "B") & (T.relaxed == False) & T.tc.notna()]
    R = pd.read_csv(f"{BASE}/referee_phot/out/per_transit.csv"); R = R[(R.sel == "strict") & R.tc.notna() & (R.epoch != -388)]
    tim = {}
    for nm, (E, tc, err) in dict(phot=(T.epoch.values, T.tc.values, T.err_cal.values),
                                 referee=(R.epoch.values, R.tc.values, R.tc_err_cal.values)).items():
        E = E.astype(int); odd = (E % 2) != 0
        # model: tc = T0 + P E + delta * (odd - 0.5)
        X = np.column_stack([np.ones(len(E)), E, odd - 0.5])
        wv = 1 / err ** 2
        cov = np.linalg.inv((X * wv[:, None]).T @ X); c = cov @ ((X * wv[:, None]).T @ tc)
        chi2 = float(np.sum(wv * (tc - X @ c) ** 2))
        dlt = c[2] * 1440; dle = np.sqrt(cov[2, 2]) * 1440
        tim[nm] = dict(odd_minus_even_min=float(dlt), err_min=float(dle), chi2=chi2, n=len(E),
                       ecosw_2P=float(c[2] * np.pi / (2 * 2 * P)), ecosw_2P_err=float(np.sqrt(cov[2, 2]) * np.pi / (2 * 2 * P)))
    out["oddeven_timing"] = tim
    print(tim)
    # ---- a2 MS pairs
    depth = 15.7
    obs = {"model-B cleaned (r01 own)": (2.67, 1.71), "referee raw per-transit": (1.57, 1.75),
           "referee bias-subtracted": (-0.60, 1.75), "search alias.json": (3.27, 1.78)}
    rows = []
    for T1, R1 in SEQ:
        for kmax in (0.36, 0.25):
            R2 = kmax * R1
            T2 = teff_of_R(R2)
            ratio = S_tess(T2) / S_tess(T1) / 1.1
            # depths: d_pri = 2 depth / (1 + ratio), d_sec = ratio d_pri
            dp = 2 * depth / (1 + ratio); ds = ratio * dp
            rows.append(dict(T1=T1, R1=R1, kmax=kmax, R2=R2, T2=float(T2), sec_over_pri=float(ratio), pred_odd_minus_even_ppt=float(dp - ds),
                             **{f"sigma_vs_{k}": float((dp - ds - v[0]) / v[1]) for k, v in obs.items()}))
    D = pd.DataFrame(rows)
    out["ms_pairs"] = D.to_dict(orient="records")
    print(D[["T1", "R1", "kmax", "T2", "sec_over_pri", "pred_odd_minus_even_ppt", "sigma_vs_model-B cleaned (r01 own)", "sigma_vs_referee raw per-transit"]].round(2).to_string())
    out["ms_pairs_min_sigma"] = {k: float(D[f"sigma_vs_{k}"].min()) for k in obs}
    # ---- b X-ray / rotation
    M = 0.37; Rx = 4.2e-4
    tau11 = 10 ** (1.16 - 1.49 * np.log10(M) - 0.54 * np.log10(M) ** 2)
    tau18 = 10 ** (2.33 - 1.50 * M + 0.31 * M ** 2)
    xr = {}
    for nm, (Ros, lRs, beta, tau) in dict(W11=(0.13, -3.13, -2.18, tau11), W18=(0.14, -3.05, -2.3, tau18)).items():
        res = {}
        for sh in (-0.3, -0.19, 0.0, 0.19, 0.3):
            lR = np.log10(Rx) + sh
            if lR >= lRs:
                res[f"{sh:+.2f}"] = dict(saturated=True, Ro_max=Ros, Prot_max_d=float(Ros * tau))
            else:
                Ro = Ros * 10 ** ((lR - lRs) / beta)
                res[f"{sh:+.2f}"] = dict(saturated=False, Ro=float(Ro), Prot_d=float(Ro * tau))
        xr[nm] = dict(tau_d=float(tau), by_logRx_shift=res)
    out["xray_rotation"] = xr
    print(json.dumps(xr, indent=0))
    jdump(out, f"{ME}/out/r08_2p_xray.json")


if __name__ == "__main__":
    main()
