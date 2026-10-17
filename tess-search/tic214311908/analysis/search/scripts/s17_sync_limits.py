"""s17: limits on photometric modulation of the M3V at the orbital period (5.6915 d) and twice it (11.383 d), i.e. on
tidally synchronised rotation, and at the TARS period (~8.7 d).

Per sector: the undetrended deblended M3V light curve (model B and model A; s14 lc_from_R: regression on background
and pointing, RRab leak removed, 1-h bins, 5.69-d windows masked) is fitted with a sinusoid + quadratic trend at
P = 5.6915, 11.383 and the TARS period of that sector; the same fit at 200 fake periods (3-14 d, avoiding +-8 % of
the tested periods) gives the empirical distribution of the fitted semi-amplitude; reported: measured amplitude,
its empirical p-value, and the 95 % upper limit (measured + 1.645 x scatter, in e-/s and % of 979 e-/s).
The pixel-injection recovery factors at P = 6 and 12 d (s14) are applied when available.
Output: out/sync_limits.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/search/scripts")
sys.path.insert(1, "tic214311908/analysis/phot/scripts")
from slib import *
from s14_rotation import lc_from_R, sin_vec
RNG = np.random.default_rng(17)


def main():
    from p05_phot import setup, photometry
    pp = json.load(open(f"{ROOT}/out/public_periodograms.json"))
    out = {}
    for s in SECTORS:
        S = setup(s)
        L = np.load(f"{PHOT}/lc_allstars_s{s:03d}.npz")
        R = photometry(s, S=S, verbose=False)
        Ptars = pp.get(f"TARS_M3V_S{s}", {}).get("top", [[8.7]])[0][0]
        res = {}
        for key, (mod, star) in {"M3V_B": ("B", "M3V"), "M3V_A": ("A", "M3V"), "pair_A": ("A", "pair")}.items():
            tb, yb = lc_from_R(R, L, mod, star)
            r = {}
            for nm, P in (("P_orb", P_B), ("2P_orb", 2 * P_B), ("P_TARS", Ptars)):
                a = float(np.hypot(*sin_vec(tb, yb, P)))
                fk = []
                while len(fk) < 200:
                    pf = RNG.uniform(3.0, 14.0)
                    if min(abs(pf - P_B), abs(pf - 2 * P_B), abs(pf - Ptars)) / pf < 0.08: continue
                    fk.append(np.hypot(*sin_vec(tb, yb, pf)))
                fk = np.array(fk)
                r[nm] = dict(P=float(P), amp_e_s=a, p_emp=float((fk >= a).mean()), fake_median=float(np.median(fk)),
                             fake_p95=float(np.percentile(fk, 95)), ul95_e_s=float(max(a, np.median(fk)) + 1.645 * np.std(fk)),
                             ul95_pct=float((max(a, np.median(fk)) + 1.645 * np.std(fk)) / 979 * 100))
            res[key] = r
        out[s] = res
        print(f"S{s}: " + " | ".join(f"{k} " + ", ".join(f"{nm} {v['amp_e_s']:.1f} e/s (p {v['p_emp']:.2f}, UL95 {v['ul95_pct']:.1f}%)" for nm, v in r.items()) for k, r in res.items()), flush=True)
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/sync_limits.json")


if __name__ == "__main__":
    main()
