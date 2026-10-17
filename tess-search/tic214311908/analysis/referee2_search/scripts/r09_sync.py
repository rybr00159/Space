"""r09: modulation of the deblended M3V flux at P_orb = 5.6915 d, 2 P_orb = 11.383 d and the TARS period (8.7 d).

Per sector: undetrended model-B flux (phot flux_raw_norm, rescaled to 979 e-/s), 5.69-d windows removed, RRab model
(r01) subtracted, 1-h bins; least-squares fit of quadratic trend + sinusoid at the trial period; semi-amplitude A in %.
Empirical calibration: the same fit at 150 fake periods drawn log-uniformly within a factor 1.4 of the trial period
(excluding +-5 % around P_orb, 2 P_orb and 8.7 d); p = fraction of fakes with amplitude >= observed; a 95 % upper
limit is the observed amplitude plus the 90th-percentile spread of the fake amplitudes (conservative, approximately).
Also: S12 and S93 at the TARS periods 8.1-8.8 d (amplitude in e-/s).
Output: out/r09_sync.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
rng = np.random.default_rng(909)


def fit(t, y, Pp):
    tt = (t - t.mean()) / (t.max() - t.min())
    X = np.column_stack([np.ones(len(t)), tt, tt ** 2, np.sin(2 * np.pi * t / Pp), np.cos(2 * np.pi * t / Pp)])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(np.hypot(c[3], c[4]))


def main():
    d = load_m3v("B")
    z = np.load(f"{ME}/data/m3v_clean_r2.npz")
    rawc = d.raw.values - z["model"]          # same row order (both from load_m3v)
    out = {}
    for sec in SECTORS:
        m = (d.s.values == sec) & ~d.intr.values
        tb, yb, nb, sb = binlc(d.t.values[m], rawc[m] - 1, d.s.values[m], d.cad.values[m], 60.0)
        span = tb.max() - tb.min()
        res = dict(span_d=float(span))
        for nm, Pp in (("Porb", P), ("2Porb", 2 * P), ("TARS_8p7", 8.7)):
            A = fit(tb, yb, Pp)
            fk = []
            while len(fk) < 150:
                pf = Pp * np.exp(rng.uniform(-np.log(1.4), np.log(1.4)))
                if min(abs(pf / P - 1), abs(pf / (2 * P) - 1), abs(pf / 8.7 - 1)) < 0.05: continue
                fk.append(fit(tb, yb, pf))
            fk = np.array(fk)
            res[nm] = dict(amp_pct=A * 100, amp_e_s=A * FREF, p_emp=float(np.mean(fk >= A)),
                           fake_median_pct=float(np.median(fk) * 100), fake_p90_pct=float(np.percentile(fk, 90) * 100),
                           ul95_pct=float((A + np.percentile(fk, 90)) * 100), periods_per_span=float(span / Pp))
        out[sec] = res
        print(sec, {k: (v if not isinstance(v, dict) else {kk: round(vv, 3) for kk, vv in v.items()}) for k, v in res.items()}, flush=True)
    jdump(out, f"{ME}/out/r09_sync.json")


if __name__ == "__main__":
    main()
