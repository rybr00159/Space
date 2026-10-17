"""Referee r15: is the timing scatter consistent with noise? Compares the real normalised O-C (r07 formal errors,
linear ephemeris refit) with the injection pull distribution of the same pipeline (r07 pixel injections, re-run here
to keep the pulls), via a two-sample KS test and a Monte Carlo of the chi2 statistic: draw 20 pulls per trial from
the per-sector injection pull pools (same sector mix as the real events), fit a line in epoch, compute chi2.
Output: out/timing_dist.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, fake_centres
from r07_pertransit import fit_one, T14, FL
from scipy.stats import ks_2samp
rng = np.random.default_rng(321)
PT = pd.read_csv(f"{R}/out/per_transit.csv")
tm = PT[(PT.sel == "strict") & (PT["cov"] > 0.5) & ~((PT.sector == 12) & (PT.epoch == -388))].copy()
pools = {}
for s in SECTORS:
    S = Sector(s); fk = fake_centres(S, 120, rng); p = []
    for tcw in fk:
        tt = tcw + rng.uniform(-0.015, 0.015)
        r = S.window(tcw, inject=[(S.P_M3V, 15.0, tt, T14, FL)])
        if r is None: continue
        f = fit_one(r["t"], r["aM"], S.exptime, tcw)
        if f is not None: p.append((f["tc"] - tt) / f["tc_err"])
    pools[s] = np.array(p)
    print(f"S{s}: {len(p)} pulls, std {np.std(p):.2f}, robust {rstd(p):.2f}, |>3| {np.mean(np.abs(p) > 3):.2f}", flush=True)
e = tm.tc_err.values; k = tm.epoch.values
A = np.column_stack([np.ones(len(k)), k]) / e[:, None]
c, *_ = np.linalg.lstsq(A, tm.tc.values / e, rcond=None)
real = (tm.tc.values - (c[0] + c[1] * k)) / e
chi_real = float(np.sum(real ** 2))
allp = np.concatenate([pools[s] for s in SECTORS])
ks = ks_2samp(np.abs(real), np.abs(allp))
# Monte Carlo chi2
sims = []
for i in range(5000):
    pulls = np.array([rng.choice(pools[s]) for s in tm.sector.values])
    tcs = pulls * e
    c2, *_ = np.linalg.lstsq(A, tcs / e, rcond=None)
    sims.append(np.sum((tcs / e - A @ c2) ** 2))
sims = np.array(sims)
out = dict(n=len(tm), chi2_real_formal=chi_real, mc_chi2_median=float(np.median(sims)), mc_p_value=float(np.mean(sims >= chi_real)),
           ks_stat=float(ks.statistic), ks_p=float(ks.pvalue), real_abs_pulls=np.round(np.sort(np.abs(real)), 2).tolist(),
           frac_real_gt2=float(np.mean(np.abs(real) > 2)), frac_inj_gt2=float(np.mean(np.abs(allp) > 2)))
save_json(out, f"{R}/out/timing_dist.json")
print(json.dumps({k_: v for k_, v in out.items() if k_ != 'real_abs_pulls'}, indent=0)); print("real |pulls| sorted:", out["real_abs_pulls"])
