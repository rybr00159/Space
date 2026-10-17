"""c04: global transit fits, per-sector and odd/even depths on the CLEANED light curves of both pipelines.

Data: data/win_inv_raw.npz (INV, undetrended flux, +-0.30 d windows, per-window offset + slope) and data/win_ref.npz
(REF strict windows). Per-point sigma = window OOT robust std x beta_sector (fake-epoch calibration, out/noise_calib.json).
Model: batman, quadratic LD with Gaussian priors u1 = 0.125 +- 0.1, u2 = 0.533 +- 0.1 (physical constraints), circular
orbit, exposure-integrated (midpoint rule). Parameters: T0 (epoch 0), P, k, b, rho* (g/cm3), u1, u2; per-window
multiplicative linear baselines profiled analytically.
Priors: 'single' = rho* ~ N(9.3, 0.9) (single M3V); 'free' = log-uniform rho* in [0.3, 100].
Sampler: emcee, 40 walkers x 6000 steps, burn 2000, thin 10.
Then, with the shape fixed at the posterior median of the 'single' fit: per-sector depths, early/late groups, odd/even
(linear fits; errors = formal with beta-scaled sigma, i.e. fake-epoch calibrated), and the same at fake ephemerides as
a check of the error scale.
Outputs: out/global_{inv,ref}_{single,free}.json, out/chains_*.npz (thinned), out/depths_sector.json
"""
import sys, os
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import CL, P1, T01, F_M3V, SECTORS, robust_std, jd
from tfit import *
from c03_windows import load_windows
import emcee
from scipy.optimize import minimize

NOISE = json.load(open(f"{CL}/out/noise_calib.json")) if False else None
import json as _json
NOISE = _json.load(open(f"{CL}/out/noise_calib.json"))
TEFF = (3350.0, 60.0)


def load(pipe):
    if pipe == "inv":
        W = load_windows(f"{CL}/data/win_inv_raw.npz"); cal = NOISE["inv_raw"]
    elif pipe == "ref":
        W = load_windows(f"{CL}/data/win_ref.npz"); cal = NOISE["ref"]
    elif pipe == "reforig":
        W = load_windows(f"{CL}/data/win_reforig.npz"); cal = NOISE["reforig"]
    beta = np.array([cal[str(s)]["beta"] for s in W.sector])
    W.s = W.s * np.maximum(beta[W.w], 1.0)
    return W


class Like:
    def __init__(self, W, prior):
        self.W = W; self.prior = prior
        self.M = W.model()
        self.wt = 1 / W.s ** 2
        self.nw = W.nw

    def chi2(self, T0, P, k, b, aR, u1, u2):
        W = self.W; m = self.M.flux(T0, P, k, b, aR, u1, u2)
        w = self.wt; tau = W.tau; f = W.f; idx = W.w
        S00 = np.bincount(idx, w * m * m, self.nw); S01 = np.bincount(idx, w * m * m * tau, self.nw)
        S11 = np.bincount(idx, w * m * m * tau * tau, self.nw)
        b0 = np.bincount(idx, w * m * f, self.nw); b1 = np.bincount(idx, w * m * tau * f, self.nw)
        det = S00 * S11 - S01 ** 2
        q = (S11 * b0 ** 2 - 2 * S01 * b0 * b1 + S00 * b1 ** 2) / det
        return float(np.sum(w * f * f) - np.sum(q))

    def lnprob(self, p):
        T0, P, k, b, lrho, u1, u2 = p
        if not (0.02 < k < 0.4 and 0 <= b < 1 + k and -1.2 < lrho < 2.0): return -np.inf
        if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0): return -np.inf
        if abs(P - P1) > 1e-3 or abs(T0 - T01) > 0.02: return -np.inf
        rho = 10 ** lrho
        lp = -0.5 * ((u1 - LD0[0]) / 0.1) ** 2 - 0.5 * ((u2 - LD0[1]) / 0.1) ** 2
        if self.prior == "single":
            lp += -0.5 * ((rho - RHO_STAR[0]) / RHO_STAR[1]) ** 2 + np.log(rho)      # Gaussian in rho (sampled in log rho)
        aR = a_over_r(rho, P)
        if b >= aR: return -np.inf
        return lp - 0.5 * self.chi2(T0, P, k, b, aR, u1, u2)


def run_fit(pipe, prior, nwalk=40, nstep=6000, burn=2000, seed=1):
    W = load(pipe); L = Like(W, prior)
    p0 = np.array([T01, P1, 0.12, 0.5, np.log10(9.3), LD0[0], LD0[1]])
    nll = lambda p: -L.lnprob(p) if np.isfinite(L.lnprob(p)) else 1e25
    best = None
    for b0 in (0.2, 0.5, 0.75):
        for lr in (np.log10(5), np.log10(9.3), np.log10(15)):
            q = minimize(nll, np.r_[p0[:3], b0, lr, p0[5:]], method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-7, fatol=1e-3))
            if best is None or q.fun < best.fun: best = q
    rng = np.random.default_rng(seed)
    sc = np.array([5e-4, 5e-6, 3e-3, 0.05, 0.03, 0.03, 0.03])
    pos = []
    while len(pos) < nwalk:
        x = best.x + sc * rng.normal(size=7)
        if np.isfinite(L.lnprob(x)): pos.append(x)
    sam = emcee.EnsembleSampler(nwalk, 7, L.lnprob)
    sam.run_mcmc(np.array(pos), nstep, progress=False)
    ch = sam.get_chain(discard=burn, thin=10, flat=True)
    lp = sam.get_log_prob(discard=burn, thin=10, flat=True)
    T0, P, k, b, lr, u1, u2 = ch.T
    rho = 10 ** lr; aR = a_over_r(rho, P)
    inc = np.degrees(np.arccos(b / aR))
    T14 = t14(P, k, b, aR) * 24; T23 = t23(P, k, b, aR) * 24
    dep = np.array([mid_depth(*x) for x in np.column_stack([k, b, aR, u1, u2])[::max(1, len(k) // 3000)]])
    rng2 = np.random.default_rng(5)
    Rs = rng2.normal(*R_STAR, len(k)); Te = rng2.normal(*TEFF, len(k))
    Rp = k * Rs / R_EARTH_RSUN; Teq = Te * np.sqrt(1 / (2 * aR))
    def summ(x):
        q = np.percentile(x, [2.5, 16, 50, 84, 97.5]); return dict(med=float(q[2]), lo=float(q[2] - q[1]), hi=float(q[3] - q[2]), p2_5=float(q[0]), p97_5=float(q[4]))
    imax = int(np.argmax(lp))
    chi2_best = L.chi2(T0[imax], P[imax], k[imax], b[imax], aR[imax], u1[imax], u2[imax])
    # V-shape / grazing fraction and flat fraction
    flat = np.where(T14 > 0, T23 / T14, 0)
    out = dict(pipe=pipe, prior=prior, n_points=int(len(W.t)), n_windows=int(W.nw), epochs=W.epoch.tolist(), sectors=W.sector.tolist(),
               acceptance=float(np.mean(sam.acceptance_fraction)),
               autocorr=[float(x) for x in np.atleast_1d(sam.get_autocorr_time(discard=burn, quiet=True))],
               chi2_best=chi2_best, map=dict(T0=float(T0[imax]), P=float(P[imax]), k=float(k[imax]), b=float(b[imax]), rho=float(rho[imax]),
                                              aR=float(aR[imax]), u1=float(u1[imax]), u2=float(u2[imax])),
               T0=summ(T0), P=summ(P), k=summ(k), b=summ(b), rho=summ(rho), aR=summ(aR), inc=summ(inc), T14_h=summ(T14),
               flat_frac=summ(flat), depth_mid_ppt=summ(dep * 1e3), Rp_RE=summ(Rp), Teq_K=summ(Teq), u1=summ(u1), u2=summ(u2),
               frac_grazing=float(np.mean(b > 1 - k)), cov_T0_P=np.cov(T0, P).tolist())
    np.savez_compressed(f"{CL}/out/chains_{pipe}_{prior}.npz", chain=ch, lnp=lp)
    jd(out, f"{CL}/out/global_{pipe}_{prior}.json")
    print(f"{pipe} {prior}: k {out['k']['med']:.4f}+{out['k']['hi']:.4f}-{out['k']['lo']:.4f} b {out['b']['med']:.2f} rho {out['rho']['med']:.1f}"
          f"+{out['rho']['hi']:.1f}-{out['rho']['lo']:.1f} T14 {out['T14_h']['med']:.2f} flat {out['flat_frac']['med']:.2f} depth {out['depth_mid_ppt']['med']:.2f}"
          f"+-{0.5*(out['depth_mid_ppt']['hi']+out['depth_mid_ppt']['lo']):.2f} Rp {out['Rp_RE']['med']:.2f}+-{0.5*(out['Rp_RE']['hi']+out['Rp_RE']['lo']):.2f} "
          f"T0 {out['T0']['med']:.5f}+-{0.5*(out['T0']['hi']+out['T0']['lo']):.5f} P {out['P']['med']:.7f}+-{0.5*(out['P']['hi']+out['P']['lo']):.7f}; "
          f"grazing {out['frac_grazing']:.3f}; acc {out['acceptance']:.2f}; chi2 {chi2_best:.0f}/{len(W.t)}", flush=True)
    return out


if __name__ == "__main__":
    which = sys.argv[1:] or ["inv:single", "ref:single", "inv:free", "ref:free"]
    for w in which:
        pipe, prior = w.split(":")
        run_fit(pipe, prior)
