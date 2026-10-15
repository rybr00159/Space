"""Convergence checks for the v6 MS chains: per-parameter integrated autocorrelation time (post-burn), split-half
comparison of medians/68% widths, and a rank-normalised split R-hat over 4 walker groups x 2 halves."""
import numpy as np, json, emcee, sys
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v7rep'
res = {}
for tag in sys.argv[1:]:
    meta = json.load(open(f'{R}/chain_{tag}_meta.json'))
    rd = emcee.backends.HDFBackend(meta['chain'], read_only=True)
    n = rd.iteration; b = int(float(__import__("os").environ.get("BURN", "0.3")) * n)
    ch = rd.get_chain(discard=b)           # steps, walkers, dim
    lp = rd.get_log_prob(discard=b)
    tau = emcee.autocorr.integrated_time(ch, tol=0)
    free = meta['free']; out = {}
    h = ch.shape[0] // 2
    for j, nme in enumerate(free):
        x1 = ch[:h, :, j].ravel(); x2 = ch[h:, :, j].ravel()
        q1 = np.percentile(x1, [16, 50, 84]); q2 = np.percentile(x2, [16, 50, 84])
        # split R-hat: 4 walker groups x 2 halves = 8 chains
        chains = []
        for g in np.array_split(np.arange(ch.shape[1]), 4):
            chains += [ch[:h, g, j].ravel(), ch[h:2 * h, g, j].ravel()]
        m = np.array([c.mean() for c in chains]); v = np.array([c.var(ddof=1) for c in chains]); L = len(chains[0])
        W = v.mean(); B = L * m.var(ddof=1); rhat = np.sqrt(((L - 1) / L * W + B / L) / W)
        out[nme] = dict(tau=float(tau[j]), n_over_tau=float(n / tau[j]), med_half1=float(q1[1]), med_half2=float(q2[1]),
                        shift_in_sigma=float((q2[1] - q1[1]) / (0.5 * (q1[2] - q1[0]))), rhat=float(rhat))
    lpm = lp.mean(axis=1)
    res[tag] = dict(iterations=int(n), burn=b, nwalkers=int(ch.shape[1]), params=out, max_tau=float(tau.max()),
                    lp_mean_first_last_quarter=[float(lpm[:len(lpm) // 4].mean()), float(lpm[-len(lpm) // 4:].mean())],
                    acceptance=float(np.mean(rd.accepted / n)))
    print(tag, 'iter', n, 'burn', b, 'max tau', round(tau.max(), 1), 'n/tau_max', round(n / tau.max(), 1), 'acc', round(res[tag]['acceptance'], 3),
          'mean lp first/last quarter', [round(x, 2) for x in res[tag]['lp_mean_first_last_quarter']])
    for nme in ('M1', 'q', 'L3', 'dmag', 'k', 'rsum', 'J', 'b0', 'esinw', 'ecosw', 'dR1', 'dR2', 'dJ', 'k', 'rsum', 'J'):
        if nme in out: print(f"   {nme:6s} tau {out[nme]['tau']:6.1f} halves {out[nme]['med_half1']:.4f} {out[nme]['med_half2']:.4f} shift {out[nme]['shift_in_sigma']:+.2f} sig  Rhat {out[nme]['rhat']:.3f}")
    print('   max Rhat', round(max(o['rhat'] for o in out.values()), 4), 'max |shift|', round(max(abs(o['shift_in_sigma']) for o in out.values()), 2))
json.dump(res, open(f'{R}/convergence_{"_".join(sys.argv[1:])}.json', 'w'), indent=1)
