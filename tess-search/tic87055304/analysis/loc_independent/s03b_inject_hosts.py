"""Injection calibration for THREE candidate hosts (PRF injection, PRF analysis): the target, the
neighbour TIC 87055315 and Gaia DR3 6720986208880830336 ('C3', G=17.96, 3.4" from the neighbour,
the only other star the single-host chi2 scan ever prefers).  Same machinery as s03 ('resampled' data
sets).  For each sector the 2-D distribution of the fitted sky offset (E,N) under each host is
approximated by a Gaussian; summing log-likelihoods over sectors gives ln L(host) for the real data.
Outputs: inj3_prf.json, inj3_prf.npz
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pipeline import *

NSIM = 600
rng = np.random.default_rng(19)
real = json.load(open(f'{HERE}/diff_prf_main.json'))
res = {}; arrays = {}
tot = {'T': 0.0, 'N': 0.0, 'C3': 0.0}
for sec in SECTORS:
    ctx = SectorContext(sec, psfkind='prf', reg_variant='local_r4', rfit=3.0)
    rw = ctx.real_windows(); nT = len(rw)
    A = real[str(sec)]['stacked']['free_amp']
    amp = A / np.mean([ctx.in_window_shape_mean(w) for w in rw])
    fw = ctx.fake_windows()
    hosts = {'T': ctx.pT, 'N': ctx.pN, 'C3': np.array([ctx.x[ctx.ic3], ctx.y[ctx.ic3]])}
    Dn, Sg, Dh = [], [], {h: [] for h in hosts}
    for w in fw:
        D, sg, _ = ctx.diff(w); Dn.append(D); Sg.append(sg)
        for h, p in hosts.items():
            Dh[h].append(ctx.inject_diff(w, p, amp))
    Dn = np.array(Dn); Sg = np.array(Sg); Dh = {h: np.array(v) for h, v in Dh.items()}
    sets = [rng.choice(len(fw), size=min(nT, len(fw)), replace=False) for _ in range(NSIM)]
    xr = np.array([real[str(sec)]['stacked']['offT_E'], real[str(sec)]['stacked']['offT_N']])
    out = {}
    for h in hosts:
        vals = []
        for k in sets:
            Ds, ss = stack(Dn[k] + Dh[h][k], Sg[k], mask=ctx.mask)
            a = ctx.analyse(Ds, ss, free=True, extra=False)
            vals.append([a['offT_E'], a['offT_N'], a['s_par'], a['perp_arcsec']])
        V = np.array(vals); arrays[f's{sec}_{h}'] = V
        ok = np.all(np.isfinite(V), axis=1); V = V[ok]
        mu = V[:, :2].mean(0); C = np.cov(V[:, :2].T)
        d = xr - mu
        lnL = float(-0.5 * d @ np.linalg.solve(C, d) - 0.5 * np.log(np.linalg.det(2 * np.pi * C)))
        out[h] = dict(mean_EN=mu.tolist(), cov_EN=C.tolist(), mahal=float(np.sqrt(d @ np.linalg.solve(C, d))), lnL=lnL,
                      s_par_mean=float(V[:, 2].mean()), s_par_std=float(V[:, 2].std()),
                      perp_mean=float(V[:, 3].mean()), perp_std=float(V[:, 3].std()))
        tot[h] += lnL
    res[str(sec)] = out
    print(f"S{sec}: real (E,N)=({xr[0]:+.2f},{xr[1]:+.2f})  " + '  '.join(
        f"{h}: mean=({o['mean_EN'][0]:+.2f},{o['mean_EN'][1]:+.2f}) mahal={o['mahal']:.1f} lnL={o['lnL']:.1f}" for h, o in out.items()), flush=True)
res['sum_lnL'] = tot
res['lnLR_N_over_T'] = tot['N'] - tot['T']
res['lnLR_N_over_C3'] = tot['N'] - tot['C3']
print('sum lnL:', tot, ' lnLR(N/T)=%.1f lnLR(N/C3)=%.1f' % (tot['N'] - tot['T'], tot['N'] - tot['C3']))
to_json(res, f'{HERE}/inj3_prf.json')
np.savez_compressed(f'{HERE}/inj3_prf.npz', **arrays)
