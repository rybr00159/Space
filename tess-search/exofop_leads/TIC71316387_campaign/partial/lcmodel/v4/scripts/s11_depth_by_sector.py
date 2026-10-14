"""Per-sector, per-product eclipse-depth scale (primary and secondary separately) with the shape fixed at the best fit
(cadence-integrated), per-event quadratic baselines profiled.  Errors: formal (Delta chi2 = 1) and empirical (scatter of
per-event depth scales / sqrt(N)); the larger is used.  Constancy chi2 across sectors per product set, and the secondary/
primary depth-ratio constancy (independent of any dilution, which scales both eclipses equally)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '2'
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
p = json.load(open(f'{FD}/results/best_geo.json'))
for s in fit3.SECS: p[f's{s}'] = 1.0
out = {}
for name in ('default', 'tc', 'qlp'):
    d = fit3.load_data(f'{FD}/results/data_{name}.pkl')
    dm = 1 - fit3.model_all(d, p); w = 1 / d['sig'] ** 2
    res = {}
    for e in range(len(d['starts'])):
        a, b = d['starts'][e], d['ends'][e]; ev = d['evtab'][e]
        dt = d['dt'][a:b]; y = d['y'][a:b]
        X = np.vstack([np.ones_like(dt), dt, dt * dt, -dm[a:b]]).T * np.sqrt(w[a:b])[:, None]
        c, *_ = np.linalg.lstsq(X, y * np.sqrt(w[a:b]), rcond=None); cov = np.linalg.inv(X.T @ X)
        res.setdefault((ev['sector'], ev['ph'], ev['prod']), []).append((c[3] / c[0], np.sqrt(cov[3, 3]) / c[0]))
    rows = {}
    for (sec, ph, prod), v in sorted(res.items()):
        v = np.array(v); wt = 1 / v[:, 1] ** 2; mu = np.sum(v[:, 0] * wt) / wt.sum(); ef = 1 / np.sqrt(wt.sum())
        emp = np.std(v[:, 0], ddof=1) / np.sqrt(len(v)); err = max(ef, emp)
        rows[(sec, ph)] = (mu, err, prod)
        out[f'{name}_S{sec}_{prod}_ph{ph}'] = dict(scale=float(mu), err_formal=float(ef), err_emp=float(emp), n=len(v))
    print(f'--- {name}')
    for sec in fit3.SECS:
        if (sec, 0.0) not in rows: continue
        mp, ep, prod = rows[(sec, 0.0)]; ms, esd, _ = rows[(sec, 0.5)]
        rat = ms / mp; erat = rat * np.hypot(ep / mp, esd / ms)
        out[f'{name}_S{sec}_ratio'] = dict(ratio=float(rat), err=float(erat))
        print(f'S{sec:3d} {prod:9s}: primary depth scale {mp:.3f} +- {ep:.3f} | secondary {ms:.3f} +- {esd:.3f} | sec/prim (relative to model) {rat:.3f} +- {erat:.3f}')
    for ph in (0.0, 0.5, 'ratio'):
        if ph == 'ratio':
            v = np.array([[out[f'{name}_S{s}_ratio']['ratio'], out[f'{name}_S{s}_ratio']['err']] for s in fit3.SECS if f'{name}_S{s}_ratio' in out])
        else:
            v = np.array([[rows[(s, ph)][0], rows[(s, ph)][1]] for s in fit3.SECS if (s, ph) in rows])
        wt = 1 / v[:, 1] ** 2; mu = np.sum(v[:, 0] * wt) / wt.sum(); chi = float(np.sum((v[:, 0] - mu) ** 2 * wt))
        out[f'{name}_{ph}_constancy'] = dict(mean=float(mu), chi2=chi, dof=len(v) - 1)
        print(f'  {name} {ph}: weighted mean {mu:.3f}; chi2 of constancy {chi:.1f} / {len(v)-1}')
json.dump(out, open(f'{FD}/results/depth_by_sector.json', 'w'), indent=1)
