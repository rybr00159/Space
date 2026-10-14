"""Per-sector eclipse-depth scale (primary and secondary separately) with the geometry fixed at the best fit; per-event
quadratic baselines profiled.  Error from delta-chi2 = 1 and from the scatter of per-event depth scales (empirical)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L, fitlib as F
p = json.load(open(L.HERE + '/results/ml_profile_A_best.json'))['free']
for s in fit3.SECS: p[f's{s}'] = 1.0
out = {}
for name, path in [('default', L.HERE + '/results/data_default.pkl'), ('tc', L.HERE + '/results/data_tc.pkl')]:
    d = fit3.load_data(path)
    m = fit3.model_all(d, p); dm = 1 - m
    w = 1 / d['sig'] ** 2
    res = {}
    for e in range(len(d['starts'])):
        a, b = d['starts'][e], d['ends'][e]; ev = d['evtab'][e]
        dt = d['dt'][a:b]; y = d['y'][a:b]
        X = np.vstack([np.ones_like(dt), dt, dt * dt, -dm[a:b]]).T * np.sqrt(w[a:b])[:, None]
        c, *_ = np.linalg.lstsq(X, y * np.sqrt(w[a:b]), rcond=None)
        cov = np.linalg.inv(X.T @ X)
        res.setdefault((ev['sector'], ev['ph']), []).append((c[3] / c[0], np.sqrt(cov[3, 3]) / c[0]))
    for (sec, ph), v in sorted(res.items()):
        v = np.array(v); v = v[np.isfinite(v[:, 0]) & (v[:, 1] < 1.0)]
        wt = 1 / v[:, 1] ** 2; mu = np.sum(v[:, 0] * wt) / wt.sum(); ef = 1 / np.sqrt(wt.sum())
        emp = np.std(v[:, 0], ddof=1) / np.sqrt(len(v))
        out[f'{name}_S{sec}_ph{ph}'] = dict(scale=float(mu), err_formal=float(ef), err_emp=float(emp), n=len(v))
        print(f'{name:8s} S{sec:3d} ph{ph}: depth scale {mu:.3f} +- {ef:.3f} (formal) +- {emp:.3f} (event scatter), N={len(v)}')
    for ph in (0.0, 0.5):
        vals = np.array([[out[f'{name}_S{s}_ph{ph}']['scale'], max(out[f'{name}_S{s}_ph{ph}']['err_formal'], out[f'{name}_S{s}_ph{ph}']['err_emp'])] for s in fit3.SECS])
        wt = 1 / vals[:, 1] ** 2; mu = np.sum(vals[:, 0] * wt) / wt.sum(); chi = np.sum((vals[:, 0] - mu) ** 2 * wt)
        print(f'  {name} ph{ph}: weighted mean {mu:.3f}, chi2 of constancy {chi:.1f} for 5 dof')
        out[f'{name}_ph{ph}_constancy'] = dict(mean=float(mu), chi2=float(chi), dof=5)
json.dump(out, open(L.HERE + '/results/depth_by_sector.json', 'w'), indent=1)
