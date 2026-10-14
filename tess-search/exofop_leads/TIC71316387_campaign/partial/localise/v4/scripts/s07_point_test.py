"""v4: calibrated tests of specific candidate positions for the eclipsing source under the chromatic model M2
(delta fixed at the candidate, colour levers free): observed dchi2 = chi2(delta fixed) - chi2_min, p-value from
parametric-bootstrap simulations with the truth at the candidate (injection-residual noise).
Candidates: VHS residual direction (PA 345) at 0.6/0.8/1.1", the Gaia neighbours within 60" (their 2016 offsets),
and the G star itself (delta = 0)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault('V4C_NS', '10')
import numpy as np
LABEL = os.environ.get('V4C_LABEL', 'baseline')
os.environ['V4C_LABEL'] = LABEL + '_pointtest_tmp'
import importlib.util
spec = importlib.util.spec_from_file_location('comb', os.path.join(os.path.dirname(os.path.abspath(__file__)), 's05_combine.py'))
comb = importlib.util.module_from_spec(spec); spec.loader.exec_module(comb)
from v4lib import V4, gaia, TARGET, to_json

g = gaia()
it = np.where(g['source_id'] == TARGET)[0][0]
cands = [('G star', 0.0, 0.0)]
for sep, pa in ((0.6, 345), (0.8, 345), (1.1, 345), (1.1, 300), (1.5, 345)):
    cands.append((f'{sep}" PA {pa}', sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))))
for i in range(len(g['ra'])):
    if i == it:
        continue
    dE = (g['ra'][i] - g['ra'][it]) * np.cos(np.radians(g['dec'][it])) * 3600; dN = (g['dec'][i] - g['dec'][it]) * 3600
    if np.hypot(dE, dN) < 60:
        cands.append((f"Gaia {g['source_id'][i]} G={g['phot_g_mean_mag'][i]:.2f}", dE, dN))
rng = np.random.default_rng(17)
_, c2min = comb.solve('M2', comb.obs)
out = []
for name, e, n in cands:
    dt = np.array([e, n])
    ptr, c2f = comb.solve('M2', comb.obs, fix=dt)
    od = c2f - c2min
    tr = dt[None, :] + np.einsum('nip,p->ni', comb.DES['M2'][:, :, 2:], ptr)
    st = []
    nsim = 3000 if od < 60 else 300
    for _ in range(nsim):
        o = tr + comb.draw()
        st.append(comb.solve('M2', o, fix=dt)[1] - comb.solve('M2', o)[1])
    st = np.array(st)
    p = float((st >= od).mean())
    out.append(dict(name=name, dE=float(e), dN=float(n), sep=float(np.hypot(e, n)), dchi2=float(od), p_emp=p, n_sim=nsim,
                    p_upper_if_zero=float(1.0 / nsim) if p == 0 else None, sim_q99=float(np.percentile(st, 99)), levers=ptr.tolist()))
    print(f"{name:40s} sep {np.hypot(e, n):6.2f}\"  dchi2 {od:8.1f}  p_emp {p:.4f} (n={nsim}, sim q99 {np.percentile(st, 99):.1f})  levers {np.round(ptr, 2)}", flush=True)
to_json(out, f'{V4}/results/s07_point_test_{LABEL}.json')
for f in (f'{V4}/results/s05_{LABEL}_pointtest_tmp.json', f'{V4}/results/s05_{LABEL}_pointtest_tmp_map.npz'):
    if os.path.exists(f):
        os.remove(f)
