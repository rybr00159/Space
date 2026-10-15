"""Injection data sets: a synthetic EB with known parameters (truths) multiplied into the REAL light curves of the
default product set, with its primary conjunction placed at orbital phase 0.25 of the real EB (so its eclipses fall
at phases 0.25 and 0.75, where the real light curve is eclipse-free), then windowed, clipped and noise-estimated by
exactly the same code as the real data.  The injected model includes exposure smearing and the per-product depth
scales of the truth.  Usage: c05_inject_build.py <truth json> <tag>"""
import numpy as np, sys, json, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import data as D, ebm
V = D.V
truth = json.load(open(sys.argv[1])); tag = sys.argv[2]
p = truth['p']
T0f = p['T0'] + 0.25 * p['P']
gm = ebm.GridModel(p['P'])
lcs = [D.load_spoc(1), D.load_spoc(28), D.load_spoc(68), D.load_tc(95), D.load_tc(105), D.load_tc(106)]
for lc in lcs:
    key = f"{lc['prod']}{lc['sector']}"
    tph = ((lc['t'] - T0f + 0.25 * p['P']) % p['P']) - 0.25 * p['P']
    f = gm.eval(tph, np.full(lc['t'].size, lc['texp']), p['P'], p['rsum'], p['k'], p['b0'], p['J'], p['ecosw'], p['esinw'],
                p['u1a'], p['usa'] - p['u1a'], p['u1b'], p['usb'] - p['u1b'])
    dep = p.get(f's_{key}', 1.0) * (1 - p['L3']) * (1 - f)
    lc['f'] = lc['f'] * (1 - dep)
d = D.build(lcs, P=p['P'], T0=T0f)
d['truth'] = dict(p, T0=T0f)
pickle.dump(d, open(f'{V}/results/data_inj{tag}.pkl', 'wb'))
json.dump(dict(p=dict(p, T0=T0f)), open(f'{V}/results/truth_inj{tag}.json', 'w'), indent=1)
print('built', tag, d['t'].size, 'points', len(d['evtab']), 'events')
