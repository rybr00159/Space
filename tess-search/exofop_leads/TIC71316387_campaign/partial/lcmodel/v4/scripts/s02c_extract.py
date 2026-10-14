"""Extract from ml_default_recent.json: best_geo.json (global ML solution, used as ETV template and depth reference) and the
injection truths (truths_A.json): the global best (total), the best total solution at L3 = 0.85, the best grazing solution with
L3 = 0 (the 'G star eclipsed' alternative), and the best partial solution at L3 = 0.90."""
import json, numpy as np
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results'
res = json.load(open(f'{R}/ml_default_recent.json'))
def pick(kind, branch, key=None, val=None):
    rr = [r for r in res if r['kind'] == kind and r['branch'] == branch and (key is None or abs(r['fixed'][key] - val) < 1e-9)]
    return min(rr, key=lambda r: r['tot'])
best = min([r for r in res if r['kind'] == 'branch'], key=lambda r: r['tot'])
keep = ['T0', 'P', 'rsum', 'k', 'cosi', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b', 's1', 's28', 's68', 's95', 's105', 's106',
        'r1', 'r2', 'b1', 'inc', 'total', 'fEB', 'tot', 'chi2']
json.dump({k: best[k] for k in keep}, open(f'{R}/best_geo.json', 'w'), indent=1)
tk = ['rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b']
truths = {'flat_best': {k: best[k] for k in tk},
          'total_L3_085': {k: pick('L3prof', 'total', 'L3', 0.85)[k] for k in tk},
          'partial_L3_0': {k: pick('branch', 'partial')[k] for k in tk},
          'partial_L3_090': {k: pick('L3prof', 'partial', 'L3', 0.90)[k] for k in tk}}
json.dump(truths, open(f'{R}/truths_A.json', 'w'), indent=1)
for n, t in truths.items():
    print(n, {k: round(v, 4) for k, v in t.items()})
print('best', {k: (round(best[k], 6) if isinstance(best[k], float) else best[k]) for k in keep})
