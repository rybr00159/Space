"""Build eclipse-window data sets (+-0.20 d) for the fits.
  default : SPOC PDCSAP S1/S28/S68 + own TESScut photometry S95/S105/S106
  tcall   : own TESScut photometry for all six sectors (check of SPOC vs own photometry)
  qlp     : SPOC S1/S28/S68 + QLP SAP S95/S105 + TESScut S106 (QLP S105 times are not reliable -> depth check only)"""
import numpy as np, sys, pickle, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import data as D
OUT = D.V + '/results'
sets = {
    'default': lambda: [D.load_spoc(1), D.load_spoc(28), D.load_spoc(68), D.load_tc(95), D.load_tc(105), D.load_tc(106)],
    'tcall': lambda: [D.load_tc(s) for s in (1, 28, 68, 95, 105, 106)],
    'qlp': lambda: [D.load_spoc(1), D.load_spoc(28), D.load_spoc(68), D.load_qlp(95), D.load_qlp(105), D.load_tc(106)],
}
summ = {}
for name, fn in sets.items():
    print(name, flush=True)
    d = D.build(fn())
    pickle.dump(d, open(f'{OUT}/data_{name}.pkl', 'wb'))
    summ[name] = dict(npts=int(d['t'].size), nev=len(d['evtab']), info=d['info'])
json.dump(summ, open(f'{OUT}/data_summary.json', 'w'), indent=1)
lc = D.load_spoc(1); print('SPOC CROWDSAP/FLFRCSAP', {s: (D.load_spoc(s)['crowd'], D.load_spoc(s)['flfr']) for s in (1, 28, 68)})
