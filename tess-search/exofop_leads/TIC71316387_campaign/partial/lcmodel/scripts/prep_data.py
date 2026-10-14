import numpy as np, sys, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit2, lclib as L
d = fit2.default_data()
pickle.dump(d, open(L.HERE + '/results/data_default.pkl', 'wb'))
print(d['t'].size, len(d['evtab']), {k: (v['nev'], round(v['sigma']*1e6), round(v['beta'],2)) for k, v in d['info'].items()})
