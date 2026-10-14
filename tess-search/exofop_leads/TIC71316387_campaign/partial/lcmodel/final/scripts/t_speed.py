import os
os.environ['NUMBA_NUM_THREADS'] = '1'
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import fit3, geofit as G
d = fit3.load_data('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/results/data_default.pkl')
p = G.to_p([], [], dict(G.DEF), 3165.6697, 'total')
t0 = time.time(); fit3.chi2(d, p); print('first', time.time() - t0)
t0 = time.time()
for i in range(20): fit3.chi2(d, p)
print('per eval', (time.time() - t0) / 20)
t0 = time.time(); o = G.solve(d, G.STARTS_TOTAL[0], branch='total'); print('fit', time.time() - t0, o['tot'], o['nfev'], o['L3'], o['k'], o['b1'])
