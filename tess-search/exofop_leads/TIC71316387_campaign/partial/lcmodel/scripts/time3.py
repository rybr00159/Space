import os; os.environ['NUMBA_NUM_THREADS']='1'
import numpy as np, sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L
d = fit3.load_data(L.HERE + '/results/data_default.pkl')
p0 = dict(T0=L.T00, P=L.P0, rsum=0.2, k=0.7, bn=0.6, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3, s95=1.0, s105=1.0, s106=1.0)
M = fit3.Model(d, list(p0), {})
x0 = np.array(list(p0.values())); M(x0)
t0 = time.time()
for i in range(30): v = M(x0)
print(v, (time.time() - t0) / 30)
