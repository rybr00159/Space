import numpy as np, sys, pickle, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit2, lclib as L
d = pickle.load(open(L.HERE + '/results/data_default.pkl', 'rb'))
p0 = dict(T0=L.T00, P=L.P0, rsum=0.2, k=0.7, cosi=0.12, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3, s95=1.0, s105=1.0, s106=1.0)
free = [k for k in p0]
post = fit2.Post(d, free, {})
x0 = np.array([p0[n] for n in free])
post(x0)
t0=time.time()
for i in range(50): v=post(x0)
print(v, (time.time()-t0)/50)
