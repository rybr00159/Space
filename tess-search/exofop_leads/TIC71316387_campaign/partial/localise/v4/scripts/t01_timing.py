import sys, time
sys.path.insert(0, '.')
from v4lib import *
S = Sector(28)
t0 = time.time(); G = PRFGrid(S); print('grid', time.time() - t0, G.M.shape, len(G.nb))
z = np.load(f'{V4}/results/win_s28.npz', allow_pickle=True)
A = z['A'][0]; sA = z['sA'][0]
t0 = time.time()
for i in range(50):
    G.fit(A, sA)
print('fit ms', (time.time() - t0) / 50 * 1000)
t0 = time.time()
for i in range(50):
    G.img(S.xt + 0.1, S.yt)
print('img ms', (time.time() - t0) / 50 * 1000)
