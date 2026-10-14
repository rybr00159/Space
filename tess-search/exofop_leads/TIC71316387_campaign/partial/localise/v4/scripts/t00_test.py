import sys, time
sys.path.insert(0, '.')
from v4lib import *
t0 = time.time()
S = Sector(95)
print('load', time.time() - t0, 'cam', S.cam, S.ccd, 'xt,yt', S.xt, S.yt, 'cut', S.ncut, 'texp', S.texp * 1440, 'J', S.J)
# PRF centroid slope check
G = PRFGrid(S)
print('grid', time.time() - t0, 'nb', G.nb)
xs = S.xt + np.linspace(-0.5, 0.5, 11)
yy, xx = np.mgrid[0:S.ny, 0:S.nx]
cen = []
for x in xs:
    m = G.img(x, S.yt)
    cen.append((m * xx).sum() / m.sum())
print('PRF moment centroid vs input: slope', np.polyfit(xs, cen, 1), 'resid', np.round(np.array(cen) - np.polyval(np.polyfit(xs, cen, 1), xs), 4))
t1 = time.time()
W = S.all_windows()
print('windows', time.time() - t1, {k: int((W['kind'] == k).sum()) for k in ('deep', 'shallow', 'fake')})
for kind in ('deep', 'shallow'):
    idx = np.where(W['kind'] == kind)[0]
    A, sA, B = stack(W, idx)
    rB = G.fit(B, sA * 0 + np.sqrt(np.abs(B) + 100.), direct=True)
    rA = G.fit(A, sA)
    rc = G.fit_cont(A, sA, rA['x'], rA['y'])
    off = S.J @ np.array([rA['x'] - rB['x'], rA['y'] - rB['y']])
    print(kind, len(idx), 'B', rB, '\n  A', rA, '\n  cont', rc, '\n  off', off, 'depth_pk_ppm', rA['a'] / rB['a'] * 1e6)
