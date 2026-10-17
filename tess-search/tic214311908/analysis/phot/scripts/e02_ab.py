import sys; sys.path.insert(0, ".")
from p05_phot import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
s = int(sys.argv[1])
R = photometry(s, verbose=False)
a, b, be = R["a_img"], R["b_img"], R["berr"]
mx, my = R["pos"]["Mira"]
ix, iy = int(round(mx)), int(round(my))
np.set_printoptions(linewidth=250, precision=0, suppress=True)
sl = (slice(iy - 3, iy + 4), slice(ix - 3, ix + 4))
print("mean image a (bkg-subtracted by median of faint px)\n", (a - np.percentile(a, 15))[sl][::-1])
print("b image x 1e4\n", (b * 1e4)[sl][::-1])
print("b err x 1e4\n", (be * 1e4)[sl][::-1])
print("a/b (only where b>0.01)\n", np.where(b > 0.01, (a - np.percentile(a, 15)) / b, np.nan)[sl][::-1])
print("Mira pos", mx, my, "M3V", R["pos"]["M3V"], "pair", R["pos"]["pair"], "sum b in ap", b[np.hypot(*np.mgrid[0:21,0:21][::-1] - np.array([mx,my])[:,None,None]) <= 2.5].sum())
