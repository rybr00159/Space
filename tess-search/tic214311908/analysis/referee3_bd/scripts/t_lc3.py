"""t_lc3: validation of lc3 binned templates against batman's own supersampling (midpoint-offset, own construction)."""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from lc3 import *
import time
for (k, b, rho) in [(0.12, 0.5, 9.4), (0.2, 0.97, 1.8), (0.4, 0.2, 20.)]:
    aR = float(aR_of_rho(rho)); cur = binned(k, b, aR)
    for c in CADS:
        t = np.linspace(-0.15, 0.15, 3001); n = 41
        ts = (t[:, None] + ((np.arange(n) + 0.5)/n - 0.5)[None, :]*c/DAY).ravel()
        p = batman.TransitParams(); p.t0 = 0; p.per = PC; p.rp = k; p.a = aR; p.inc = float(np.degrees(np.arccos(b/aR))); p.ecc = 0; p.w = 90; p.u = [U1, U2]; p.limb_dark = 'quadratic'
        ref = batman.TransitModel(p, ts).light_curve(p).reshape(-1, n).mean(axis=1) - 1
        mine = np.interp(t, *cur[c])
        print(k, b, rho, c, 'max |diff| ppm', round(np.max(np.abs(mine - ref))*1e6, 2), 'depth ppt', round(-ref.min()*1e3, 2))
r = load_ref(); i = load_inv()
print('REF N', len(r.t), 'nw', r.nw, ' INV N', len(i.t), 'nw', i.nw)
t0 = time.time(); bank = TemplateBank(np.r_[np.linspace(0.09, 0.15, 7), [0.19, 0.20, 0.22, 0.25, 0.28, 0.32]], nb=30); print('bank', len(bank.items), time.time() - t0)
t0 = time.time(); sc = bank.scan(r); print('scan', time.time() - t0, min(sc)[:4])
