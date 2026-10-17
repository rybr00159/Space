import sys; sys.path.insert(0, ".")
from p05_phot import *
from wotan import flatten
def tmpl(t, P=P_EPH, T0=T0_EPH):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5; a = np.abs(ph * P); T14 = 0.067; T23 = 0.78 * T14
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
def dip(R, mk, nm, P=P_EPH, T0=T0_EPH):
    t = R["t"]; ok = R["qok"] & (R["bkg"] < R["bkg_lim"]) & R["jit_ok"]
    names = R["FIT"][mk]["names"]; f = 1000 + R["FIT"][mk]["F"][:, names.index(nm)]
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    m = ok & np.isfinite(f) & ~(np.abs(ph * P) < 0.1) & ~(np.abs((((t - T0_EPH) / P_EPH + 0.5) % 1 - 0.5) * P_EPH) < 0.1)
    tr = flatten(t[m], f[m], method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
    kk = np.isfinite(tr); r = f - np.interp(t, t[m][kk], tr[kk])
    y = tmpl(t, P, T0); sel = ok & np.isfinite(r)
    A = np.column_stack([np.ones(sel.sum()), -y[sel]]); c, *_ = np.linalg.lstsq(A, r[sel], rcond=None)
    s_ = robust_std(r[m]); e = s_ / np.sqrt(np.sum((y[sel] - y[sel].mean()) ** 2))
    return c[1], e, s_
for s in (39, 104):
    S = setup(s)
    for rf in (1.5, 2.0, 2.5, 3.0, 4.0):
        R = photometry(s, S=S, rfit=rf, verbose=False)
        out = [f"{mk}:{nm} {dip(R, mk, nm)[0]:.1f}+-{dip(R, mk, nm)[1]:.1f} (rstd {dip(R, mk, nm)[2]:.1f})" for mk, nm in (("B", "M3V"), ("A", "M3V"), ("A", "pair"), ("A", "Mira"))]
        print(f"S{s} rfit {rf}: " + "; ".join(out), flush=True)
