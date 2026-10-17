import sys; sys.path.insert(0, ".")
from p05_phot import *
for s in (39, 104):
    S = setup(s); C, W, J, g, prf, rBA = S
    R = photometry(s, S=S, verbose=False)
    fm = R["fitmask"].ravel(); sig = R["sig"]
    T = R["templates"]; yy, xx = np.mgrid[0:21, 0:21]; mx, my = R["pos"]["Mira"]
    u = ((xx - mx) / 3).ravel()[fm]; v = ((yy - my) / 3).ravel()[fm]
    for lab, cols in (("A Mira,M3V,pair,plane", [T["Mira"].ravel()[fm], T["M3V"].ravel()[fm], T["pair"].ravel()[fm], np.ones(fm.sum()), u, v]),
                      ("B Mira,M3V,plane", [T["Mira"].ravel()[fm], T["M3V"].ravel()[fm], np.ones(fm.sum()), u, v]),
                      ("C Mira,M3V,pair,const", [T["Mira"].ravel()[fm], T["M3V"].ravel()[fm], T["pair"].ravel()[fm], np.ones(fm.sum())]),
                      ("D Mira,M3V,const", [T["Mira"].ravel()[fm], T["M3V"].ravel()[fm], np.ones(fm.sum())])):
        A = np.column_stack(cols) / sig[:, None]; cov = np.linalg.pinv(A.T @ A); e = np.sqrt(np.diag(cov)); cor = cov / np.outer(e, e)
        print(f"S{s} {lab}: errors {np.round(e[:3],1)}; corr M3V-Mira {cor[0,1]:+.2f}" + (f", M3V-pair {cor[1,2]:+.2f}, Mira-pair {cor[0,2]:+.2f}" if 'pair' in lab else ""))
