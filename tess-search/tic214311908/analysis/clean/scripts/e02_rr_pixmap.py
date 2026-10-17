"""e02 (exploratory): where does the residual 0.6485-d signal live in the cleaned INV difference images?
Per pixel within 4 px of the M3V: 3-harmonic fit at the Gaia RRab period (plus a 1-d running-median high-pass) of the
difference image D (regression nuisance removed, Mira variation kept) after subtracting the per-cadence model-B Mira
term (b x F_Mira) and plane; first-harmonic cos/sin amplitude images are then fitted by the M3V PRF, the Mira b image
and the Mira's x/y gradient images (pointing-like dipole). Printed per mode."""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from scipy.ndimage import median_filter
s = int(sys.argv[1]) if len(sys.argv) > 1 else 93
for mode in ("none", "all"):
    S = setup(s, mode)
    R = photometry(s, S, return_D=True)
    L = make_lc(s, R, models=(("B", "M3V"),))
    D = R["D"]; F = R["FIT"]["B"]["F"]
    ok = L["quality_ok"] & ~L["intransit_window"] & np.all(np.isfinite(D.reshape(len(D), -1)), 1)
    yy, xx = np.mgrid[0:21, 0:21]
    mx, my = R["pos"]["Mira"]; m3x, m3y = R["pos"]["M3V"]
    u = (xx - mx) / 3; v = (yy - my) / 3
    model_noM3V = F[:, 0, None, None] * R["b_img"][None] + F[:, 2, None, None] + F[:, 3, None, None] * u[None] + F[:, 4, None, None] * v[None]
    Rz = (D - model_noM3V)[ok]
    t = L["t_tdb"][ok]
    reg = np.hypot(xx - m3x, yy - m3y) <= 3.0
    n = max(5, int(1 / R["exptime"]) | 1)
    ph = 2 * np.pi * (t - TMAX_RR) / P_RR
    X = np.column_stack([np.ones(len(t))] + [f(h * ph) for h in (1, 2, 3) for f in (np.cos, np.sin)])
    A = np.zeros((2,) + reg.shape); E = np.zeros(reg.shape)
    for (i, j) in zip(*np.where(reg)):
        y = Rz[:, i, j]; y = y - median_filter(y, size=n, mode="nearest")
        c, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ c
        A[0, i, j], A[1, i, j] = c[1], c[2]; E[i, j] = np.std(r) * np.sqrt(2 / len(t))
    T3 = S["prf"].image(m3x, m3y); b = R["b_img"]
    gx = np.zeros_like(b); gy = np.zeros_like(b); gx[:, 1:-1] = (b[:, 2:] - b[:, :-2]) / 2; gy[1:-1] = (b[2:] - b[:-2]) / 2
    for lab, tm in (("M3V", [T3]), ("Mira_b", [b]), ("Mira_grad", [gx, gy]), ("M3V+grad", [T3, gx, gy])):
        Xm = np.column_stack([q[reg] for q in tm] + [np.ones(reg.sum())]) / E[reg][:, None]
        chi = 0; amps = []
        for k in (0, 1):
            c, *_ = np.linalg.lstsq(Xm, A[k][reg] / E[reg], rcond=None); r = A[k][reg] / E[reg] - Xm @ c; chi += r @ r; amps.append(c[0])
        null = np.sum((A[:, reg] / E[reg]) ** 2)
        print(f"S{s} {mode}: template {lab:9s}: chi2 {chi:.1f} (null {null:.1f}, npix {reg.sum()}x2); first coef cos/sin {amps[0]:+.2f} {amps[1]:+.2f}", flush=True)
