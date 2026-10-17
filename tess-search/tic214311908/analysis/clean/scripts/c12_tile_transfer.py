"""c12: transfer function of the tiled REF light curve for a signal at the RRab period: a sinusoid of 3.0 e-/s
semi-amplitude at P_RR (phase-locked to the Gaia maximum) is injected at the M3V (REF PRF) into the cube, the tiled
REF-method curve is rebuilt (cleaned aux) and the first-harmonic amplitude is re-measured (minus the no-injection value,
vector difference). Gives the factor by which tiled-REF RR amplitudes and upper limits must be divided.
Output: out/tile_transfer.json"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from refclean import *
from cleanlib import P_RR, TMAX_RR, F_M3V, harm_fit, P1, T01
out = {}
for s in (int(a) for a in (sys.argv[1:] or ["39", "93", "103"])):
    S = SectorC(s, "all"); flat0 = S.flat.copy()
    res = []
    for inj in (False, True):
        if inj:
            amp = 3.0 * np.cos(2 * np.pi * (S.ttdb - TMAX_RR) / P_RR); amp[~np.isfinite(amp)] = 0
            S.flat = flat0 + amp[:, None] * S.P_M3V.ravel()[None]
        t, a, sg, w = tiled(S)
        k = np.abs((((t - T01) / P1 + 0.5) % 1 - 0.5) * P1) > 0.1
        _, _, c = harm_fit(t[k], a[k], 1 / sg[k] ** 2, P_RR, TMAX_RR, 3)
        res.append(c[1:3])
    S.flat = flat0
    d = res[1] - res[0]
    out[s] = dict(recovered_cos=float(d[0]), recovered_sin=float(d[1]), transfer=float(np.hypot(*d) / 3.0))
    print(s, out[s], flush=True)
jd(out, f"{CL}/out/tile_transfer.json")
