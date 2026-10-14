"""Final O-C figure: per-sector weighted means (primary, secondary) after the best linear ephemeris, with the quadratic
fit band and example LTTE curves for a 1.03-Msun tertiary at 10, 30 and 100 au (edge-on circular, phase chosen to
maximise curvature over 2018-2026) after removing their own best linear ephemeris."""
import numpy as np, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/'
D = json.load(open(H + 'results/etvfit_prelim.json'))
E = np.array(D['E']); oc = np.array(D['oc']); es = np.array(D['es']); sec = np.array(D['sec']); ph = np.array(D['ph']); t = np.array(D['t'])
tyr = (t - 2900) / 365.25
A = np.array([np.ones_like(E), E / 1000, (ph == 0.5).astype(float)]).T
Aw = A / es[:, None]; c, *_ = np.linalg.lstsq(Aw, oc / es, rcond=None); cov = np.linalg.inv(Aw.T @ Aw)
res = oc - A @ c
print('T0 shift %.3f +- %.3f min ; P shift %.4f +- %.4f s/cycle ; secondary offset %.3f +- %.3f min' % (c[0], np.sqrt(cov[0, 0]), c[1] / 1000 * 60, np.sqrt(cov[1, 1]) / 1000 * 60, c[2], np.sqrt(cov[2, 2])))
fig, ax = plt.subplots(figsize=(9, 4.8))
yrs = lambda tt: 2000 + (tt + 2457000 - 2451544.5) / 365.25
for f, col, lab, dx in [(0.0, 'k', 'primary eclipses (sector mean)', -0.04), (0.5, 'C3', 'secondary eclipses (sector mean)', 0.04)]:
    for s in np.unique(sec):
        m = (sec == s) & (ph == f)
        w = 1 / es[m] ** 2; mu = np.sum(res[m] * w) / w.sum(); er = 1 / np.sqrt(w.sum())
        ax.errorbar(yrs(np.mean(t[m])) + dx, mu, er, fmt='o' if f == 0 else 's', color=col, ms=6, capsize=2, label=lab if s == 1 else None)
tt = np.linspace(t.min() - 30, t.max() + 30, 400); ty = (tt - 2900) / 365.25
G = 6.674e-11; MS = 1.989e30; AU = 1.496e11; C = 2.998e8
for a_au, ls in [(10, ':'), (30, '--'), (100, '-.')]:
    Mt = 1.03 + 0.64; Pout = np.sqrt(a_au ** 3 / Mt); amp = a_au * 1.03 / Mt * AU / C / 60
    best = None
    for phi in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        tau = amp * np.sin(2 * np.pi * tyr / Pout + phi)
        cc = np.polyfit(tyr, tau, 1); r = tau - np.polyval(cc, tyr)
        if best is None or np.ptp(r) > best[0]: best = (np.ptp(r), phi)
    tau = amp * np.sin(2 * np.pi * ty / Pout + best[1]); tauo = amp * np.sin(2 * np.pi * tyr / Pout + best[1])
    cc = np.polyfit(tyr, tauo, 1)
    ax.plot(yrs(tt), tau - np.polyval(cc, ty), ls, color='C0', lw=1.2, label=f'example LTTE, 1.03 Msun at {a_au} au (max-curvature phase, i=90)')
ax.axhline(0, color='0.6', lw=0.8)
ax.set_ylim(-6, 6); ax.set_xlabel('year'); ax.set_ylabel('O - C (min), linear ephemeris removed')
ax.set_title('TIC 71316387 EB (P = 1.3036721 d): eclipse timings 2018-2026')
ax.legend(fontsize=7, loc='lower left'); plt.tight_layout(); plt.savefig(H + 'plots/etv_final.png', dpi=120)
