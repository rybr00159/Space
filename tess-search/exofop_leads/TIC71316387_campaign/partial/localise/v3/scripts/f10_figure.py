"""Summary figure: (a) raw per-sector difference-image offsets vs the direction away from the camera boresight;
(b) chromatic-corrected joint localisation (profile dchi2 map, calibrated contours) with Gaia/VHS constraints."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v3lib import V3, SECTORS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

S1, S2, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#d9d8d4'
r = json.load(open(f'{V3}/results/s05_baseline.json'))
z = np.load(f'{V3}/results/s05_baseline_map.npz')
fig, ax = plt.subplots(1, 2, figsize=(13, 6.2))
a = ax[0]
for i, lab in enumerate(r['labels']):
    o = np.array(r['obs'][i]); u = np.array(r['u_to_boresight'][i]); c = np.array(r['cov'][i])
    col = S1 if 'deep' in lab else S2
    a.errorbar(o[0], o[1], xerr=np.sqrt(c[0, 0]), yerr=np.sqrt(c[1, 1]), fmt='o', ms=8, color=col, mec='white', mew=1.5, elinewidth=1.2, capsize=0)
    if 'deep' in lab:
        a.annotate('', xy=(-3.0 * u[0], -3.0 * u[1]), xytext=(0, 0), arrowprops=dict(arrowstyle='->', color=MUTED, lw=1.2))
        a.text(-3.35 * u[0], -3.35 * u[1], lab.split('_')[0], color=MUTED, fontsize=9, ha='center', va='center')
    a.text(o[0] - 0.18, o[1] + 0.18, lab.split('_')[0], color=INK, fontsize=8)
a.plot(0, 0, '*', ms=16, color=INK)
a.plot([], [], 'o', color=S1, label='deep eclipse (bias-corrected)'); a.plot([], [], 'o', color=S2, label='shallow eclipse')
a.plot([], [], '-', color=MUTED, label='direction away from camera boresight')
a.set_xlim(6.5, -6.5); a.set_ylim(-6.5, 6.5); a.set_aspect('equal'); a.grid(color=GRID, lw=0.6)
a.set_xlabel('dE = dRA cos(dec) [arcsec]  (East left)'); a.set_ylabel('dN [arcsec]')
a.set_title('(a) Raw TESS PRF offsets of the eclipsing light from the G star', fontsize=11, loc='left')
a.legend(loc='lower left', fontsize=8, frameon=False)
b = ax[1]
g = z['g']; M = z['M']
lv = [r['thresholds'][k] for k in ('68.27', '95.45', '99.73')]
b.contourf(g, g, M, levels=[0] + lv, colors=['#2a78d6', '#7fb0ea', '#cfe0f7'], alpha=0.9)
b.contour(g, g, M, levels=lv, colors=[S1], linewidths=1.0)
b.plot(0, 0, '*', ms=16, color=INK, label='G star (Gaia DR3)')
d = r['M2_delta']; b.plot(d[0], d[1], 'x', ms=10, mew=2, color=INK, label='TESS best fit (chromatic model)')
t = np.linspace(0, 2 * np.pi, 300)
b.plot(0.6 * np.sin(t), 0.6 * np.cos(t), '--', color=S2, lw=1.5, label='VHS Ks: dKs<=2.5 excluded outside 0.6"')
b.plot(1.0 * np.sin(t), 1.0 * np.cos(t), ':', color=S2, lw=1.5, label='VHS Ks: dKs<=3.0 excluded outside 1.0"')
pa = np.radians(348); b.plot(1.13 * np.sin(pa), 1.13 * np.cos(pa), 'D', ms=7, color=S2, mec='white', label='VHS residual feature (dKs~4.2, too faint for EB)')
b.set_xlim(2.0, -2.0); b.set_ylim(-2.0, 2.0); b.set_aspect('equal'); b.grid(color=GRID, lw=0.6)
b.set_xlabel('dE [arcsec]  (East left)'); b.set_ylabel('dN [arcsec]')
b.set_title('(b) Joint 6-sector position, 68/95/99.7% (injection-calibrated)', fontsize=11, loc='left')
b.legend(loc='lower left', fontsize=8, frameon=False)
plt.tight_layout(); plt.savefig(f'{V3}/plots/fig_localisation_summary.png', dpi=110)
print('saved')
