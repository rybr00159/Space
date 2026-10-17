"""f4: host.png - is the host a single metal-rich M3V or a near-equal twin?
(a) absolute M_Ks vs BP-RP of Gaia-clean nearby M dwarfs (referee_star/data/mdw_cmd.csv: plx > 12 mas, plx/err > 20,
BP-RP 2.55-2.95, 2MASS) with the KDE-mode ridge of referee2_fpp r02 and the target;
(b) the G-Ks colour residual at fixed BP-RP (r) against the M_Ks offset from the ridge (x) for the 2,126 clean stars
within 0.1 mag of the target's BP-RP (the r02 sample): running median, r02's fitted single-star (metallicity) trend,
and the track of unresolved binaries built from r02's saved Mann/Mamajek shifts (normal primary + companion of mass
ratio q). The ridge, residual definitions and sample cuts are those of referee2_fpp/scripts/r02_cmd_colour.py,
recomputed here; the fitted trend and binary shifts are read from referee2_fpp/out/r02_cmd_colour.json."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd
from scipy import stats

style()
TBPRP, TG, TK, MKS_OBS = 2.756149, 14.281043, 10.482, 6.4406     # referee2_fpp rlib / r02


def main():
    R = json.load(open(f'{B}/referee2_fpp/out/r02_cmd_colour.json'))
    d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
    d = d[np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A') & (d.ks_msigcom < 0.05) & np.isfinite(d.ruwe)].copy()
    d['MKs'] = d.ks_m + 5 * np.log10(d.parallax / 100.)
    d['GK'] = d.phot_g_mean_mag - d.ks_m
    cb = np.arange(2.55, 2.951, 0.05); modes, cens = [], []
    for lo in cb[:-1]:
        s = d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05) & (d.ruwe < 1.1)].MKs.values
        xs = np.linspace(np.percentile(s, 1), np.percentile(s, 99), 3000)
        modes.append(xs[np.argmax(stats.gaussian_kde(s, 0.12)(xs))]); cens.append(np.median(d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05)].bp_rp))
    SL, IC = np.polyfit(cens, modes, 1)
    d['x'] = d.MKs - (IC + SL * d.bp_rp)
    core = d[(d.ruwe < 1.1) & (np.abs(d.x) < 0.25)]
    CG = np.polyfit(core.bp_rp, core.GK, 1)
    d['r'] = d.GK - np.polyval(CG, d.bp_rp)
    x_t = MKS_OBS - (IC + SL * TBPRP); r_t = (TG - TK) - np.polyval(CG, TBPRP)
    print(f'ridge slope {SL:.4f} (r02 {R["ridge"]["slope"]:.4f}), intercept {IC:.4f} (r02 {R["ridge"]["intercept"]:.4f}); '
          f'target x {x_t:.3f} (r02 {R["target"]["x"]:.3f}), r {r_t:.3f} (r02 {R["target"]["r"]:.3f})')
    clean = (d.ruwe < 1.15) & (d.ipd_frac_multi_peak == 0) & (d.ipd_gof_harmonic_amplitude < 0.04) & (d.non_single_star == 0)
    win = clean & (np.abs(d.bp_rp - 2.756) < 0.10) & (d.x > -1.6) & (d.x < 1.0)
    S = d[win]
    print('clean stars', int(clean.sum()), 'window sample', len(S), '(r02 n', R['n'], ')')
    frac_bright = float(np.mean(S.x <= x_t))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.6, 4.9), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.28))
    # ---------------- (a) CMD
    C = d[clean]
    ax.scatter(C.bp_rp, C.MKs, s=2.0, color=MUTED, alpha=0.35, lw=0, rasterized=True, label=f'Gaia-clean M dwarfs within 83 pc (n = {len(C):,})')
    bb = np.linspace(2.55, 2.95, 10)
    ax.plot(bb, IC + SL * bb, '-', color=BLUE, lw=1.8, label='single-star ridge (density peak)')
    ax.plot(bb, IC + SL * bb - 0.753, '--', color=ORANGE, lw=1.6, label='equal-mass twin of a ridge star (0.753 mag brighter)')
    ax.axvspan(2.656, 2.856, color=GRID, alpha=0.5, lw=0, zorder=0, label='colour window used in panel b')
    ax.plot(TBPRP, MKS_OBS, marker='*', ms=17, color=INK, mec='white', mew=1.0, zorder=5, ls='none', label='TIC 214311908 (M3V)')
    ax.annotate(f'$\\Delta M_{{K_s}}$ = {x_t:+.2f} mag'.replace('-', '\u2212') + f'\n(local ridge of referee2_host: $-$0.75)', xy=(TBPRP, MKS_OBS), xytext=(2.79, 5.55),
                fontsize=8, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.7))
    ax.set_xlim(2.55, 2.95); ax.set_ylim(9.2, 5.2)
    ax.set_xlabel('Gaia BP $-$ RP (mag)'); ax.set_ylabel('Absolute $M_{K_s}$ (mag)')
    ax.legend(loc='lower left', fontsize=7.4, markerscale=1.0, handlelength=1.8)
    ax.set_title('$\\bf{a}$   Over-luminous for its colour', loc='left')
    # ---------------- (b) colour residual vs offset
    bx.scatter(S.x, S.r, s=3.0, color=MUTED, alpha=0.4, lw=0, rasterized=True, label=f'clean stars with |BP$-$RP $-$ 2.756| < 0.1 (n = {len(S):,})')
    eds = np.array([-1.5, -1.1, -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0])
    xm, ym, ye = [], [], []
    for a_, b_ in zip(eds[:-1], eds[1:]):
        g = S[(S.x >= a_) & (S.x < b_)]
        if len(g) < 12: continue
        xm.append(g.x.median()); ym.append(g.r.median()); ye.append(1.2533 * g.r.std() / np.sqrt(len(g)))
    bx.errorbar(xm, ym, ye, fmt='o', ms=4.5, color=BLUE, mec='white', mew=0.5, elinewidth=1.1, zorder=4, label='running median (single stars dominate)')
    prow = [r for r in R['profiles']['F_twin=0.1']['rows'] if abs(r['sig_b'] - 0.40) < 1e-6][0]['params']
    xx = np.linspace(-1.5, 1.0, 100)
    tr = prow['a'] + prow['c'] * xx + prow['c2'] * xx ** 2
    bx.fill_between(xx, tr - prow['sr'], tr + prow['sr'], color=BLUE, alpha=0.12, lw=0, zorder=1)
    bx.plot(xx, tr, '-', color=BLUE, lw=1.6, zorder=3, label=f'fitted single-star trend $\\pm$ scatter ({prow["sr"]:.3f})')
    q = np.array(R['shifts']['q']); dx = np.array(R['shifts']['dx']); dr = np.array(R['shifts']['dr'])
    bx.plot(dx, prow['a'] + dr, '--', color=ORANGE, lw=1.7, zorder=3, label='unresolved binary: normal primary + companion')
    for qq, ddx, ddy, ha in ((0.6, 0.09, -0.037, 'left'), (0.8, 0.04, -0.006, 'left'), (0.9, 0.03, -0.022, 'left'), (1.0, -0.035, 0.0, 'right')):
        j = int(np.argmin(np.abs(q - qq)))
        bx.plot(dx[j], prow['a'] + dr[j], 'o', ms=4.5, color=ORANGE, mec='white', mew=0.5, zorder=4)
        bx.text(dx[j] + ddx, prow['a'] + dr[j] + ddy, f'q = {qq:g}', fontsize=7.4, color='#b2471b', ha=ha, va='center')
        print('q', qq, 'dx', round(dx[j], 3), 'r', round(prow['a'] + dr[j], 3))
    bx.plot(x_t, r_t, marker='*', ms=17, color=INK, mec='white', mew=1.0, zorder=6, ls='none', label='TIC 214311908')
    rs = prow['a'] + prow['c'] * x_t + prow['c2'] * x_t ** 2
    bx.annotate(f'TIC 214311908: r = {r_t:+.3f}\nsingle-star trend here: {rs:+.3f} $\\pm$ {prow["sr"]:.3f}\nnear-equal twin: about 0.00',
                xy=(x_t, r_t), xytext=(-0.55, 0.205), fontsize=7.8, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.7))
    bx.axvline(0, color=BASE, lw=0.8, zorder=0); bx.axhline(0, color=BASE, lw=0.8, zorder=0)
    bx.set_xlim(-1.5, 1.0); bx.set_ylim(-0.2, 0.26)
    bx.set_xlabel('$M_{K_s}$ offset from the ridge, x (mag; negative = brighter)')
    bx.set_ylabel('$G - K_s$ residual at fixed BP$-$RP, r (mag)')
    bx.legend(loc='lower left', fontsize=7.2, handlelength=1.8, ncol=1)
    bx.set_title('$\\bf{b}$   Red in G$-K_s$ like metal-rich single stars', loc='left')
    save(fig, 'host.png')
    print('fraction of window stars at least as bright as target', round(frac_bright, 4), 'single trend at target', round(rs, 3))


if __name__ == '__main__':
    main()
