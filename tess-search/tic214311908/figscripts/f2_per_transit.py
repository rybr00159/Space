"""f2: per_transit.png - (a) per-transit O-C against the adopted ephemeris for both cleaned pipelines
(referee3_clean q05: times re-measured on the delivered clean/ curves, errors from light-curve-level injections, one
value per sector and pipeline); (b) per-sector depths of the PRF photometries (INV, REF) and of the raw-pixel
apertures (3x3 box, r <= 2.5 px), all measured the same way by referee3_clean q14 (fixed shape, linear baseline,
fake-epoch errors), with the adopted depth band."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd
from astropy.time import Time

style()
SECS = (12, 39, 66, 93, 103, 104)


def main():
    q = pd.read_csv(f'{B}/referee3_clean/out/q05_per_transit.csv')
    q['oc'] = (q.tc - (T0_ADOPT + q.epoch * P_ADOPT)) * 1440
    fig = plt.figure(figsize=(10.6, 7.0))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.15], hspace=0.42, wspace=0.08)
    axs = [fig.add_subplot(gs[0, i]) for i in range(6)]
    for i, (ax, s) in enumerate(zip(axs, SECS)):
        g = q[q.sector == s]
        e0, e1 = g.epoch.min(), g.epoch.max()
        ee = np.linspace(e0 - 0.6, e1 + 0.6, 50)
        sig = np.sqrt(T0_ERR ** 2 + (ee * P_ERR) ** 2 + 2 * ee * RHO_T0P * T0_ERR * P_ERR) * 1440
        ax.fill_between(ee, -sig, sig, color=GRID, lw=0, zorder=0, label='adopted ephemeris $\\pm1\\sigma$' if i == 0 else None)
        ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
        clean_line = ((T0_CLEAN + ee * P_CLEAN) - (T0_ADOPT + ee * P_ADOPT)) * 1440
        ax.plot(ee, clean_line, '--', color=INK2, lw=0.9, zorder=1, label='clean/ combined ephemeris' if i == 0 else None)
        for pl, col, mk, dx in (('inv', C_INV, 'o', -0.13), ('ref', C_REF, 's', 0.13)):
            h = g[g.pl == pl]
            ax.errorbar(h.epoch + dx, h.oc, h.err_min, fmt=mk, ms=4.5, color=col, mec='white', mew=0.5, elinewidth=1.1, zorder=3,
                        label=('INV' if pl == 'inv' else 'REF') + ' (cleaned)' if i == 0 else None)
        ax.set_xlim(e0 - 0.7, e1 + 0.7)
        ax.set_xticks(sorted(g.epoch.unique()))
        ax.tick_params(axis='x', labelsize=7.0)
        ax.set_ylim(-48, 48)
        yr = Time(2457000 + g.tc.mean(), format='jd').to_value('decimalyear')
        ax.set_title(f'S{s} ({yr:.1f})', fontsize=9, color=INK)
        if i: plt.setp(ax.get_yticklabels(), visible=False)
        else: ax.set_ylabel('O - C (min)')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower right', bbox_to_anchor=(0.905, 0.905), ncol=4, fontsize=7.8, handlelength=1.8, columnspacing=1.4)
    fig.text(0.5, 0.515, 'Epoch (transit number; epoch 0 = BTJD 3840.915, 2025 June)', ha='center', fontsize=9, color=INK)
    axs[0].text(-0.02, 1.25, 'a', transform=axs[0].transAxes, fontsize=11, fontweight='bold', ha='right', va='bottom')
    axs[0].text(0.0, 1.25, 'Transit times minus the adopted ephemeris (no TTVs)', transform=axs[0].transAxes, fontsize=10, color=INK,
                ha='left', va='bottom')
    # ---------------- (b) per-sector depths
    ax = fig.add_subplot(gs[1, :])
    c14 = json.load(open(f'{B}/referee3_clean/out/q14_contrast_all.json'))
    sets = [('INV_delivered', 'INV PRF photometry', C_INV, 'o', 'full', -0.27),
            ('REF_delivered', 'REF PRF photometry', C_REF, 's', 'full', -0.09),
            ('AP_box3', 'raw pixels, 3x3 box', AQUA_DARK, '^', 'full', 0.09),
            ('AP_m25', 'raw pixels, r $\\leq$ 2.5 px', AQUA_DARK, 'v', 'none', 0.27)]
    D, Ds, Dh = 15.8, 0.9, (1.1, 1.5)
    lo, hi = D - np.hypot(Ds, Dh[0]), D + np.hypot(Ds, Dh[1])
    ax.axhspan(lo, hi, color='#eeeeea', lw=0, zorder=0)
    ax.axhspan(D - Ds, D + Ds, color='#dcdcd5', lw=0, zorder=0)
    ax.axhline(D, color=INK2, lw=1.0, zorder=1)
    xs = np.arange(len(SECS))
    for key, lab, col, mk, fill, dx in sets:
        ps = c14[key]['per_sector']
        y = np.array([ps[str(s)][0] for s in SECS]); e = np.array([ps[str(s)][1] for s in SECS])
        cst = c14[key]['contrast']; ch = c14[key]['chi2']
        ax.errorbar(xs + dx, y, e, fmt=mk, ms=6, color=col, mfc=col if fill == 'full' else 'white', mec=col if fill != 'full' else 'white',
                    mew=1.0 if fill != 'full' else 0.6, elinewidth=1.2, zorder=3,
                    label=f'{lab}: late$-$early {cst[0]:+.1f} $\\pm$ {cst[1]:.1f} ppt, sector $\\chi^2$ {ch[0]:.1f}/{ch[1]}')
    pn = c14['REF_delivered_noE58']['per_sector']['103']
    ax.errorbar(xs[4] - 0.09, pn[0], pn[1], fmt='s', ms=6, mfc='white', mec=C_REF, color=C_REF, mew=1.0, elinewidth=0.9, alpha=0.9, zorder=2)
    ax.annotate('REF S103 without the\nscattered-light event E58', xy=(xs[4] - 0.09, pn[0] - 1.5), xytext=(xs[4] - 0.95, 6.2), fontsize=7.4,
                color=INK2, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    ax.set_xticks(xs)
    ax.set_xticklabels([f'S{s}\n' + Time(2457000 + q[q.sector == s].tc.mean(), format='jd').to_value('iso')[:7] for s in SECS])
    ax.set_xlim(-0.6, len(SECS) - 0.4)
    ax.set_ylim(4, 32)
    ax.set_ylabel('Transit depth per sector (ppt of M3V flux)')
    ax.text(5.55, 30.8, 'Line and dark band: adopted depth 15.8 $\\pm$ 0.9 (stat) ppt\nlight band: including the systematic +1.5/$-$1.1 ppt',
            ha='right', va='top', fontsize=7.8, color=INK2)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fontsize=7.8, ncol=2, columnspacing=1.6, handletextpad=0.4)
    ax.set_title('Depth by sector: the late/early difference of the PRF photometry is absent in the raw pixels', loc='left', fontsize=10)
    ax.text(-0.02, 1.02, 'b', transform=ax.transAxes, fontsize=11, fontweight='bold', ha='right', va='bottom')
    ax.grid(axis='x', visible=False)
    save(fig, 'per_transit.png')
    for key in c14:
        print(key, {s: [round(v, 1) for v in c14[key]['per_sector'][str(s)]] for s in SECS}, c14[key]['contrast'][:2], c14[key]['chi2'])
    print('O-C rms (min):', q.groupby('pl').oc.std().round(1).to_dict(), 'max |O-C|/err', (q.oc.abs() / q.err_min).max().round(2))


if __name__ == '__main__':
    main()
