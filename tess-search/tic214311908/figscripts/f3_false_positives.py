"""f3: false_positives.png - (a) sign-fixed scenario posteriors for the REF and INV light curves
(referee3_bd/out/s02_fpp_fixed.csv, adopted setting: photo-eccentric sign fixed, Kipping planet e prior, short-period
BD e law, f_BD = 1.3e-4); brown-dwarf channels shown both as computed by the sign-fixed grid (thin outline bar) and
after referee3_bd's two calibrations (filled bar): the Gaia RVS factor 0.49 (median, applied to BDs on the M3V and
on the primary of a binary) and the injection-calibrated shape likelihood ratio (x0.076 REF, x0.29 INV, all BD
channels) - the same recipe as referee3_bd/out/s12_summary.json, split by channel;
(b) the follow-up FPP ladder recomputed with the adopted composition by writeup/scripts/w05_fpp_ladder.py
(writeup/fpp_ladder_adopted.json; ranges over the three brown-dwarf calibrations, the Zorro contrast curves for the
speckle-only step, the colour models and the 5/10% depth-ratio precisions), with the follow-up referee's published
ladder (referee3_followup/out/r05_hri_ladder.json; f_BD = 1e-3, uncalibrated shape grid) as open markers.
The scenario called PTP in the FPP code (a planet on an unresolved companion) is labelled with its TRICERATOPS name,
STP, because TRICERATOPS uses 'PTP' for a planet on the primary."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd

style()


def main():
    S2 = pd.read_csv(f'{B}/referee3_bd/out/s02_fpp_fixed.csv')
    S12 = json.load(open(f'{B}/referee3_bd/out/s12_summary.json'))
    rvmed = S12['rv_factor_empirical_masses_min_med_max'][1]
    rows, shapecal, fpps = {}, {}, {}
    for ds in ('ref_beta', 'inv_beta'):
        r = S2[S2.tag == f'{ds} planet-e kipping BD-e bdsp f_BD 0.00013'].iloc[0]
        rows[ds] = r
        t = [x for x in S12['table'] if x['ds'] == ds and x['planet_e'] == 'kipping' and x['bd_e'] == 'bdsp'][0]
        shapecal[ds] = t['shape_corr_factor']
        fpps[ds] = (r.FPP, t['BD_total_RVmed_shapecal'], t['BD_total_noRV'])
        # check: per-channel calibrated sum equals s12's total
        tot = (r.BD_on_M3V + r.BD_on_primary) * rvmed * t['shape_corr_factor'] + r.BD_on_companion * t['shape_corr_factor']
        print(ds, 'FPP', round(r.FPP, 5), 'BD cal total', f"{tot:.3g}", 's12', f"{t['BD_total_RVmed_shapecal']:.3g}", 'P_any_star', r.P_planet_any_star)
    cats = [('TP_single', 'Planet on the M3V (single star)', 'planet'),
            ('TP_diluted', 'Planet on the M3V, diluted by a companion', 'planet'),
            ('TP_twin', 'Planet on a companion within 0.75 mag (twin)', 'planet'),
            ('PTP', 'Planet on a fainter bound companion ("STP")', 'fp'),
            ('NTP_pairA', 'Planet / EB on the 12" pair, star A (NTP)', 'fp'),
            ('WD', 'White dwarf totally occulted by the M3V', 'fp'),
            ('BG', 'Background-star scenarios', 'fp'),
            ('HEB', 'Hierarchical eclipsing binary (bound)', 'fp'),
            ('BD_on_M3V', 'Brown dwarf grazing the M3V (eccentric)', 'bd'),
            ('BD_on_companion', 'Brown dwarf on a bound companion', 'bd'),
            ('BD_on_primary', 'Brown dwarf on the primary of a binary', 'bd')]
    fig = plt.figure(figsize=(10.6, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.32, 1], wspace=0.62)
    ax = fig.add_subplot(gs[0, 0]); axl = fig.add_subplot(gs[0, 1])
    y = np.arange(len(cats))[::-1].astype(float)
    hgt = 0.36
    for j, (ds, col, lab) in enumerate((('ref_beta', C_REF, 'REF light curve'), ('inv_beta', C_INV, 'INV light curve'))):
        r = rows[ds]
        for i, (key, name, kind) in enumerate(cats):
            v = r[key]
            yy = y[i] + (0.19 if j == 0 else -0.19)
            if kind == 'bd':
                fac = shapecal[ds] * (rvmed if key in ('BD_on_M3V', 'BD_on_primary') else 1.0)
                vc = v * fac
                ax.barh(yy, v, height=hgt, left=1e-9, color='none', edgecolor=col, lw=1.0, zorder=2)
                ax.barh(yy, vc, height=hgt, left=1e-9, color=col, edgecolor='white', lw=0.5, zorder=3)
            else:
                ax.barh(yy, v, height=hgt, left=1e-9, color=col, edgecolor='white', lw=0.5, zorder=3,
                        label=lab if i == 0 else None)
    ax.set_xscale('log')
    ax.set_xlim(1e-7, 1.6)
    ax.set_yticks(y)
    ax.set_yticklabels([c[1] for c in cats], fontsize=8.2)
    for i, (key, name, kind) in enumerate(cats):
        if kind == 'planet':
            ax.get_yticklabels()[i].set_fontweight('bold')
    ax.axhline(y[2] - 0.5, color=BASE, lw=0.8)
    ax.axhline(y[7] - 0.5, color=BASE, lw=0.8)
    ax.set_xlabel('Posterior probability (log scale)')
    ax.grid(axis='y', visible=False)
    from matplotlib.patches import Patch
    hl, ll = ax.get_legend_handles_labels()
    hl += [Patch(facecolor='none', edgecolor=INK2, lw=1.0), Patch(facecolor=INK2, edgecolor='white')]
    ll += ['brown dwarf, sign-fixed grid', 'brown dwarf, after Gaia RV and\ncalibrated shape factors']
    ax.legend(hl, ll, loc='lower right', fontsize=7.6, handlelength=1.4, borderaxespad=0.2)
    fr, fi = rows['ref_beta'].FPP, rows['inv_beta'].FPP
    ax.set_title('$\\bf{a}$   Scenario posteriors (sign-fixed Bayesian FPP)', loc='left')
    pa = {ds: rows[ds].P_planet_any_star + rows[ds].NTP_pairA for ds in rows}
    ax.text(1.25, 0.5 * (y[6] + y[7]), f'FPP = P(not a planet on the M3V\nor its twin) = {fr:.4f} (REF), {fi:.4f} (INV)\n\n'
            f'P(planet on the M3V or an unresolved\ncompanion) = {rows["ref_beta"].P_planet_any_star:.4f} (REF), '
            f'{rows["inv_beta"].P_planet_any_star:.4f} (INV);\nincluding the bound pair star A: '
            f'{pa["ref_beta"]:.4f} / {pa["inv_beta"]:.4f}', fontsize=7.8, color=INK2, va='center', ha='right')

    # ---------------- (b) ladder: adopted composition (filled) and the published referee ladder (open)
    LA = json.load(open(f'{W}/fpp_ladder_adopted.json'))
    PUB = 'published composition (f_BD 1e-3, uncalibrated shape grid)'
    names = [('now', 'Now (TESS + archival data)'),
             ('speckle', '+ 8-m speckle imaging\n(Gemini-S/Zorro 832 nm, null)'),
             ('two-band only', '+ two-band ground transit\n(5-10% depth ratio), no imaging'),
             ('speckle + two-band', 'speckle (conservative curve)\n+ two-band transit'),
             ('speckle + two-band + BD closed', 'speckle + two-band transit\n+ brown-dwarf channel closed\n(ground transit shape or RV)')]
    def steps(ds):
        return [(lab, LA['ranges_all_bd'][ds][k], LA['ranges'][ds][PUB][k]) for k, lab in names]
    yl = np.arange(5)[::-1].astype(float)
    for j, (ds, col, mk) in enumerate((('ref_beta', C_REF, 's'), ('inv_beta', C_INV, 'o'))):
        for i, (name, (lo, hi), (plo, phi)) in enumerate(steps(ds)):
            yy = yl[i] + (0.20 if j == 0 else -0.08)
            axl.plot([lo, hi], [yy, yy], '-', color=col, lw=3.2, solid_capstyle='round', alpha=0.85, zorder=2)
            axl.plot([lo, hi] if hi / lo > 1.03 else [lo], [yy, yy] if hi / lo > 1.03 else [yy], mk, color=col, ms=5.5, mec='white', mew=0.6, zorder=3,
                     label=(('REF' if j == 0 else 'INV') + ', adopted inputs') if i == 0 else None)
            axl.text(hi * 1.25, yy, f'{lo:.1e}' if hi / lo < 1.03 else f'{lo:.1e}$-${hi:.1e}', fontsize=7.0, va='center', color=INK2)
            yp = yy - 0.12
            if abs(plo / lo - 1) > 0.10 or abs(phi / hi - 1) > 0.10:      # draw the published ladder only where it differs
                axl.plot([plo, phi], [yp, yp], '-', color=col, lw=0.9, alpha=0.7, zorder=2)
                axl.plot([plo, phi] if phi / plo > 1.03 else [plo], [yp, yp] if phi / plo > 1.03 else [yp], mk, mfc='white', mec=col,
                         ms=4.2, mew=0.9, zorder=3, label=(('REF' if j == 0 else 'INV') + ', published ladder (f$_{BD}$ = 10$^{-3}$)') if i == 1 else None)
    axl.set_yticks(yl)
    axl.set_yticklabels([s_[0] for s_ in steps('ref_beta')], fontsize=8.0)
    axl.set_xscale('log')
    axl.set_xlim(1e-5, 6e-2)
    axl.axvline(0.015, color=MUTED, lw=0.9, ls=':')
    axl.text(0.0145, yl[0] + 0.48, 'FPP 0.015', rotation=0, fontsize=7.2, color=MUTED, ha='right', va='center')
    axl.set_xlabel('FPP (log scale)')
    axl.grid(axis='y', visible=False)
    axl.set_ylim(-0.6, 4.75)
    axl.legend(loc='upper center', bbox_to_anchor=(0.42, -0.11), ncol=2, fontsize=7.2, handlelength=1.2, columnspacing=1.0)
    axl.set_title('$\\bf{b}$   Expected FPP after follow-up', loc='left')

    save(fig, 'false_positives.png')
    for ds in ('ref_beta', 'inv_beta'):
        print(ds, [(n.split('\n')[0], f'{v[0]:.2e}', f'{v[1]:.2e}', 'published', f'{p[0]:.2e}', f'{p[1]:.2e}') for n, v, p in steps(ds)])
        r = rows[ds]
        print({k: f'{r[k]:.2e}' for k, _, _ in cats})


if __name__ == '__main__':
    main()
