"""f7: followup_2027.png - transit windows observable from the southern sites in the 2027 and 2028 seasons, from the
write-up's own per-site tables (../observability/obs_<site>.csv, adopted ephemeris, written by scripts/w02_predictions.py):
window = T14 1.70 h + 2 x 3 sigma + 2 x 1.5 h baseline at airmass < 2; 'full' = Sun < -18 deg throughout,
'twilight' = Sun < -12 deg throughout; Moon illumination at mid-transit. Network counts from scripts/w02_checks.json.
Chile row = CTIO/Cerro Pachon, Las Campanas and La Silla, which have identical window lists (checked below)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

style()


def main():
    O = f'{W}/observability'
    files = {'Dome_C': 'Dome_C', 'Siding_Spring': 'Siding_Spring', 'SAAO': 'SAAO', 'Paranal': 'Paranal', 'CTIO': 'CTIO_SOAR_GeminiS',
             'Las_Campanas': 'LCO', 'La_Silla': 'La_Silla'}
    parts = []
    for site, f in files.items():
        x = pd.read_csv(f'{O}/obs_{f}.csv'); x['site'] = site; parts.append(x)
    d = pd.concat(parts, ignore_index=True)
    d['date'] = pd.to_datetime(d.utc_mid)
    d['full18'] = d.full_window.astype(bool); d['full12'] = d.full_window_sun12.astype(bool); d['frac18'] = d.frac_window_usable
    chile = ['CTIO', 'Las_Campanas', 'La_Silla']
    ref = d[d.site == 'CTIO'].set_index('epoch')
    for s_ in chile[1:]:
        o = d[d.site == s_].set_index('epoch')
        ix = ref.index.intersection(o.index)
        assert (o.loc[ix].full18 == ref.loc[ix].full18).all() and (o.loc[ix].full12 == ref.loc[ix].full12).all(), s_
        assert set(o[o.full18].index) == set(ref[ref.full18].index), s_
    rows = [('Dome_C', 'Dome C\n(ASTEP+)'), ('Siding_Spring', 'Siding Spring\n(MuSCAT4 on FTS, LCOGT)'), ('SAAO', 'SAAO Sutherland\n(LCOGT)'),
            ('Paranal', 'Paranal\n(SPECULOOS-S, VLT)'), ('CTIO', 'Chile: CTIO, Cerro Pachon,\nLas Campanas, La Silla')]
    chk = json.load(open(f'{W}/scripts/w02_checks.json'))
    seasons = [('2027-03-10', '2027-09-30'), ('2028-03-10', '2028-09-30')]
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.6), sharey=True, gridspec_kw=dict(wspace=0.05))
    cmap = plt.get_cmap('cividis')
    hi = {(126, 'Siding_Spring'), (127, 'CTIO'), (137, 'CTIO')}
    for ax, (a, b) in zip(axs, seasons):
        a, b = pd.Timestamp(a), pd.Timestamp(b)
        for i, (site, lab) in enumerate(rows):
            s = d[(d.site == site) & (d.date >= a) & (d.date <= b)]
            part = s[(~s.full12) & (s.frac18 >= 0.5)]
            ax.plot(part.date, [i] * len(part), '|', color=MUTED, ms=7, mew=1.0, zorder=2)
            f12 = s[s.full12 & ~s.full18]
            f18 = s[s.full18]
            ax.scatter(f12.date, [i] * len(f12), marker='^', s=58, c=f12.moon_illum, cmap=cmap, vmin=0, vmax=1, edgecolors=INK, linewidths=0.7, zorder=3)
            sc = ax.scatter(f18.date, [i] * len(f18), marker='o', s=70, c=f18.moon_illum, cmap=cmap, vmin=0, vmax=1, edgecolors=INK, linewidths=0.7, zorder=4)
            for _, r in f18.iterrows():
                if (int(r.epoch), site) in hi:
                    hhmm = (r.date + pd.Timedelta(seconds=30)).strftime('%H:%M')      # nearest minute, as in the README table
                    ax.annotate(f"E{int(r.epoch)}, {r.date.day} {r.date.strftime('%b')}\n{hhmm} UTC, Moon {r.moon_illum:.0%}", xy=(r.date, i),
                                xytext=(r.date, i + 0.36), fontsize=7.0, ha='center', va='bottom',
                                color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6), zorder=6,
                                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.9))
        ax.set_xlim(a, b)
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.set_title(f'{a.year} season', loc='left', fontsize=10)
        ax.grid(axis='y', visible=False)
        ax.set_ylim(-0.6, len(rows) + 0.05)
    axs[0].set_yticks(range(len(rows)))
    axs[0].set_yticklabels([r[1] for r in rows], fontsize=8)
    cb = fig.colorbar(sc, ax=axs, fraction=0.025, pad=0.015)
    cb.set_label('Moon illumination at mid-transit', fontsize=8.5)
    leg = [Line2D([], [], marker='o', ls='none', mfc='#c8c8c8', mec=INK, ms=8, label='whole window observable, Sun < $-$18$\\degree$'),
           Line2D([], [], marker='^', ls='none', mfc='#c8c8c8', mec=INK, ms=8, label='whole window only with Sun < $-$12$\\degree$'),
           Line2D([], [], marker='|', ls='none', color=MUTED, ms=8, mew=1.0, label='partial window ($\\geq$ 50% dark)')]
    fig.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.46, -0.075), ncol=3, fontsize=8)
    n27 = len(chk['full_epochs_2027_excl_DomeC']); n27d = len(chk['full_epochs_2027_incl_DomeC'])
    n28 = len(chk['full_epochs_2028_excl_DomeC']); n28d = len(chk['full_epochs_2028_incl_DomeC'])
    wmin, wmax = d.window_h.min(), d.window_h.max()
    fig.suptitle('TIC 214311908: transit windows from southern sites, 2027 and 2028 seasons', x=0.125, y=1.04, ha='left', fontsize=10.5)
    fig.text(0.125, 0.965, f'Window {wmin:.1f}-{wmax:.1f} h = T14 + $\\pm$3$\\sigma$ timing + 1.5-h baselines, all at altitude > 30$\\degree$ (adopted ephemeris).  '
             f'Distinct epochs with a full dark window: {n27} in 2027 ({n27d} including Dome C), {n28} in 2028 ({n28d}).',
             ha='left', fontsize=8.4, color=INK2)
    save(fig, 'followup_2027.png')
    for site, _ in rows:
        s = d[(d.site == site) & d.full18]
        print(site, list(zip(s.epoch, s.utc_mid, s.moon_illum.round(2))))


if __name__ == '__main__':
    main()
