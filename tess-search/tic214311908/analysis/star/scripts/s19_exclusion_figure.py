"""Combined exclusion map for companions/blends of TIC 214311908 in (projected separation or semi-major axis [au],
Delta T [TESS mag]), from: Gaia astrometry (s12c; masses converted to Delta T with the Mamajek MS table and the
Stassun+19 G->T relation), Gaia IPD diagnostics (s12b, red-calibrator variant; contrast in I ~ T), DECaPS2 residual
contrast (s13c; mean of i and z, >= 1"), and the DSS look-back test (s13b; unbound background stars only, within 3"
of the TESS-epoch positions). Horizontal lines: the faintest source that can supply the 15.4 e-/s flux loss with a
100% / 50% eclipse (Delta T = 4.46 / 3.71). Output: figs/companion_exclusion.png, out/exclusion_summary.json"""
import json, numpy as np
B = 'tic214311908/analysis/star'
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
dist = 64.31
ast = json.load(open(f'{B}/out/astrometric_limits.json')); ipd = json.load(open(f'{B}/out/ipd_limits.json')); dc = json.load(open(f'{B}/out/decaps_contrast.json'))
dss = json.load(open(f'{B}/out/dss_limits.json'))
tab = []
for line in open(f'{B}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] in 'KM' and p[0].endswith('V'):
        try: tab.append((float(p[30]), float(p[13]), float(p[11])))
        except Exception: pass
tab = np.array(sorted(tab))
def MT(M):
    mg = np.interp(M, tab[:, 0], tab[:, 1]); c = np.interp(M, tab[:, 0], tab[:, 2])
    return mg - (0.00522555 * c ** 3 - 0.0891337 * c ** 2 + 0.633923 * c - 0.0324473)
MT1 = 10.233 - (0.00522555 * 2.756 ** 3 - 0.0891337 * 2.756 ** 2 + 0.633923 * 2.756 - 0.0324473)
fig, ax = plt.subplots(figsize=(8, 5))
M2 = np.array(ast['M2']); a = np.array(ast['a_au']); L = np.array(ast['L_combined'])
dT = np.array([MT(m) - MT1 if m >= 0.078 else 12.0 for m in M2])
summ = dict(astrometric=[])
for i in range(len(M2)):
    ex = a[L[i] < 0.05]
    if len(ex):
        ax.plot([ex.min(), ex.max()], [dT[i], dT[i]], color='C0', lw=6, alpha=0.6, solid_capstyle='butt')
        summ['astrometric'].append(dict(M2=float(M2[i]), dT=float(dT[i]), a_excluded_au=[float(ex.min()), float(ex.max())]))
ax.plot([], [], color='C0', lw=6, alpha=0.6, label='Gaia RUWE + DR2-DR3 PM (bound; L<0.05)')
col = {'excluded': 'C3', 'disfavoured': 'C1'}
for c in ipd['red_variant']['cells']:
    if c['status'] in col:
        ax.add_patch(plt.Rectangle((c['sep'][0] * dist, c['dm'][0]), (c['sep'][1] - c['sep'][0]) * dist, c['dm'][1] - c['dm'][0], color=col[c['status']], alpha=0.45, lw=0))
ax.add_patch(plt.Rectangle((0, 0), 0, 0, color='C3', alpha=0.45, label='Gaia IPD/RUWE image diagnostics: excluded'))
ax.add_patch(plt.Rectangle((0, 0), 0, 0, color='C1', alpha=0.45, label='... disfavoured (LR 0.1-0.33)'))
ci = {r['sep']: r['delta_mag_5sig'] for r in dc['i']['curve'] if 'delta_mag_5sig' in r}; cz = {r['sep']: r['delta_mag_5sig'] for r in dc['z']['curve'] if 'delta_mag_5sig' in r}
ss = sorted(s for s in ci if s >= 1.0 and s in cz)
ax.fill_between(np.array(ss) * dist, 0, [0.5 * (ci[s] + cz[s]) for s in ss], color='C2', alpha=0.3, label='DECaPS2 2018 residual image (5 sigma)')
summ['decaps_dT_5sig'] = {f'{s}"': 0.5 * (ci[s] + cz[s]) for s in ss}
for y, ls, lab in ((4.46, '--', '100% eclipse limit'), (3.71, ':', '50% eclipse limit')):
    ax.axhline(y, color='k', ls=ls, lw=1); ax.text(0.012, y - 0.1, lab, fontsize=7)
ax.set_xscale('log'); ax.set_xlim(0.01, 500); ax.set_ylim(8, 0)
ax.set_xlabel('semi-major axis (astrometry) / projected separation (imaging), au'); ax.set_ylabel(r'$\Delta T$ (companion - M3V)')
ax.set_title('TIC 214311908: companions excluded by archival data (bound and unbound)', fontsize=10)
s26 = {k: v['grid_T_lim_min'] for k, v in dss.items()}
ax.text(0.012, 7.6, 'DSS 1979/1992 look-back: no unbound star with T < 18.4-19.1\nwithin 3" of the 2019-2026 positions (Delta T > 5.3)', fontsize=7)
ax.legend(fontsize=7, loc='lower right')
fig.tight_layout(); fig.savefig(f'{B}/figs/companion_exclusion.png', dpi=130)
summ['dss_grid_T_lim_min'] = s26; summ['MT_primary'] = float(MT1)
json.dump(summ, open(f'{B}/out/exclusion_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=0)[:2500])
