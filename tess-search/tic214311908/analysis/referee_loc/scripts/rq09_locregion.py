"""Localisation regions from the combined real delta-chi2 map (baseline run, M-dwarf source colour) with thresholds
from M3V injections (quantiles of chi2(true M3V) - chi2(free)).  Output: out/locregion.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
tag = sys.argv[1] if len(sys.argv) > 1 else 'base'
Z = np.load(f'{ROOT}/out/inj_{tag}.npz')
GE, GN = Z['GE'], Z['GN']
m = Z['real_gsum'].reshape(len(GN), len(GE))
m = m - m.min()
thr_inj = np.percentile(Z['M3V__comb_true_minus_free'], [68.27, 95.45, 99.73])
EE, NN = np.meshgrid(GE, GN)
out = dict(tag=tag, min_pos=[float(EE.flat[np.argmin(m)]), float(NN.flat[np.argmin(m)])], M3V_value=float(m[np.argmin(np.abs(GN)), np.argmin(np.abs(GE))]))
g = gaia(); m3 = g[g.source_id == ID['M3V']].iloc[0]
# catalogued stars' 2025 offsets from the M3V (PM-corrected, epoch 2025.0) - map is co-moving with the M3V
dt = 9.0
rel = {k: ((g[g.source_id == v].iloc[0].dE - m3.dE) + (np.nan_to_num(g[g.source_id == v].iloc[0].pmra) - m3.pmra) * dt / 1000,
           (g[g.source_id == v].iloc[0].dN - m3.dN) + (np.nan_to_num(g[g.source_id == v].iloc[0].pmdec) - m3.pmdec) * dt / 1000)
       for k, v in ID.items() if k not in ('ref',)}
for lab, thr in zip(('68', '95', '99.7'), thr_inj):
    for kind, t in (('injection', thr), ('gaussian', {'68': 2.30, '95': 6.18, '99.7': 11.83}[lab])):
        inside = m <= t
        d = np.hypot(EE[inside], NN[inside])
        val = {k: float(m[np.argmin(np.abs(GN - p[1])), np.argmin(np.abs(GE - p[0]))]) for k, p in rel.items()}
        out[f'{lab}_{kind}'] = dict(threshold=float(t), area_arcsec2=float(inside.sum() * 0.25), max_dist_from_M3V=float(d.max()),
                                    E_range=[float(EE[inside].min()), float(EE[inside].max())], N_range=[float(NN[inside].min()), float(NN[inside].max())],
                                    stars_inside=[k for k, v in val.items() if v <= t])
        print(lab, kind, out[f'{lab}_{kind}'])
out['map_value_at_star_positions_Mdwarf_colour_frame'] = {k: float(m[np.argmin(np.abs(GN - p[1])), np.argmin(np.abs(GE - p[0]))]) for k, p in rel.items()}
out['star_offsets_from_M3V_2025'] = {k: [float(p[0]), float(p[1])] for k, p in rel.items()}
print(out['map_value_at_star_positions_Mdwarf_colour_frame'])
jdump(out, f'{ROOT}/out/locregion_{tag}.json')
