"""Localisation map of the real signal assuming a normal-colour source (BP-RP 1.3, like the G = 17.5 star 3.2" from
the M3V) instead of the M-dwarf colour 2.76: where would a blue(r) uncatalogued source have to be?  Real data only;
thresholds taken from the M3V injections (out/post_*_prf_all.json).
Output: out/colourmap.json"""
import sys, json
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
out = {}
g = stars(); r3 = g[g.source_id == IDS['M3V']].iloc[0]
for meth in ('DIFF', 'AMP'):
    T = json.load(open(f'{BASE}/out/post_{meth.lower()}_prf_all.json'))['loc_thresholds_used']
    for col in (2.76, 1.3):
        gs = None
        for sec in SECTORS:
            M = SectorModel(sec, meth, grid_colour=col, grid_half=30.0, grid_step=0.5)
            cg, *_ = M.fit(M.d_real, 'G'); gs = cg if gs is None else gs + cg
            gr = M.grid
        d = gs - gs.min(); i0 = int(np.argmin(gs))
        res = dict(min_offset_from_M3V=gr[i0].tolist())
        for lab, tv in zip(('68', '95', '99.7'), T):
            m = d <= tv
            res[lab] = dict(max_dist=float(np.hypot(*gr[m].T).max()), E=[float(gr[m, 0].min()), float(gr[m, 0].max())], N=[float(gr[m, 1].min()), float(gr[m, 1].max())])
        # delta-chi2 at the faint star's sky position (mean epoch) in this colour frame
        s = g[g.source_id == 5925373272827260544].iloc[0]; ep = 2023.6 - 2016
        pE = (s.dE - r3.dE) + (s.pmra - r3.pmra) * ep / 1000; pN = (s.dN - r3.dN) + (s.pmdec - r3.pmdec) * ep / 1000
        j = int(np.argmin(np.hypot(gr[:, 0] - pE, gr[:, 1] - pN)))
        res['dchi2_at_260544_skypos'] = float(d[j]); res['dchi2_at_M3V'] = float(d[int(np.argmin(np.hypot(*gr.T)))])
        out[f'{meth}_C{col}'] = res
        print(f"{meth} colour {col}: minimum at {np.round(gr[i0], 1).tolist()}\" from M3V; 68/95/99.7% max distance "
              f"{res['68']['max_dist']:.1f}/{res['95']['max_dist']:.1f}/{res['99.7']['max_dist']:.1f}\"; delta-chi2 at M3V {res['dchi2_at_M3V']:.1f}, "
              f"at the sky position of the G17.5 star (3.2\") {res['dchi2_at_260544_skypos']:.1f}", flush=True)
to_json(out, f'{BASE}/out/colourmap.json')
