"""Referee: recompute the DSS look-back T limits from the investigator's own measured quantities
(star/out/dss_limits.json: per-grid-point baseline S/N, noise sigma and offset, magnitude-amplitude calibration),
with a correct upper-limit rule.
Investigator rule (s13b): amplitude a source must ADD to reach detection = off + 3 sigma - base  -> T_lim.
  For grid points that already have a positive residual (base > off) this makes the limit fainter, and for
  base > off + 3 sigma it returns T_lim = 24 (no constraint shown, but the point is silently treated as 'very deep').
Referee rule: 3-sigma upper limit on the amplitude of a source at that position = max(base - off, 0) + 3 sigma.
Points with base > off + 3 sigma are listed as 'positive residual' (a possible source / PSF-wing residual).
Output: out/r03_dss_limits.json"""
import json, numpy as np
S = 'tic214311908/analysis/star/out/dss_limits.json'
O = 'tic214311908/analysis/referee_star/out/r03_dss_limits.json'
d = json.load(open(S))
def G_of_T(T, c): return T + 0.00522555 * c ** 3 - 0.0891337 * c ** 2 + 0.633923 * c - 0.0324473
out = {}
for lab, v in d.items():
    if not isinstance(v, dict) or 'grid_limits' not in v: continue
    cc = v['calib']['coef']; sig = v['noise_amp']; off = v['noise_offset']
    Ts = np.linspace(10, 24, 2801)
    def amp(T, c):
        G = G_of_T(T, c); return 10 ** (cc[0] + cc[1] * (G - 18) + cc[2] * (G - 18) ** 2 + cc[3] * (c - 1.2))
    res = {}
    for c in (1.0, 2.5):
        a = amp(Ts, c)
        # amplitude must be monotonic decreasing over Ts used
        ok = np.r_[True, np.diff(a) < 0]
        rows = [x for x in v['grid_limits'] if x['bp_rp'] == c]
        lims, pos = [], []
        for x in rows:
            bs = x['baseline_snr']
            up = (max(bs, 0) + 3.0) * sig
            Tl = float(np.interp(-np.log10(up), -np.log10(a[ok]), Ts[ok]))
            lims.append(dict(epoch=x['epoch'], r=x['r'], th=x['th'], baseline_snr=bs, T_lim_referee=Tl, T_lim_investigator=x['T_lim'], sep_from_star_then=x['sep_from_star_then']))
            if bs > 3: pos.append(lims[-1])
        L = np.array([l['T_lim_referee'] for l in lims]); R = np.array([l['r'] for l in lims])
        Li = np.array([l['T_lim_investigator'] for l in lims])
        res[f'bp_rp={c}'] = dict(min_within3_referee=float(L.min()), min_within2_referee=float(L[R <= 2].min()),
                                 min_within3_investigator=float(Li.min()), min_within2_investigator=float(Li[R <= 2].min()),
                                 median_referee=float(np.median(L)), n_points=len(L), positive_residual_points=pos)
    out[lab] = dict(calib_sd_dex=v['calib']['robust_sd_dex'], sigma=sig, offset=off, limits=res)
    for c in (1.0, 2.5):
        r_ = res[f'bp_rp={c}']
        print(f"{lab:16s} BP-RP {c}: within 3\": investigator {r_['min_within3_investigator']:.2f}  referee {r_['min_within3_referee']:.2f}; "
              f"within 2\": inv {r_['min_within2_investigator']:.2f} ref {r_['min_within2_referee']:.2f}; positive-residual points (>3 sigma): "
              + '; '.join(f"r={p['r']:.0f}\" ep {p['epoch']} S/N {p['baseline_snr']:.1f} -> T_lim {p['T_lim_referee']:.2f} (inv {p['T_lim_investigator']:.1f}), {p['sep_from_star_then']:.1f}\" from star then" for p in r_['positive_residual_points']))
json.dump(out, open(O, 'w'), indent=1)
