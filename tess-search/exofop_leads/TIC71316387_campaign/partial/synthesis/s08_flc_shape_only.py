# Synthesis s08: sensitivity of the light-fraction factor F_LC to dropping the main-sequence
# consistency term (shape-only L(f), referee r09 'dchi2_lc'), e.g. if the M dwarfs were young/inflated.
import numpy as np, json
from astropy.table import Table
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
rng = np.random.default_rng(3)
t = Table.read(B + 'synthesis/field_Mdwarfs_1p5deg.ecsv')
ra0, de0 = 330.887948, -42.835511
ra = np.radians(np.array(t['ra'])); de = np.radians(np.array(t['dec']))
sep = np.degrees(np.arccos(np.clip(np.sin(np.radians(de0)) * np.sin(de) + np.cos(np.radians(de0)) * np.cos(de) * np.cos(ra - np.radians(ra0)), -1, 1)))
br = np.array(t['bp_rp'], float); G = np.array(t['phot_g_mean_mag'], float)
Tm = G - 0.00522555 * br**3 + 0.0891337 * br**2 - 0.633923 * br + 0.0324473
f = 10**(-0.4 * (Tm - 11.358))
r9 = json.load(open(B + 'referee_all/r09_msprior.json'))
out = {}
for col in ('dchi2_lc', 'total'):
    for s in ('A_Mann', 'B_PM_table'):
        tab = r9[s]; ff = np.array([x['f'] for x in tab]); c = np.array([x[col] for x in tab]); dist = np.array([x['dist_pc'] for x in tab])
        Lf = lambda x: np.exp(-0.5 * (np.interp(np.log(np.clip(x, 1e-4, 1)), np.log(ff), c, left=1e3, right=1e3) - c.min()))
        lnf297 = np.interp(np.log(297.4), np.log(dist[::-1]), np.log(ff[::-1]))
        Lb = float(Lf(np.exp(rng.normal(lnf297, 0.15 / 1.0857, 100000))).mean())
        res = {}
        for lab, (c1, c2) in {'2.1-2.5': (2.1, 2.5), '1.8-2.7': (1.8, 2.7)}.items():
            m = (sep > 1 / 60) & (br >= c1) & (br <= c2) & (f >= 0.022) & (f <= 0.065)
            res[lab] = dict(n=int(m.sum()), mean_L_hosts=round(float(Lf(f[m]).mean()), 3), F_LC=round(Lb / float(Lf(f[m]).mean()), 2))
        out[f'{col}_{s}'] = dict(f_at_297=round(float(np.exp(lnf297)), 4), L_bound=round(Lb, 3), **res)
json.dump(out, open(B + 'synthesis/s08_flc_shape_only.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
