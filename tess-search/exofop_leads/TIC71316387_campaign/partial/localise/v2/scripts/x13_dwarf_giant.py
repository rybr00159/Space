"""Split the red end of the lateral-colour curve into dwarfs and giants (Gaia parallax, M_G)."""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from astropy.table import Table
from loclib import R, to_json
import x12_red_end as X   # recomputes frames (prints)
g = Table.read(f'{R}/gaia/cone3.5deg_G15.fits').to_pandas().set_index('source_id')
rows = []
for fr, f in zip(X.frames, sorted(glob.glob(f'{R}/results/x10w_*.npz'))):
    a = np.load(f)['stars']
    a = a[(a['chi2r'] < 20) & (a['ruwe'] < 1.4) & (a['peak'] < 1.5e5) & np.isfinite(a['C'])]
    assert len(a) == len(fr['C'])
    plx = g.reindex(a['sid']).parallax.values
    MG = a['G'] + 5 * np.log10(np.clip(plx, 1e-3, None) / 100.)
    rows.append(pd.DataFrame(dict(C=a['C'], rad=fr['rad'] * X.kref / fr['k_rad'], keep=fr['keep'], MG=MG, plx=plx, sec=fr['sec'], sid=a['sid'])))
d = pd.concat(rows)
d['cls'] = np.where(d.MG > 4.0 + 1.5 * (d.C - 0.8).clip(lower=0), 'dwarf', 'giant')
out = []
for lo, hi in [(1.25, 1.5), (1.5, 1.75), (1.75, 2.0), (2.0, 2.3), (2.3, 2.6), (2.6, 3.0), (3.0, 3.6), (3.6, 4.5)]:
    for cls in ['dwarf', 'giant']:
        s = d.keep & (d.C >= lo) & (d.C < hi) & (d.cls == cls)
        if s.sum() < 2:
            continue
        med = np.median(d.rad[s]); mad = 1.4826 * np.median(np.abs(d.rad[s] - med))
        n_u = d.sid[s].nunique()
        out.append(dict(lo=lo, hi=hi, cls=cls, n=int(s.sum()), n_unique=int(n_u), C_med=float(np.median(d.C[s])), rad_mpx=float(med * 1e3),
                        err_mpx=float(1.2533 * mad / np.sqrt(max(n_u, 1)) * 1e3)))
        o = out[-1]
        print(f"  C {lo:.2f}-{hi:.2f} {cls:5s}: n={o['n']:4d} (unique {o['n_unique']:3d}) med C {o['C_med']:.2f} rad {o['rad_mpx']:+7.1f} +- {o['err_mpx']:5.1f} mpx")
to_json(out, f'{R}/results/x13_dwarf_giant.json')
