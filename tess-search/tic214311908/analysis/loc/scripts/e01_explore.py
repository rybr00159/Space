"""Exploration: geometry per sector (camera, field angle, boresight PA, Jacobian), pixel positions of key stars,
median image, and per-pixel noise in 200-s data."""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
g = stars()
key = g[g.source_id.isin(IDS.values())]
for sec in SECTORS:
    S = Sector(sec)
    xy = S.star_xy(key, colour=False)
    print(f'S{sec}: cam{S.cam} ccd{S.ccd} theta={S.theta:.2f} PA_bs={S.pa_bs:.1f} colscale={S.colscale:.3f} exp={S.exptime*1440:.1f}min '
          f'ngood={S.good.sum()} J=[{S.J[0,0]:.2f} {S.J[0,1]:.2f}; {S.J[1,0]:.2f} {S.J[1,1]:.2f}] epoch={S.epoch:.2f}')
    for (i, r), p in zip(key.iterrows(), xy):
        nm = [k for k, v in IDS.items() if v == r.source_id][0]
        print(f'    {nm:6s} x={p[0]:.2f} y={p[1]:.2f}  T={r["T"]:.2f} C={r.C:.2f} F={r.F:.0f}')
    med = np.median(S.flux[S.good], 0)
    iy, ix = np.unravel_index(np.argmax(med), med.shape)
    print(f'    median image: max {med.max():.0f} at (x={ix}, y={iy}); 5x5 around Mira sum {med[int(round(S.ym0))-2:int(round(S.ym0))+3, int(round(S.xm0))-2:int(round(S.xm0))+3].sum():.0f}; min {med.min():.0f}')
    # flux of the central 3x3 over the sector: min/max (Mira variability)
    yc, xc = int(round(S.ym0)), int(round(S.xm0))
    lc = S.flux[S.good][:, yc - 1:yc + 2, xc - 1:xc + 2].sum((1, 2))
    print(f'    central 3x3 flux: 5th pct {np.percentile(lc, 5):.0f}, 95th {np.percentile(lc, 95):.0f}')
