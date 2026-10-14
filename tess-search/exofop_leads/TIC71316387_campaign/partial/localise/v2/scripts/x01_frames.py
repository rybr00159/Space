import json, numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/scripts')
from loclib import *
from astropy.io import fits
import tess_stars2px as t2p
for sec in SECTORS:
    d = json.load(open(f'{R}/results/real_s{sec}_prf.json'))
    J = np.array(d['jac']); Ji = np.linalg.inv(J)
    with fits.open(cutout_file(sec)) as h:
        col0 = int(h[1].header['1CRV4P']); row0 = int(h[1].header['2CRV4P'])
    # roll: PA of +x and +y pixel axes on sky (deg E of N)
    pax = np.degrees(np.arctan2(J[0, 0], J[1, 0])); pay = np.degrees(np.arctan2(J[0, 1], J[1, 1]))
    # focal-plane coords: use tess-point to get camera/ccd/col/row; then position relative to camera boresight
    res = t2p.tess_stars2px_function_entry(1, 330.887948, -42.835511, trySector=sec)
    outID, outEclLon, outEclLat, outSec, outCam, outCcd, outColPix, outRowPix, scinfo = res
    # boresight direction: compute sky PA from star towards camera centre using tess-point's camera pointing
    camra = scinfo.ra_ab if hasattr(scinfo,'ra_ab') else None
    line = f"S{sec:3d} cam{outCam[0]} ccd{outCcd[0]} col {outColPix[0]:.0f} row {outRowPix[0]:.0f}  PA(+x)={pax:6.1f} PA(+y)={pay:6.1f} scale={np.sqrt(abs(np.linalg.det(J))):.2f}\"/px"
    print(line)
    for k in ['deep', 'shallow', 'all']:
        r = d[k]
        dx, dy = Ji @ np.array([r['dE'], r['dN']])
        print(f"   {k:7s} sky ({r['dE']:+.2f},{r['dN']:+.2f})\"  pix ({dx:+.3f},{dy:+.3f}) px   [{dx*21:+.2f},{dy*21:+.2f}]\" in pixel axes")
