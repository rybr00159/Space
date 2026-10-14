import json, numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/scripts')
from loclib import *
import tess_stars2px as t2p
from astropy.coordinates import SkyCoord
import astropy.units as u
out = {}
ra0, de0 = 330.887948, -42.835511
for sec in SECTORS:
    d = json.load(open(f'{R}/results/real_s{sec}_prf.json'))
    J = np.array(d['jac']); Ji = np.linalg.inv(J)
    res = t2p.tess_stars2px_function_entry(1, ra0, de0, trySector=sec)
    sc = res[-1]; cam = res[4][0]
    i = list(sc.sectors).index(sec)
    cra, cde = sc.camRa[cam - 1][i], sc.camDec[cam - 1][i]
    s0 = SkyCoord(ra0 * u.deg, de0 * u.deg); c = SkyCoord(cra * u.deg, cde * u.deg)
    pa = s0.position_angle(c).deg; sep = s0.separation(c).deg
    # unit vector toward boresight in sky (E,N) and pixel frames
    uE, uN = np.sin(np.radians(pa)), np.cos(np.radians(pa))
    up = Ji @ np.array([uE, uN]); up /= np.hypot(*up)
    rows = []
    for k in ['deep', 'shallow', 'all']:
        r = d[k]; v = np.array([r['dE'], r['dN']])
        rad = v @ np.array([uE, uN]); tan = v @ np.array([uN, -uE])
        rows.append((k, rad, tan))
    out[sec] = dict(cam=int(cam), field_angle_deg=sep, pa_to_boresight=pa, dir_pix_to_boresight=up.tolist(),
                    comps={k: dict(radial_toward_bs=a, tangential=b) for k, a, b in rows})
    print(f"S{sec:3d} cam{cam} field angle {sep:5.2f} deg, PA->boresight {pa:6.1f}, pixel dir->bs ({up[0]:+.2f},{up[1]:+.2f}) | " +
          '  '.join(f"{k}: rad {a:+.2f} tan {b:+.2f}" for k, a, b in rows))
to_json(out, f'{R}/results/x02_boresight.json')
