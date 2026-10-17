"""Time-system check: my target-direction BJD vs TESScut TIME, and vs the loc agent's correction."""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import l0lib
out = {}
for s in SECTORS:
    S = Sec(s, cube=False)
    L = l0lib.Sector(s, load_cube=False)
    ok = S.good & np.isfinite(L.t)
    d = (S.t[ok] - L.t[ok]) * 86400
    out[s] = dict(my_minus_tesscut_s=[float(np.nanmin(S.dcorr_s[ok])), float(np.nanmax(S.dcorr_s[ok]))],
                  my_minus_loc_s=[float(d.min()), float(d.max())], cadence_s=S.dt * 86400, epoch=S.epoch,
                  cam=S.cam, ccd=S.ccd, theta=S.theta, u=S.u.tolist())
    print(s, out[s])
jdump(out, f'{ROOT}/out/time_check.json')
