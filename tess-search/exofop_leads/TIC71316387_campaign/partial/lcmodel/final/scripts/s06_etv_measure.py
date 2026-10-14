"""Per-event eclipse times (primary and secondary, every sector) with injection-calibrated empirical errors.
Each event: template (best-fit physical model, cadence-integrated) shifted in time, free depth scale and quadratic baseline
(linear least squares at each trial shift on a 0.1-min grid over +-12 min; parabolic refinement of chi2; formal error from
Delta chi2 = 1 with the sector's OOE scatter).
Calibration: the same template is injected into the real light curve at OOE phases 0.20/0.25/0.30/0.70/0.75/0.80 (both
primary and secondary shapes, random sub-cadence offsets) and re-measured; per (sector, eclipse type) we record the bias and
the robust scatter of (measured - injected)/formal  -> error scale factor (floored at 1).
Usage: s06_etv_measure.py <params.json> <tag> [default|tc|qlp]"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '2')
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import lclib as L, etvlib as E
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final'
pbest = json.load(open(sys.argv[1])); tag = sys.argv[2]; prodset = sys.argv[3] if len(sys.argv) > 3 else 'default'
P, T0 = pbest['P'], pbest['T0']
def lcs_for(ps):
    if ps == 'default': return L.default_set()
    if ps == 'tc': return [L.load_tc(s) for s in (1, 28, 68, 95, 105, 106)]
    if ps == 'qlp': return [L.load_qlp(s) for s in (1, 28, 68, 95, 105)]
rng = np.random.default_rng(7)
rows = []; inj = []
for lc in lcs_for(prodset):
    sec = lc['sector']; lc = L.clip_outliers(lc, P=P, T0=T0)
    scale = pbest.get(f's{sec}', 1.0) if prodset == 'default' else 1.0
    ws = L.windows(lc, P=P, T0=T0, hw=0.20)
    res = []
    for w in ws:
        i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
        c = np.polyfit(dt[oot], y[oot], 2); res.append(y[oot] / np.polyval(c, dt[oot]) - 1)
    r = np.concatenate(res); sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    for w in ws:
        i = w['idx']; tcm = T0 + (w['n'] + w['ph']) * P
        o = E.time_event(lc['t'][i], lc['f'][i], np.full(i.size, sig), tcm, pbest, lc['texp'], lc['nsub'], scale)
        o.update(sector=sec, prod=lc['prod'], n=int(w['n']), ph=float(w['ph']), tpred=float(tcm), sig=float(sig))
        rows.append(o)
    ninj = 0
    for ph0 in (0.20, 0.25, 0.30, 0.70, 0.75, 0.80):
        for w in L.windows(lc, P=P, T0=T0, hw=0.20, phases=(ph0,)):
            for kind in (0.0, 0.5):
                i = w['idx']; tinj = w['tc'] + rng.uniform(-0.5, 0.5) * lc['texp']
                tshift = (T0 + (w['n'] + kind) * P) - tinj
                m = E.template(pbest, lc['t'][i] + tshift, lc['texp'], lc['nsub'], scale)
                y = lc['f'][i] * m
                o = E.time_event(lc['t'][i] + tshift, y, np.full(i.size, sig), T0 + (w['n'] + kind) * P, pbest, lc['texp'], lc['nsub'], scale)
                o.update(sector=sec, kind=kind, n=int(w['n']), ph_inj=ph0); inj.append(o); ninj += 1
    print(f"S{sec} {lc['prod']}: sigma {sig*1e6:.0f} ppm; events {len(ws)}; injections {ninj}", flush=True)
json.dump(dict(rows=rows, inj=inj, P=P, T0=T0, params=pbest), open(f'{FD}/results/etv_{tag}.json', 'w'), indent=1)
for sec in sorted(set(r['sector'] for r in rows)):
    for ph in (0.0, 0.5):
        rr = [r for r in rows if r['sector'] == sec and r['ph'] == ph and not r['edge']]
        ii = [r for r in inj if r['sector'] == sec and r['kind'] == ph and not r['edge']]
        if not rr: continue
        x = np.array([r['dt_min'] for r in rr]); e = np.array([r['err_min'] for r in rr])
        xi = np.array([r['dt_min'] for r in ii]); ei = np.array([r['err_min'] for r in ii])
        z = xi / ei; zs = 1.4826 * np.median(np.abs(z - np.median(z)))
        print(f"S{sec} ph{ph}: N={len(rr)} median formal {np.median(e):.2f} min, real robust rms {1.4826*np.median(np.abs(x-np.median(x))):.2f} | "
              f"inj N={len(xi)} bias {np.median(xi):+.2f} robust sd {1.4826*np.median(np.abs(xi-np.median(xi))):.2f} min, z-scale {zs:.2f}")
