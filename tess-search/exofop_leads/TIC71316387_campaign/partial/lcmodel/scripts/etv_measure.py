"""Measure per-event eclipse times for primaries and secondaries in all sectors, plus injection-calibrated
empirical timing errors (template injected at OOE phases 0.25/0.75 in the real data)."""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np, sys, json, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import lclib as L, fit2, etvlib as E
pbest = json.load(open(sys.argv[1]))
tag = sys.argv[2] if len(sys.argv) > 2 else 'default'
prodset = sys.argv[3] if len(sys.argv) > 3 else 'default'
P, T0 = pbest['P'], pbest['T0']
def lcs_for(prodset):
    if prodset == 'default': return L.default_set()
    if prodset == 'tc': return [L.load_tc(s) for s in (1, 28, 68, 95, 105, 106)]
    if prodset == 'qlp': return [L.load_qlp(s) for s in (1, 28, 68, 95, 105)]
rng = np.random.default_rng(7)
rows = []; inj = []
for lc in lcs_for(prodset):
    sec = lc['sector']; lc = L.clip_outliers(lc, P=P, T0=T0)
    scale = pbest.get(f's{sec}', 1.0)
    # sector noise: OOE scatter from windows
    ws = L.windows(lc, P=P, T0=T0, hw=0.20)
    res = []
    for w in ws:
        i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
        c = np.polyfit(dt[oot], y[oot], 2); res.append(y[oot] / np.polyval(c, dt[oot]) - 1)
    r = np.concatenate(res); sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    for w in ws:
        i = w['idx']
        tcm = T0 + (w['n'] + w['ph']) * P   # template conjunction time for this event
        # template's own eclipse centre for secondary includes ecosw offset -> measure relative to template
        o = E.time_event(lc['t'][i], lc['f'][i], np.full(i.size, sig), tcm, pbest, lc['texp'], lc['nsub'], scale)
        o.update(sector=sec, prod=lc['prod'], n=int(w['n']), ph=float(w['ph']), tpred=float(tcm), sig=float(sig))
        rows.append(o)
    # injections at OOE phases (0.25, 0.75) with a random sub-cadence jitter of the injected time
    wsi = L.windows(lc, P=P, T0=T0, hw=0.20, phases=(0.25, 0.75))
    for w in wsi:
        for kind in (0.0, 0.5):
            i = w['idx']; tinj = w['tc'] + rng.uniform(-0.5, 0.5) * lc['texp']
            # model eclipse of type 'kind' placed at tinj: evaluate template at conjunction (kind) shifted
            tshift = (T0 + (w['n'] + kind) * P) - tinj
            m = E.template(pbest, lc['t'][i] + tshift, lc['texp'], lc['nsub'], scale)
            y = lc['f'][i] * m
            o = E.time_event(lc['t'][i] + tshift, y, np.full(i.size, sig), T0 + (w['n'] + kind) * P, pbest, lc['texp'], lc['nsub'], scale)
            o.update(sector=sec, kind=kind, n=int(w['n']))
            inj.append(o)
    print(f"S{sec} {lc['prod']}: sigma {sig*1e6:.0f} ppm; events {len(ws)}; injections {len(wsi)*2}", flush=True)
json.dump(dict(rows=rows, inj=inj, P=P, T0=T0), open(L.HERE + f'/results/etv_{tag}.json', 'w'), indent=1)
# summary
import collections
for sec in sorted(set(r['sector'] for r in rows)):
    for ph in (0.0, 0.5):
        rr = [r for r in rows if r['sector'] == sec and r['ph'] == ph and not r['edge']]
        ii = [r for r in inj if r['sector'] == sec and r['kind'] == ph and not r['edge']]
        if not rr: continue
        x = np.array([r['dt_min'] for r in rr]); e = np.array([r['err_min'] for r in rr])
        xi = np.array([r['dt_min'] for r in ii]); ei = np.array([r['err_min'] for r in ii])
        emp = 1.4826 * np.median(np.abs(xi - np.median(xi))) if len(xi) > 3 else np.nan
        print(f"S{sec} ph{ph}: N={len(rr)} median formal {np.median(e):.2f} min, real O-C rms(robust) {1.4826*np.median(np.abs(x-np.median(x))):.2f} min, wmean {np.average(x,weights=1/e**2):+.2f}; "
              f"inj N={len(xi)} bias {np.median(xi):+.2f} robust-sd {emp:.2f} sd {np.std(xi):.2f} (median formal {np.median(ei):.2f}) -> ratio {emp/np.median(ei):.2f}")
