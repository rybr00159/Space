"""Light-curve-level positive control for the ETV analysis: in every real eclipse window divide out the best-fit template
and multiply by the template shifted by dt_inj = c2 * tyr^2 (tyr = years from BTJD 2900), i.e. a quadratic O-C written
into the real photometry with its real noise; re-measure every eclipse time with the same pipeline (s06 settings) and fit
linear+quadratic with the per-event errors scaled by the injection-calibrated z-scales of the real run.
Usage: s07b_etv_lcinj.py <params.json> <etv_real.json>"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import lclib as L, etvlib as E
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final'
pbest = json.load(open(sys.argv[1])); P, T0 = pbest['P'], pbest['T0']
real = json.load(open(sys.argv[2]))
zsc = {}
for sec in set(r['sector'] for r in real['inj']):
    for ph in (0.0, 0.5):
        ii = [r for r in real['inj'] if r['sector'] == sec and r['kind'] == ph and not r['edge']]
        z = np.array([r['dt_min'] / r['err_min'] for r in ii]); zsc[(sec, ph)] = max(1.0, 1.4826 * np.median(np.abs(z - np.median(z))))
lcs = [L.clip_outliers(lc, P=P, T0=T0) for lc in L.default_set()]
out = []
for c2 in [0.0, 0.1, 0.2, 0.4]:
    E_, oc, es, ph_ = [], [], [], []
    for lc in lcs:
        sec = lc['sector']; scale = pbest.get(f's{sec}', 1.0)
        ws = L.windows(lc, P=P, T0=T0, hw=0.20)
        res = []
        for w in ws:
            i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
            c = np.polyfit(dt[oot], y[oot], 2); res.append(y[oot] / np.polyval(c, dt[oot]) - 1)
        r = np.concatenate(res); sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        for w in ws:
            i = w['idx']; t = lc['t'][i]; tcm = T0 + (w['n'] + w['ph']) * P
            sh = c2 * ((tcm - 2900.0) / 365.25) ** 2 / 1440.0
            m0 = E.template(pbest, t, lc['texp'], lc['nsub'], scale); m1 = E.template(pbest, t - sh, lc['texp'], lc['nsub'], scale)
            o = E.time_event(t, lc['f'][i] / m0 * m1, np.full(i.size, sig), tcm, pbest, lc['texp'], lc['nsub'], scale)
            if o['edge']: continue
            E_.append(w['n'] + w['ph']); oc.append(o['dt_min']); es.append(o['err_min'] * zsc[(sec, w['ph'])]); ph_.append(w['ph'])
    E_ = np.array(E_); oc = np.array(oc); es = np.array(es); ph_ = np.array(ph_)
    tyr = (T0 + E_ * P - 2900) / 365.25
    A = np.array([np.ones_like(E_), E_ / 1000, (ph_ == 0.5).astype(float), tyr ** 2]).T
    Aw = A / es[:, None]; c, *_ = np.linalg.lstsq(Aw, oc / es, rcond=None); cov = np.linalg.inv(Aw.T @ Aw)
    r = (oc - A @ c) / es; s = max(1, np.sqrt(r @ r / (len(oc) - 4)))
    out.append(dict(c2_inj=c2, c2_rec=float(c[3]), c2_err=float(np.sqrt(cov[3, 3]) * s), n=len(oc), chi2=float(r @ r)))
    print(f'c2 injected {c2:.2f} min/yr^2 -> recovered {c[3]:+.3f} +- {np.sqrt(cov[3,3])*s:.3f}  (N={len(oc)}, chi2 {r@r:.1f})', flush=True)
json.dump(out, open(f'{FD}/results/etv_lcinj.json', 'w'), indent=1)
