"""Model-light shape test: is the SECONDARY eclipse flat-bottomed (total occultation of the smaller star) or V-shaped (grazing)?
Each eclipse type gets its own symmetric trapezoid: depth, total duration T14, flat fraction f = T23/T14 in [0,1); the secondary
centre is free (offset from phase 0.5).  Cadence integration by supersampling, per-event quadratic baselines profiled (as in the
physical fit), per-sector depth scales s95/s105/s106 fixed at the physical best fit.  Statistic: Delta chi2 between the
best V-shaped secondary (f_sec = 0) and the best free f_sec; also the profile in f_sec.  The same statistic is computed on the
injection data sets of s04 (real noise, truths = total-occultation best fit / grazing L3=0 / partial L3=0.90), so the observed value
is calibrated empirically against both hypotheses.
Usage: s14_trapezoid.py [inject tag]"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, json, glob, pickle, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, lclib as L
from scipy.optimize import least_squares
from multiprocessing import get_context
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
best = json.load(open(f'{FD}/results/best_geo.json'))
SC = {s: best.get(f's{s}', 1.0) for s in fit3.SECS}

def trap(dt, depth, T14, f):
    """unit-normalised trapezoid profile (1 - depth*shape), dt from centre"""
    T23 = f * T14; a = np.abs(dt)
    ing = 0.5 * (T14 - T23)
    sh = np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max(ing, 1e-9)))
    return depth * sh

def model(d, th, P, T0):
    dp, T14p, fp, ds, T14s, fs, dts = th
    m = np.empty(d['t'].size)
    for sec, idx, texp, nsub in d['groups']:
        t = d['t'][idx]; acc = np.zeros(idx.size)
        for j in range(nsub):
            tt = t + texp * ((j + 0.5) / nsub - 0.5)
            ph = ((tt - T0) / P + 0.25) % 1 - 0.25          # -0.25..0.75
            dtp = ph * P; dtsec = (ph - 0.5) * P - dts / 1440
            acc += np.where(np.abs(ph) < 0.25, trap(dtp, dp, T14p, fp), trap(dtsec, ds, T14s, fs))
        m[idx] = 1 - SC[sec] * acc / nsub
    return m

def fit(d, P, T0, fs_fix=None, starts=((0.4, 0.4), (0.1, 0.1), (0.7, 0.6))):
    bestr = None
    for fp0, fs0 in starts:
        x0 = [0.0058, 0.064, fp0, 0.0027, 0.064, fs0 if fs_fix is None else fs_fix, 0.0]
        lb = [0.001, 0.03, 0.0, 0.0005, 0.03, 0.0, -20]; ub = [0.02, 0.10, 0.95, 0.01, 0.10, 0.95, 20]
        free = [i for i in range(7) if not (i == 5 and fs_fix is not None)]
        def rf(x):
            th = np.array(x0, float); th[free] = x
            return fit3.profile_resid(d['y'], model(d, th, P, T0), d['dt'], 1 / d['sig'], d['starts'], d['ends'])
        xs = np.array(x0)[free]; l = np.array(lb)[free]; u = np.array(ub)[free]
        xs = np.clip(xs, l + 1e-6, u - 1e-6)
        r = least_squares(rf, xs, bounds=(l, u), x_scale=np.array([5e-4, 2e-3, 0.05, 5e-4, 2e-3, 0.05, 0.5])[free], diff_step=1e-3, max_nfev=400)
        th = np.array(x0, float); th[free] = r.x
        c2 = float(2 * r.cost)
        if bestr is None or c2 < bestr[0]: bestr = (c2, th.tolist())
    return bestr

def analyse(d, P, T0):
    free = fit(d, P, T0)
    v0 = fit(d, P, T0, fs_fix=0.0)
    prof = {fs: fit(d, P, T0, fs_fix=fs, starts=((0.4, fs),))[0] - free[0] for fs in (0.1, 0.2, 0.3, 0.4, 0.5)}
    return dict(chi2_free=free[0], th_free=free[1], chi2_V=v0[0], th_V=v0[1], D_flat=v0[0] - free[0], prof=prof)

def job(a):
    path, P, T0, lab = a
    d = fit3.load_data(path)
    o = analyse(d, P, T0); o['label'] = lab
    return o

if __name__ == '__main__':
    t0 = time.time()
    tag = sys.argv[1] if len(sys.argv) > 1 else 'A'
    jobs = [(f'{FD}/results/data_default.pkl', best['P'], best['T0'], 'real')]
    for f in sorted(glob.glob(f'{FD}/results/inj_cache/{tag}_*.pkl')):
        nm = os.path.basename(f)[:-4]; phi = float(nm.split('_')[-1]); tn = nm[len(tag) + 1:-(len(nm.split('_')[-1]) + 1)]
        jobs.append((f, L.P0, L.T00 + phi * L.P0, f'{tn}|{phi:.2f}'))
    with get_context('spawn').Pool(int(os.environ.get('NPROC', '4'))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    json.dump(res, open(f'{FD}/results/trapezoid_{tag}.json', 'w'), indent=1)
    real = res[0]
    th = real['th_free']
    print(f"REAL: free fit: prim depth {th[0]*1e6:.0f} ppm T14 {th[1]*24:.3f} h flat {th[2]:.2f} | sec depth {th[3]*1e6:.0f} ppm T14 {th[4]*24:.3f} h "
          f"flat {th[5]:.2f} offset {th[6]:+.2f} min | D_flat = chi2(V sec) - chi2(free) = {real['D_flat']:+.2f}; profile {real['prof']}")
    by = {}
    for r in res[1:]:
        tn = r['label'].split('|')[0]; by.setdefault(tn, []).append(r)
    for tn, rr in by.items():
        D = np.array([r['D_flat'] for r in rr]); fsv = np.array([r['th_free'][5] for r in rr])
        print(f"  truth {tn:15s}: D_flat median {np.median(D):+.2f} range [{D.min():+.2f}, {D.max():+.2f}]; frac D >= real: {np.mean(D >= real['D_flat']):.2f};"
              f" recovered f_sec median {np.median(fsv):.2f} [{fsv.min():.2f}, {fsv.max():.2f}]")
    print('time %.0f s' % (time.time() - t0))
