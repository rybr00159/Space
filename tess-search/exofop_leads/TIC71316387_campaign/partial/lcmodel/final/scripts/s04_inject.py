"""Injection-recovery calibration of the geometric EB fit with REAL noise.
For each truth configuration and each phase offset phi (13 offsets, 0.07-0.43 in steps of 0.03 => injected eclipse cores
never overlap the real ones, and successive realisations put the eclipse cores on different stretches of real data):
  1. take the default light curves, delete the real eclipses (|dt| < 0.6 T14 around both conjunctions),
  2. multiply in the synthetic EB (cadence-integrated) with conjunction at T0 + phi P,
  3. build eclipse windows exactly as for the real data and run the identical ML fits (3 total-branch + 4 partial-branch
     starts, none of them at the truth),
  4. record chi2(best partial) - chi2(best total) and the recovered parameters.
This gives the empirical distribution of the branch statistic under each truth, i.e. a calibrated significance for the
value seen in the real data, and the empirical bias/scatter of L3, k, b1, r1, J.
Usage: s04_inject.py <truths.json> <out tag>"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import fit3, fitlib as F, lclib as L, geofit as G
from ebmodel import eb_flux
from multiprocessing import get_context
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/results'
TRUTHS = json.load(open(sys.argv[1])); TAG = sys.argv[2]
PHIS = [round(0.07 + 0.03 * i, 2) for i in range(13)]
CACHE = f'{R}/inj_cache'; os.makedirs(CACHE, exist_ok=True)

def make(tname, phi):
    path = f'{CACHE}/{TAG}_{tname}_{phi:.2f}.pkl'
    if os.path.exists(path): return pickle.load(open(path, 'rb'))
    tr = TRUTHS[tname]; T0i = L.T00 + phi * L.P0
    u1a, u2a = fit3.q2u(tr['q1a'], tr['q2a']); u1b, u2b = fit3.q2u(tr['q1b'], tr['q2b'])
    lcs = []
    for lc in L.default_set():
        t = lc['t']; ph = ((t - L.T00) / L.P0) % 1
        dd = np.minimum(np.minimum(ph, 1 - ph), np.abs(ph - 0.5)) * L.P0
        keep = dd > 0.6 * L.T14
        lc = {k: (v[keep] if isinstance(v, np.ndarray) and v.size == t.size else v) for k, v in lc.items()}
        cosi = tr['bn'] * tr['rsum']
        m = eb_flux(lc['t'], T0i, L.P0, tr['rsum'], tr['k'], cosi, tr['J'], tr['L3'], tr['ecosw'], tr['esinw'],
                    u1a, u2a, u1b, u2b, lc['texp'], fit3.NSUB[lc['sector']], fit3.XG, fit3.WG)
        lc['f'] = lc['f'] * m
        lcs.append(lc)
    d = F.prepare(lcs, hw=0.20, P=L.P0, T0=T0i, verbose=False)
    pickle.dump(d, open(path, 'wb'))
    return d

_dc = {}
def task(a):
    tname, phi, branch, si = a
    key = (tname, phi)
    if key not in _dc:
        _dc.clear()
        d = make(tname, phi)
        d['groups'] = [(sec, idx, texp, fit3.NSUB[sec]) for (sec, idx, texp, nsub) in d['groups']]
        _dc[key] = d
    d = _dc[key]
    st = (G.STARTS_TOTAL if branch == 'total' else G.STARTS_PARTIAL)[si]
    o = G.solve(d, st, branch=branch, T0ref=L.T00 + phi * L.P0)
    o.update(truth=tname, phi=phi, start=si)
    return o

if __name__ == '__main__':
    t0 = time.time()
    # build caches serially first (cheap, avoids races)
    for tn in TRUTHS:
        for phi in PHIS: make(tn, phi)
    print('data built %.0f s' % (time.time() - t0), flush=True)
    jobs = [(tn, phi, b, si) for tn in TRUTHS for phi in PHIS for b, si in
            [('total', 0), ('total', 1), ('total', 2), ('partial', 0), ('partial', 1), ('partial', 2), ('partial', 3)]]
    with get_context('spawn').Pool(int(os.environ.get('NPROC', '4'))) as pool:
        res = pool.map(task, jobs, chunksize=7)
    json.dump(dict(truths=TRUTHS, phis=PHIS, fits=res), open(f'{R}/inject_{TAG}.json', 'w'), indent=1)
    for tn, tr in TRUTHS.items():
        print(f'=== truth {tn}: ' + ' '.join(f'{k}={tr[k]:.4g}' for k in ('rsum', 'k', 'bn', 'J', 'L3', 'esinw')))
        dl = []
        for phi in PHIS:
            rr = [r for r in res if r['truth'] == tn and r['phi'] == phi]
            bt = min([r for r in rr if r['branch'] == 'total'], key=lambda r: r['tot'])
            bp = min([r for r in rr if r['branch'] == 'partial'], key=lambda r: r['tot'])
            b = bt if bt['tot'] <= bp['tot'] else bp
            dl.append(bp['tot'] - bt['tot'])
            print(f"  phi={phi:.2f}: chi2(partial)-chi2(total) = {bp['tot']-bt['tot']:+7.2f} | best: {b['branch']:7s} L3={b['L3']:.3f} k={b['k']:.3f} "
                  f"b1={b['b1']:.2f} r1={b['r1']:.4f} J={b['J']:.3f} esinw={b['esinw']:+.3f}", flush=True)
        dl = np.array(dl)
        print(f'  branch statistic: median {np.median(dl):+.2f}, min {dl.min():+.2f}, max {dl.max():+.2f}')
    print('time %.0f s' % (time.time() - t0))
