"""Summarise the profile likelihoods: real data (L3 and k) and the injection-recovery profiles (L3).
For each profile: minimum, the range of L3 (and of k, b0, rsum, r1, implied f_EB, undiluted primary depth) with
delta chi2 < 1, 4, 9, 25; for injections, delta chi2 at the true L3 (interpolated) and the recovered minimum."""
import numpy as np, sys, json, os
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import ebm, derive
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'

def undil_depth(p):
    dsec = ebm.secondary_offset(p['P'], p['ecosw'], p['esinw'])
    tt = np.array([0.0, dsec])
    f = ebm.flux_grid(tt, p['P'], p['rsum'], p['k'], p['b0'], p['J'], p['ecosw'], p['esinw'], p['u1a'], p['usa'] - p['u1a'], p['u1b'], p['usb'] - p['u1b'], ebm.XG, ebm.WG)
    return 1 - f

out = {}
for tag, which in (('default', 'L3'), ('default', 'k'), ('injA', 'L3'), ('injB', 'L3'), ('injC', 'L3')):
    fn = f'{V}/results/profile_{which}_{tag}.json'
    if not os.path.exists(fn): continue
    pr = sorted(json.load(open(fn))['profile'], key=lambda r: r['val'])
    c0 = min(r['chi2'] for r in pr)
    rows = []
    for r in pr:
        p = r['p']; dp = undil_depth(p); r1 = p['rsum'] / (1 + p['k'])
        rows.append(dict(val=r['val'], dchi2=r['chi2'] - c0, L3=p['L3'], k=p['k'], b0=p['b0'], rsum=p['rsum'], r1=r1, J=p['J'], esinw=p['esinw'],
                         inc=float(np.degrees(np.arccos(p['b0'] * p['rsum']))), undil_prim=float(dp[0]), undil_sec=float(dp[1]),
                         rho1_q0p5=float(derive.rho1(r1, 0.5))))
    res = dict(rows=rows, chi2_min=c0)
    for thr in (1, 4, 9, 25):
        ok = [x for x in rows if x['dchi2'] < thr]
        res[f'range_dchi2_lt_{thr}'] = {k: [min(x[k] for x in ok), max(x[k] for x in ok)] for k in ('L3', 'k', 'b0', 'rsum', 'r1', 'undil_prim', 'rho1_q0p5', 'J')}
    if tag.startswith('inj'):
        tr = json.load(open(f'{V}/results/truth_{tag}.json'))['p']
        xs = np.array([x['L3'] for x in rows]); ys = np.array([x['dchi2'] for x in rows])
        res['truth_L3'] = tr['L3']; res['dchi2_at_truth'] = float(np.interp(tr['L3'], xs, ys))
        res['L3_at_min'] = float(xs[np.argmin(ys)])
        tchi = None
    out[f'{tag}_{which}'] = res
    print(f'== {tag} profile in {which}: chi2_min {c0:.2f}')
    for x in rows:
        print(f"  {which}={x['val']:<6} dchi2 {x['dchi2']:7.2f}  L3 {x['L3']:.4f} k {x['k']:.3f} b0 {x['b0']:.3f} i {x['inc']:.2f} rsum {x['rsum']:.4f} r1 {x['r1']:.4f} J {x['J']:.3f} "
              f"esinw {x['esinw']:+.4f} undiluted depths {x['undil_prim']:.4f}/{x['undil_sec']:.4f} rho1(q=0.5) {x['rho1_q0p5']:.2f}")
    for thr in (4, 9):
        print(f'  delta chi2 < {thr}:', {k: [round(v[0], 4), round(v[1], 4)] for k, v in res[f'range_dchi2_lt_{thr}'].items()})
    if tag.startswith('inj'):
        print(f"  truth L3 {res['truth_L3']:.3f}: delta chi2 at truth {res['dchi2_at_truth']:.2f}; profile minimum at L3 = {res['L3_at_min']:.3f}")
json.dump(out, open(f'{V}/results/profile_summary.json', 'w'), indent=1)
