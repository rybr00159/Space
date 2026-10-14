"""Summarise injection-recovery runs: per truth, the empirical distribution of the branch statistic
D = chi2(best partial) - chi2(best total) and of the recovered parameters (overall best fit), compared with the real data.
Usage: s04b_inject_summary.py <tag> [<tag2> ...]"""
import numpy as np, json, sys
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final'
real = json.load(open(f'{FD}/results/ml_default_recent.json'))
rb = {b: min([r for r in real if r['kind'] == 'branch' and r['branch'] == b], key=lambda r: r['tot']) for b in ('total', 'partial')}
Dreal = rb['partial']['tot'] - rb['total']['tot']
print(f'REAL DATA: D = chi2(partial) - chi2(total) = {Dreal:+.2f}; best total L3={rb["total"]["L3"]:.4f} r1={rb["total"]["r1"]:.4f} k={rb["total"]["k"]:.3f} b1={rb["total"]["b1"]:.2f} J={rb["total"]["J"]:.3f}')
out = {'D_real': Dreal}
for tag in sys.argv[1:]:
    R = json.load(open(f'{FD}/results/inject_{tag}.json'))
    for tn, tr in R['truths'].items():
        rows = []
        for phi in R['phis']:
            rr = [r for r in R['fits'] if r['truth'] == tn and r['phi'] == phi]
            bt = min([r for r in rr if r['branch'] == 'total'], key=lambda r: r['tot'])
            bp = min([r for r in rr if r['branch'] == 'partial'], key=lambda r: r['tot'])
            b = bt if bt['tot'] <= bp['tot'] else bp
            rows.append(dict(phi=phi, D=bp['tot'] - bt['tot'], L3=b['L3'], r1=b['r1'], k=b['k'], b1=b['b1'], J=b['J'], branch=b['branch'],
                             L3_tot=bt['L3'], r1_tot=bt['r1'], k_tot=bt['k'], J_tot=bt['J']))
        D = np.array([x['D'] for x in rows])
        r1t = tr['rsum'] / (1 + tr['k']); b1t = tr['bn'] * tr['rsum'] / r1t
        tot_true = tr['bn'] <= abs(1 - tr['k']) / (1 + tr['k'])
        s = dict(truth=tr, true_branch='total' if tot_true else 'partial', D=D.tolist(), D_med=float(np.median(D)), D_min=float(D.min()), D_max=float(D.max()),
                 frac_D_ge_real=float(np.mean(D >= Dreal)), frac_D_le_minus_real=float(np.mean(D <= -Dreal)), rows=rows)
        for key in ('L3', 'r1', 'k', 'J', 'b1'):
            v = np.array([x[key] for x in rows]); s[key] = [float(np.percentile(v, q)) for q in (0, 16, 50, 84, 100)]
        out[f'{tag}:{tn}'] = s
        print(f"=== {tag}:{tn} (true branch {s['true_branch']}; L3={tr['L3']:.3f} r1={r1t:.4f} k={tr['k']:.3f} b1={b1t:.2f} J={tr['J']:.3f})")
        print(f"   D: median {np.median(D):+.2f}, range [{D.min():+.2f}, {D.max():+.2f}]; fraction with D >= real ({Dreal:+.1f}): {np.mean(D >= Dreal):.2f}")
        for key, tv in (('L3', tr['L3']), ('r1', r1t), ('k', tr['k']), ('J', tr['J']), ('b1', b1t)):
            print(f"   recovered {key}: median {s[key][2]:.4f}, 16-84% [{s[key][1]:.4f}, {s[key][3]:.4f}], range [{s[key][0]:.4f}, {s[key][4]:.4f}]  (true {tv:.4f})")
        print('   per realisation (phi: D, branch, L3, k, b1):', ' '.join(f"{x['phi']:.2f}:{x['D']:+.1f},{x['branch'][0]},{x['L3']:.3f},{x['k']:.2f},{x['b1']:.2f}" for x in rows))
json.dump(out, open(f'{FD}/results/inject_summary_{"_".join(sys.argv[1:])}.json', 'w'), indent=1)
