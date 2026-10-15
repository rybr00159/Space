"""Empirical check of the TESScut time system (target-specific Roemer delay, a01_timefix.py) against TESS-SPOC.
In S1, S28 and S68 the same eclipses are timed in SPOC PDCSAP (SPOC target-specific TIMECORR) and in our own TESScut
photometry (FFI TIMECORR replaced by the Horizons target delay).  Paired per-event differences (same epoch and
eclipse type) share most of the photon noise (same pixels), so their mean measures any residual time-system offset.
Error of the mean = sd(differences)/sqrt(N).  Also the size of the correction that was applied, for comparison:
an uncorrected TESScut series would be off by the a01 'add_s' values (-15 to +21 s).
Usage: t01_timecheck.py <etv json default> <etv json tcall>"""
import numpy as np, sys, json
A = json.load(open(sys.argv[1])); Bt = json.load(open(sys.argv[2]))
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
tf = json.load(open(f'{V}/results/timefix.json'))
out = {}
for sec in (1, 28, 68):
    a = {(e['n'], e['ph']): e for e in A['events'] if e['key'] == f'SPOC{sec}' and np.isfinite(e['delta_min'])}
    b = {(e['n'], e['ph']): e for e in Bt['events'] if e['key'] == f'TC{sec}' and np.isfinite(e['delta_min'])}
    common = sorted(set(a) & set(b))
    for lab, sel in (('prim', 0.0), ('sec', 0.5), ('both', None)):
        d = np.array([(b[c]['t_obs'] - a[c]['t_obs']) * 86400 for c in common if sel is None or c[1] == sel])
        if d.size < 3: continue
        out[f'S{sec}_{lab}'] = dict(n=int(d.size), mean_s=float(d.mean()), sem_s=float(d.std(ddof=1) / np.sqrt(d.size)), median_s=float(np.median(d)),
                                    sd_s=float(d.std(ddof=1)))
    o = out[f'S{sec}_both']
    print(f"S{sec}: TC - SPOC eclipse times over {o['n']} paired events: mean {o['mean_s']:+.1f} +- {o['sem_s']:.1f} s (median {o['median_s']:+.1f}, sd {o['sd_s']:.0f} s); "
          f"correction applied to TESScut: {tf[f'tesscut_s{sec}']['add_s_mean']:+.1f} s", flush=True)
m = np.array([out[f'S{s}_both']['mean_s'] for s in (1, 28, 68)]); e = np.array([out[f'S{s}_both']['sem_s'] for s in (1, 28, 68)])
w = 1 / e ** 2
out['combined'] = dict(mean_s=float((w * m).sum() / w.sum()), err_s=float(1 / np.sqrt(w.sum())), chi2=float((((m - (w * m).sum() / w.sum()) / e) ** 2).sum()))
print('combined TC - SPOC offset: %+.1f +- %.1f s (chi2 %.1f / 2)' % (out['combined']['mean_s'], out['combined']['err_s'], out['combined']['chi2']))
json.dump(out, open(f'{V}/results/timecheck.json', 'w'), indent=1)
