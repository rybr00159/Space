"""Collect seed solutions for the profile/variant fits from the multi-start ML results (first_fit_<tag>.json) and
write best_<tag>.json (lowest -log posterior) and seeds_<tag>.json (all distinct solutions plus two hand-made
family members spanning the L3 degeneracy: near-central large-k and grazing small-k).
Usage: c02b_seeds.py <tag> [extra json with {'p':...} ...]"""
import json, sys
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
tag = sys.argv[1]
F = json.load(open(f'{V}/results/first_fit_{tag}.json'))
res = sorted(F['results'], key=lambda r: r['nlp'])
best = res[0]
json.dump(dict(source=f'first_fit_{tag} best', chi2=best['chi2'], nlp=best['nlp'], p=best['p']), open(f'{V}/results/best_{tag}.json', 'w'), indent=1)
seeds = [r['p'] for r in res]
b = dict(best['p'])
for over in (dict(L3=0.972, k=0.45, b0=0.05, rsum=0.161, J=0.52, esinw=0.012), dict(L3=0.5, k=0.116, b0=0.77, rsum=0.233, J=0.39, esinw=-0.006),
             dict(L3=0.0, k=0.083, b0=0.83, rsum=0.265, J=0.38, esinw=-0.005)):
    s = dict(b); s.update(over); seeds.append(s)
for fn in sys.argv[2:]:
    seeds.append(json.load(open(fn))['p'])
json.dump(dict(seeds=seeds), open(f'{V}/results/seeds_{tag}.json', 'w'), indent=1)
print('best', {k: round(v, 5) for k, v in best['p'].items()}, 'chi2', best['chi2'])
print(len(seeds), 'seeds written')
