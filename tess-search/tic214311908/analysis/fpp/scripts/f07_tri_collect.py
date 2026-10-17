"""f07: collect the TRICERATOPS runs (out/tri/*.json, written by f06 via run_tri_batch*.sh) into run-to-run statistics.
Groups: A main (inv200, conservative IPD/DECaPS2 curve, neighbours restricted to the 3.1"/5.5" stars), C no contrast
curve, D all TIC neighbours (TRICERATOPS default treatment), E1 REF light curve, E2 600-s binning incl. S39, E3 optimistic
curve, E4 propagated errors, E5 1.5-px aperture with all neighbours, F/G = A/C with circular EB orbits.
Per run it also derives:
* FPP_no_x2P: x2P scenarios removed and the rest renormalised (TRICERATOPS' x2P likelihood fits one eclipse shape
  only and ignores the second eclipse's phase and depth; our odd/even depths and times agree);
* FPP_system: STP (planet on the unresolved companion) counted as planet, like TP/PTP/DTP;
* NFPP_loc: every neighbour's NTP/NEB probability multiplied by its pixel likelihood ratio, then renormalised:
  injection-calibrated Student-t LRs (f01) for pair A (10^-3.22), pair B, the 3.1" (10^-3.87), 5.5" and 12" stars; for the
  other neighbours exp(-Delta chi2/(2 x 4.8)) from the referee's per-star fits, or 3e-7 for unlisted stars > 40" away
  (search agent's whole-cutout scan), capped at pair A's LR (see star_lr).
  NFPP_loc_gauss uses the Gaussian-fit injection LRs instead (pair A 10^-6.7).
Validation thresholds (Giacalone et al. 2021): FPP < 0.015 and NFPP < 1e-3.
Output: out/tri_runs.csv, out/tri_summary.json"""
import json, glob, sys, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
LOC = jload(f'{F}/out/loc_likelihoods.json')
LR = {'1515709327': 10**LOC['catalogued']['G175']['log10LR_t'], '1515709330': 10**LOC['catalogued']['G177']['log10LR_t']}
LR_OTHER = 10**LOC['catalogued']['pairA']['log10LR_t']
LRG = {'1515709327': 10**LOC['catalogued']['G175']['log10LR_gauss'], '1515709330': 10**LOC['catalogued']['G177']['log10LR_gauss']}
LRG_OTHER = 10**LOC['catalogued']['pairA']['log10LR_gauss']
# per-star pixel weights for the all-neighbour runs: injection-calibrated LRs for pair A (TIC 214311910), pair B
# (TIC 1515709322), the 3.1"/5.5" stars and the G17.0 star (TIC 214311904); every other neighbour gets
# exp(-Delta chi2/(2 x 4.8)) from its own hypothesis fit in referee_loc/out/fit_base.json (comb_dchi2_all; kappa 4.8 is
# the most conservative injection temperature), capped at pair A's LR; stars not among the referee's 93 hypotheses
# (beyond 80" or T > 17.6) get pair A's LR. TIC -> Gaia DR3 from data/tri/tic_neighbours.csv (MAST TIC query).
FIT = jload('referee_loc/out/fit_base.json')['real']['comb_dchi2_all']
TICG = pd.read_csv(f'{F}/data/tri/tic_neighbours.csv', dtype={'ID': str, 'GAIA': str}).set_index('ID').GAIA.to_dict()
INJ = {'214311910': 'pairA', '1515709322': 'pairB', '1515709327': 'G175', '1515709330': 'G177', '214311904': 'G170'}
# Stars outside the referee's list but > 40" (2 px) from the M3V: the search agent's whole-cutout scan (search/README.md)
# found every position more than 2 px away worse by Delta chi2 = 24-121 per sector, i.e. >= 144 over six sectors ->
# exp(-144/9.6) = 3e-7.
def star_lr(tic, kind='t', sep=None):
    key = 'log10LR_t' if kind == 't' else 'log10LR_gauss'
    if tic in INJ: return 10**LOC['catalogued'][INJ[tic]][key]
    cap = 10**LOC['catalogued']['pairA'][key]
    g = TICG.get(tic)
    if g is not None and g in FIT: return min(cap, float(np.exp(-FIT[g]/(2*4.8))))
    if sep is not None and sep > 40: return min(cap, float(np.exp(-144/9.6)))
    return cap
PLANET = ('TP', 'PTP', 'DTP')
rows, scen = [], []
for f in sorted(glob.glob(f'{F}/out/tri/*.json')):
    tag = f.split('/')[-1][:-5]
    if '_s' not in tag: continue
    d = json.load(open(f))
    df = pd.DataFrame(d['scenarios']); df['prob'] = df.prob.astype(float); df['ID'] = df.ID.astype(str)
    tgt = df.ID == '214311908'
    p = df.prob.values
    nb = ~tgt
    x2p = df.scenario.str.endswith('x2P').values
    # no-x2P renormalisation
    pn = np.where(x2p, 0, p); pn = pn/pn.sum()
    fpp_nox2p = 1 - pn[(tgt & df.scenario.isin(PLANET)).values].sum()
    fpp_sys = 1 - p[(tgt & df.scenario.isin(PLANET + ('STP',))).values].sum()
    # localisation-weighted neighbours
    SEP = {str(c['ID']): float(c['sep (arcsec)']) for c in d['nfp_candidates']}
    w = np.array([1.0 if i == '214311908' else star_lr(i, 't', SEP.get(i)) for i in df.ID.values])
    pl = p*w; pl = pl/pl.sum()
    nfpp_loc = pl[nb.values].sum()
    wg = np.array([1.0 if i == '214311908' else star_lr(i, 'g', SEP.get(i)) for i in df.ID.values])
    pg = p*wg; pg = pg/pg.sum(); nfpp_loc_g = pg[nb.values].sum()
    fpp_loc = 1 - pl[(tgt & df.scenario.isin(PLANET)).values].sum()
    pln = np.where(x2p, 0, pl); pln = pln/pln.sum()
    fpp_loc_nox2p = 1 - pln[(tgt & df.scenario.isin(PLANET)).values].sum()
    grp = tag.split('_s')[0]
    rows.append(dict(run=tag, group=grp, seed=d['seed'], N=d['N'], lc=d['lc_tag'], cc=d['cc'], nearby=d['nearby'],
                     ap=d['aperture_radius_px'], err=d['errmode'], eb_circular=d.get('eb_circular', False),
                     FPP=d['FPP'], NFPP=d['NFPP'], FPP_no_x2P=fpp_nox2p, FPP_system=fpp_sys,
                     NFPP_loc=nfpp_loc, NFPP_loc_gauss=nfpp_loc_g, FPP_loc=fpp_loc, FPP_loc_no_x2P=fpp_loc_nox2p,
                     n_nfp_candidates=len(d['nfp_candidates']) - 1, runtime_s=d['runtime_s']))
    for _, r in df.iterrows():
        scen.append(dict(group=grp, seed=d['seed'], ID=r.ID, scenario=r.scenario, prob=r.prob))
R = pd.DataFrame(rows); R.to_csv(f'{F}/out/tri_runs.csv', index=False)
S = pd.DataFrame(scen)
out = dict(method=__doc__, thresholds=dict(FPP=0.015, NFPP=0.001), groups={})
desc = dict(A='main: INV 200-s sectors, conservative IPD/DECaPS2 contrast curve, neighbours = 3.1"/5.5" stars only',
            C='as A, no contrast curve', D='as A, all TIC neighbours kept (TRICERATOPS default), N=3e4',
            E1='REF light curve', E2='INV incl. S39, 600-s bins', E3='optimistic contrast curve', E4='propagated (larger) errors',
            E5='1.5-px aperture, all neighbours, N=3e4', F='as A with circular EB orbits', G='as C with circular EB orbits')
for g, sub in R.groupby('group'):
    st = {}
    for c in ('FPP', 'NFPP', 'FPP_no_x2P', 'FPP_system', 'NFPP_loc', 'NFPP_loc_gauss', 'FPP_loc', 'FPP_loc_no_x2P'):
        v = sub[c].values
        st[c] = dict(mean=float(v.mean()), std=float(v.std(ddof=1)) if len(v) > 1 else None, median=float(np.median(v)),
                     min=float(v.min()), max=float(v.max()), sem=float(v.std(ddof=1)/np.sqrt(len(v))) if len(v) > 1 else None)
    ss = S[S.group == g].groupby(['ID', 'scenario']).prob.mean().sort_values(ascending=False)
    st['top_scenarios_mean_prob'] = {f'{i}:{s}': float(v) for (i, s), v in ss.head(10).items()}
    st['n_runs'] = int(len(sub)); st['N'] = int(sub.N.iloc[0]); st['description'] = desc.get(g, '')
    st['frac_runs_FPP_lt_0p015'] = float(np.mean(sub.FPP < 0.015)); st['frac_runs_NFPP_lt_0p001'] = float(np.mean(sub.NFPP < 0.001))
    out['groups'][g] = st
    print(g, desc.get(g, ''), 'n', len(sub), '| FPP %.4f +- %.4f (median %.4f, %.4f-%.4f) | NFPP %.2e | FPP_no_x2P %.4f | FPP_system %.4f | NFPP_loc %.2e (Gauss %.1e)' % (
        st['FPP']['mean'], st['FPP']['std'] or 0, st['FPP']['median'], st['FPP']['min'], st['FPP']['max'], st['NFPP']['mean'],
        st['FPP_no_x2P']['mean'], st['FPP_system']['mean'], st['NFPP_loc']['mean'], st['NFPP_loc_gauss']['mean']))
    print('   top:', {k: round(v, 4) for k, v in list(st['top_scenarios_mean_prob'].items())[:6]})
jdump(out, 'out/tri_summary.json')
