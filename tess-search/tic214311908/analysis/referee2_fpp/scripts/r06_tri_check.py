"""r06: independent re-aggregation of the analyst's TRICERATOPS runs (fpp/out/tri/*.json, read-only) and of any
referee-2 re-runs (referee2_fpp/out/tri/*.json, written by r06b). For each run: FPP, NFPP, the scenario probabilities
summed by class; FPP recomputed from the stored scenario table (1 - TP - PTP - DTP on the target); 'twin-aware' variants.
Also tests whether TRICERATOPS' STP (planet on an unresolved companion, counted as FP by TRICERATOPS) is the analog of
the Bayesian PTP or of PTP + TP-twin. Output: out/r06_tri_check.json"""
import json, glob, numpy as np, pandas as pd, sys, os
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
rows = []
for src, pat in (('analyst', f'{B}/fpp/out/tri/*.json'), ('referee2', f'{R2}/out/tri/*.json')):
    for f in sorted(glob.glob(pat)):
        tag = os.path.basename(f)[:-5]
        try: d = json.load(open(f))
        except Exception: continue
        if 'scenarios' not in d: continue
        df = pd.DataFrame(d['scenarios']); df['prob'] = df.prob.astype(float); df['ID'] = df.ID.astype(str)
        tgt = df.ID == '214311908'
        bysc = df[tgt].groupby('scenario').prob.sum().to_dict()
        planet = sum(bysc.get(s, 0) for s in ('TP', 'PTP', 'DTP'))
        rows.append(dict(src=src, run=tag, group=tag.split('_s')[0], seed=d.get('seed'), N=d.get('N'), FPP_stored=d['FPP'], NFPP_stored=d['NFPP'],
                         FPP_recomputed=1 - planet, NFPP_recomputed=float(df[~tgt].prob.sum()), eb_circular=d.get('eb_circular', False),
                         lc=d.get('lc_tag'), cc=d.get('cc'), nearby=d.get('nearby'), err=d.get('errmode'), flux_err=d.get('flux_err_0'),
                         **{f'P_{k}': float(v) for k, v in bysc.items()}))
R = pd.DataFrame(rows).fillna(0.0)
R.to_csv(f'{R2}/out/r06_tri_runs.csv', index=False)
out = {}
for (src, g), sub in R.groupby(['src', 'group']):
    v = sub.FPP_recomputed.values
    # bootstrap of the seed mean
    rng = np.random.default_rng(1)
    bm = [np.mean(rng.choice(v, len(v))) for _ in range(4000)] if len(v) > 1 else [v[0]]
    out[f'{src}:{g}'] = dict(n=int(len(v)), N=int(sub.N.iloc[0]), FPP_mean=float(v.mean()), FPP_sd=float(v.std(ddof=1)) if len(v) > 1 else None,
                             FPP_median=float(np.median(v)), FPP_min=float(v.min()), FPP_max=float(v.max()),
                             FPP_mean_boot68=[float(np.percentile(bm, 16)), float(np.percentile(bm, 84))],
                             frac_lt_0p015=float(np.mean(v < 0.015)), max_abs_stored_minus_recomputed=float(np.max(np.abs(sub.FPP_stored - sub.FPP_recomputed))),
                             NFPP_mean=float(sub.NFPP_recomputed.mean()),
                             mean_P={c: float(sub[c].mean()) for c in sub.columns if c.startswith('P_')})
    print(src, g, json.dumps(out[f'{src}:{g}'], default=float)[:600])
jdump(out, 'out/r06_tri_check.json')
