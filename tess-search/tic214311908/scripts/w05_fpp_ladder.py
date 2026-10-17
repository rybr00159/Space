"""w05: follow-up FPP ladder recomputed with the ADOPTED false-positive composition (writeup/fpp_ladder_adopted.json).

Why: the ladder in referee3_followup/out/r05_hri_ladder.json starts every step after 'now' from a brown-dwarf (BD) posterior
computed with f_BD = 1e-3 and the uncalibrated r03 shape grid (its composition_now BD = 1.33e-4 REF / 5.12e-4 INV), while
the adopted FPP (referee3_bd) uses f_BD = 1.3e-4 and the Ferreira+26 short-period BD eccentricity law. This script keeps
referee3_followup's bookkeeping rules and factors (r05_hri_ladder.py: speckle removes the fraction x of PTP and HEB and
multiplies the BD term by the Zorro BD factor; a two-band transit multiplies PTP and HEB by the colour Bayes factor and
removes NTP (pair A resolved) and WD (chromatic); 'BD closed' sets BD = 0) and only swaps in the adopted composition:
PTP, NTP from referee3_bd/out/s12_summary.json (Kipping planet e, short-period BD e law, f_BD 1.3e-4); WD, HEB, BG
from referee3_bd/out/s02_fpp_fixed.csv (same row). Three BD versions: without the Gaia RV factor (the one the adopted
total FPP uses), with the median Gaia RV factor (the README table), and with the RV factor and the injection-calibrated
shape likelihood ratio. As in the referee's ladder, the speckle-only step runs over the three Zorro contrast curves
(analyst nominal, conservative B and C) and the combined steps use conservative curve B; the two-band steps run over the
referee's colour models and depth-ratio precisions (empirical g/z, empirical r/z, BT-Settl r/z; 5% and 10%). Ranges are
given per BD version ('ranges') and over all three ('ranges_all_bd'). The background term (3e-5) is carried over
unchanged, as in the referee's ladder.
Cross-check: with the referee's own composition_now the script reproduces the published steps to < 1e-12."""
import json, os
W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the write-up folder (writeup/ or tic214311908/)
B = os.path.dirname(W) if os.path.isdir(os.path.join(os.path.dirname(W), 'clean')) else os.path.join(W, 'analysis')
import pandas as pd

L = json.load(open(f'{B}/referee3_followup/out/r05_hri_ladder.json'))
S12 = json.load(open(f'{B}/referee3_bd/out/s12_summary.json'))
S2 = pd.read_csv(f'{B}/referee3_bd/out/s02_fpp_fixed.csv')
ZORRO = {'nominal': 'Zorro832 analyst nominal', 'consB': 'Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")',
         'consC': 'Zorro832 cons C (dI=1.2 at 0.05")'}
COLOURS = ('empirical g/z', 'empirical r/z', 'BT-Settl r/z')
tot = lambda c: float(sum(c.values()))


def ladder(comp, ds):
    hri = L['hri'][f'{ds}__signfixed']; bdf = L['ladder'][ds]['bd_factor_from_zorro']; tb = L['twoband'][f'{ds}__signfixed']
    st = {'now': tot(comp)}
    for zk, zn in ZORRO.items():
        x = hri[zn]; c = dict(comp); c['PTP'] *= 1 - x; c['HEB'] *= 1 - x; c['BD'] *= bdf
        st[f'speckle {zk}'] = tot(c)
    for lab in COLOURS:
        for s in ('0.05', '0.1'):
            bf = tb[lab][f'BF_s{s}']
            c = dict(comp); c['PTP'] *= bf; c['HEB'] *= bf; c['NTP'] = 0; c['WD'] = 0
            st[f'two-band {lab} {s}'] = tot(c)
            for zk, zn in ZORRO.items():
                x = hri[zn]; c2 = dict(c); c2['PTP'] *= 1 - x; c2['HEB'] *= 1 - x; c2['BD'] *= bdf
                st[f'speckle {zk} + two-band {lab} {s}'] = tot(c2)
                c3 = dict(c2); c3['BD'] = 0
                st[f'speckle {zk} + two-band {lab} {s} + BD closed'] = tot(c3)
    return st


STEPS = ('now', 'speckle', 'two-band only', 'speckle + two-band', 'speckle + two-band + BD closed')
def ranges(st):
    g = lambda f: [min(v for k, v in st.items() if f(k)), max(v for k, v in st.items() if f(k))]
    return {'now': g(lambda k: k == 'now'),
            'speckle': g(lambda k: k.startswith('speckle') and 'two-band' not in k),
            'two-band only': g(lambda k: k.startswith('two-band')),
            'speckle + two-band': g(lambda k: k.startswith('speckle consB + two-band') and 'BD closed' not in k),
            'speckle + two-band + BD closed': g(lambda k: k.startswith('speckle consB + two-band') and k.endswith('BD closed'))}


out = dict(method=__doc__, ladders={}, ranges={}, ranges_all_bd={}, published_check={})
for ds in ('ref_beta', 'inv_beta'):
    t = [r for r in S12['table'] if r['ds'] == ds and r['planet_e'] == 'kipping' and r['bd_e'] == 'bdsp'][0]
    a = S2[(S2.ds == ds) & (S2.ecc == 'kipping+bdsp') & (S2.f_BD == 1.3e-4) & (S2.sign == 'fixed')].iloc[0]
    base = dict(PTP=t['PTP'], NTP=t['NTP'], WD=float(a.WD), HEB=float(a.HEB), BG=float(a.BG))
    comps = {'BD without Gaia RV factor (as in the adopted total FPP)': dict(base, BD=t['BD_total_noRV']),
             'BD with the median Gaia RV factor (README table)': dict(base, BD=t['BD_total_RVmed']),
             'BD with Gaia RV factor and calibrated shape LR': dict(base, BD=t['BD_total_RVmed_shapecal'])}
    out['ladders'][ds] = {k: dict(composition=c, steps=ladder(c, ds)) for k, c in comps.items()}
    out['ranges'][ds] = {k: ranges(v['steps']) for k, v in out['ladders'][ds].items()}
    out['ranges_all_bd'][ds] = {k: [min(r[k][0] for r in out['ranges'][ds].values()), max(r[k][1] for r in out['ranges'][ds].values())]
                                for k in STEPS}
    # published ladder (f_BD 1e-3 composition) recomputed with the same code: must reproduce referee3_followup
    pub = ladder(dict(L['ladder'][ds]['composition_now']), ds)
    ps = L['ladder'][ds]['steps']
    chk = {'speckle consB': ps['speckle: Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")'],
           'speckle consB + two-band BT-Settl r/z 0.05': ps['Zorro consB + two-band BT-Settl r/z sigma 0.05'],
           'speckle consB + two-band empirical r/z 0.1 + BD closed': ps['Zorro consB + two-band empirical r/z sigma 0.1 + low-res RV (BD closed)']}
    out['published_check'][ds] = {k: dict(published=v, recomputed=pub[k], absdiff=abs(v - pub[k])) for k, v in chk.items()}
    out['ranges'][ds]['published composition (f_BD 1e-3, uncalibrated shape grid)'] = ranges(pub)

json.dump(out, open(f'{W}/fpp_ladder_adopted.json', 'w'), indent=1)
for ds, d in out['ranges'].items():
    for bd, rr in d.items():
        print(ds, '|', bd); [print(f'   {k:32s} {v[0]:.2e} - {v[1]:.2e}') for k, v in rr.items()]
    print(ds, '| all BD versions'); [print(f'   {k:32s} {v[0]:.2e} - {v[1]:.2e}') for k, v in out['ranges_all_bd'][ds].items()]
print('published check max |diff|:', max(c['absdiff'] for d in out['published_check'].values() for c in d.values()))
