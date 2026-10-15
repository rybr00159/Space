"""Rebuild the default (strict cuts + AllWISE) cleaned W1/W2 series exactly as in wise3.py."""
import numpy as np
import lib3 as L
import json
import pandas as pd

T = json.load(open(L.V3 + '/results/tess3.json'))
E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
sod = T['s_over_d']
EX = 7.7 / 86400
DATA = L.CH + '/data/wise/'


def load(cut='strict', allwise=True):
    n = pd.read_csv(DATA + 'neowise.csv')
    if allwise:
        a = pd.read_csv(DATA + 'allwise_mep.csv').rename(columns={'w1mpro_ep': 'w1mpro', 'w1sigmpro_ep': 'w1sigmpro', 'w2mpro_ep': 'w2mpro',
                                                                 'w2sigmpro_ep': 'w2sigmpro'})
        a['qual_frame'] = 10 * (a.qi_fact.values > 0)
        n = pd.concat([n, a], ignore_index=True)
    m = np.isfinite(n.w1mpro.values) & (n.nb.values == 1) & (n.cc_flags.astype(str).isin(['0', '0000']))
    if cut == 'strict':
        m &= (n.qual_frame.values > 0) & (n.saa_sep.values > 5) & (n.moon_masked.astype(str).isin(['0', '00', '0000']))
    elif cut == 'loose':
        m &= (n.qual_frame.values > 0)
    n = n[m]
    t = L.jd_utc_to_bjd_tdb(n.mjd.values + 2400000.5)
    return n, t


def cleaned(band, cut='strict', allwise=True):
    n, t = load(cut, allwise)
    fl = 10 ** (-0.4 * (n[band + 'mpro'].values - 10.3))
    ok = np.isfinite(fl)
    t, fl = t[ok], fl[ok]
    seg = L.segments(t, 5.0)
    med = np.array([np.median(fl[seg == j]) for j in range(seg.max() + 1)])
    f = fl / med[seg]
    r, cl = L.run_survey(t, f, seg, shp, EX, 0, nfake=50, ninj=2, inj_depths=(0.006,), s_over_d=sod, minpts=5, floor=3e-3, seed=1)
    return cl['t'], cl['f'], cl['w'], cl['seg'], r
