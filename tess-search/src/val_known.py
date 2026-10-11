import sys, time, numpy as np, pandas as pd
from fetch_qlp import load
from bls import search_robust as search, detrend, diagnostics
v = pd.read_csv('../meta/val_tois.csv').sort_values('Depth (ppm)').head(int(sys.argv[1]))
rows = []
for _, r in v.iterrows():
    x = load(int(r['TIC ID']), 104)
    if x is None: print('no data', r['TIC ID']); continue
    t, f = detrend(x[0], x[1])
    t1 = time.time(); res, (tc, fc) = search(t, f); d = res; dt = time.time() - t1

    Pt = r['Period (days)']; ratio = res['P'] / Pt
    rows.append(dict(toi=r['TOI'], P_toi=round(Pt, 4), P_bls=round(res['P'], 4), ratio=round(ratio, 3), depth_toi=int(r['Depth (ppm)']),
                     depth_bls=int(d['depth'] * 1e6), dur_toi=round(r['Duration (hours)'], 2), dur_bls=round(res['dur'] * 24, 2),
                     snr_w=round(res['snr_w'], 1), snr_red=round(d['snr_red'], 1), sde=round(res['sde'], 1), n_tr=d['n_tr'], sec=round(d['n_tr']), secsig=round(d['sec_sig'], 1),
                     oe=round(d['oddeven_sig'], 1), maxev=round(d['max_event_frac'], 2), sec_t=round(dt, 2)))
print(pd.DataFrame(rows).to_string())
