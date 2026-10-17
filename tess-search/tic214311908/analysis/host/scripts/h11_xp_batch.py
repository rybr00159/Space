"""h11: batch download of Gaia DR3 XP continuous spectra (DataLink, DATAMODEL_GAIA, csv; anonymous limit 50 ids per request) for a list of source_ids
(data/xp/xp_ids.csv, built by h12a), chunk by chunk; each chunk is calibrated at once with GaiaXPy (sampled flux on
336-1020 nm, 2 nm) and only the sampled flux/error, the 2x55 coefficients and the XP summary numbers are kept
(data/xp/xp_chunk_XXX.npz); the raw download is deleted. Resumable (skips existing chunks).
Usage: python h11_xp_batch.py"""
import os, sys, numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
B = 'tic214311908/analysis/host'
sys.path.insert(0, f'{B}/scripts')
from h00_xputil import fetch_xp
from gaiaxpy import calibrate
ids = pd.read_csv(f'{B}/data/xp/xp_ids.csv').source_id.astype(np.int64).values
CH = 400
wl = np.arange(336, 1021, 2.0)
for c in range(0, len(ids), CH):
    fn = f'{B}/data/xp/xp_chunk_{c // CH:03d}.npz'
    if os.path.exists(fn):
        continue
    sub = ids[c:c + CH]
    for attempt in range(4):
        try:
            t = fetch_xp(sub, chunk=50); break
        except Exception as e:
            print('retry', c, str(e)[:120], flush=True); time.sleep(30); t = None
    if t is None or len(t) == 0:
        print('failed chunk', c, flush=True); continue
    cal, w = calibrate(t, sampling=wl, save_file=False, truncation=False)
    cal = cal.set_index('source_id').loc[t.source_id.values]
    np.savez_compressed(fn, source_id=t.source_id.values.astype(np.int64), flux=np.vstack(cal.flux.values).astype(np.float32),
                        flux_error=np.vstack(cal.flux_error.values).astype(np.float32),
                        bp=np.vstack(t.bp_coefficients.values).astype(np.float32), rp=np.vstack(t.rp_coefficients.values).astype(np.float32),
                        bp_std=t.bp_standard_deviation.values, rp_std=t.rp_standard_deviation.values,
                        bp_nrel=t.bp_n_relevant_bases.values, rp_nrel=t.rp_n_relevant_bases.values,
                        bp_nrej=t.bp_n_rejected_measurements.values, rp_nrej=t.rp_n_rejected_measurements.values,
                        bp_nmeas=t.bp_n_measurements.values, rp_nmeas=t.rp_n_measurements.values)
    print('chunk', c // CH, len(t), flush=True)
print('done')
