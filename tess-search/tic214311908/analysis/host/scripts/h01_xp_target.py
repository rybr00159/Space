"""h01: Retrieve Gaia DR3 XP continuous spectra for the M3V (TIC 214311908) and for the wide companion entry
TIC 214311910 A (whose BP/RP light includes B at 0.83"), via the Gaia archive DataLink, and calibrate them with
GaiaXPy (sampled 336-1020 nm, 2 nm). Output: data/xp/xp_cont_target.pkl, out/xp_target_sampled.csv, out/xp_target_meta.csv"""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
B = 'tic214311908/analysis/host'
sys.path.insert(0, f'{B}/scripts')
from h00_xputil import fetch_xp
ids = [5925373272850784768, 5925373272850782208]
t = fetch_xp(ids)
print(t.shape); print(list(t.columns))
t.to_pickle(f'{B}/data/xp/xp_cont_target.pkl')
meta = ['source_id', 'bp_n_parameters', 'rp_n_parameters', 'bp_chi_squared', 'rp_chi_squared', 'bp_degrees_of_freedom',
        'rp_degrees_of_freedom', 'bp_n_measurements', 'rp_n_measurements', 'bp_n_rejected_measurements', 'rp_n_rejected_measurements',
        'bp_standard_deviation', 'rp_standard_deviation', 'bp_n_relevant_bases', 'rp_n_relevant_bases', 'bp_n_transits', 'rp_n_transits',
        'bp_n_contaminated_transits', 'rp_n_contaminated_transits', 'bp_n_blended_transits', 'rp_n_blended_transits']
meta = [m for m in meta if m in t]
print(t[meta].T)
t[meta].to_csv(f'{B}/out/xp_target_meta.csv', index=False)
from gaiaxpy import calibrate
cal, wl = calibrate(t, sampling=np.arange(336, 1021, 2.0), save_file=False, truncation=False)
rows = []
for _, r in cal.iterrows():
    rows.append(pd.DataFrame(dict(source_id=r.source_id, wl_nm=wl, flux=r.flux, flux_error=r.flux_error)))
pd.concat(rows).to_csv(f'{B}/out/xp_target_sampled.csv', index=False)
print('saved')
