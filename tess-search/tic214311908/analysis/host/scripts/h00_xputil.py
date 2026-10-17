"""h00: helper - download Gaia DR3 XP_CONTINUOUS (DataLink, csv, RAW) for a list of source_ids in chunks,
returning one pandas DataFrame with numpy-array coefficient columns (GaiaXPy-ready)."""
import os, glob, zipfile, tempfile, shutil, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')

ARR = ['bp_coefficients', 'bp_coefficient_errors', 'bp_coefficient_correlations', 'rp_coefficients',
       'rp_coefficient_errors', 'rp_coefficient_correlations']

def _parse(s):
    if isinstance(s, str):
        s = s.strip().strip('()[]')
        if not s:
            return np.array([])
        return np.array([float(x) for x in s.split(',')])
    return s

def fetch_xp(ids, chunk=500, tmpdir=None):
    from astroquery.gaia import Gaia
    ids = [int(i) for i in ids]
    out = []
    for k in range(0, len(ids), chunk):
        sub = ids[k:k + chunk]
        d = tempfile.mkdtemp(dir=tmpdir)
        try:
            res = Gaia.load_data(ids=sub, data_release='Gaia DR3', retrieval_type='XP_CONTINUOUS', data_structure='DATAMODEL_GAIA',
                                 format='csv', output_file=os.path.join(d, 'xp.zip'), dump_to_file=False, verbose=False)
            for key, lst in res.items():
                for tab in lst:
                    df = tab.to_pandas()
                    out.append(df)
        finally:
            shutil.rmtree(d, ignore_errors=True)
        print(f'  XP chunk {k // chunk + 1}: {len(sub)} ids', flush=True)
    if not out:
        return pd.DataFrame()
    df = pd.concat(out, ignore_index=True)
    for c in ARR:
        if c in df:
            df[c] = df[c].apply(_parse)
    return df
