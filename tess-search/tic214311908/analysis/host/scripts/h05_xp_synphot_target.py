"""h05: XP synthetic photometry (GaiaXPy generate, standardised systems JKC_Std, SDSS_Std, PanSTARRS1_Std, DECam)
for the M3V and the A(+B) entry. Output: out/xp_synphot_target.csv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from gaiaxpy import generate, PhotometricSystem
B = 'tic214311908/analysis/host'
t = pd.read_pickle(f'{B}/data/xp/xp_cont_target.pkl')
systems = [PhotometricSystem.JKC_Std, PhotometricSystem.SDSS_Std, PhotometricSystem.PanSTARRS1_Std, PhotometricSystem.DECam]
out = generate(t, photometric_system=systems, save_file=False)
out.to_csv(f'{B}/out/xp_synphot_target.csv', index=False)
pd.set_option('display.width', 250)
cols = [c for c in out.columns if 'mag' in c and 'error' not in c]
print(out[['source_id'] + cols].T)
