"""q03: NEOWISE-R single-exposure photometry of the Mira (Gaia DR3 5925373268528684288) and of the M3V, read-only IRSA TAP.
Output: data/neowise_mira.csv"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import pandas as pd
from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
D=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/data'
g=pd.read_csv(f'{D}/gaia_cone25.csv'); m=g[g.source_id==5925373268528684288].iloc[0]
print('Mira', m.ra, m.dec)
q=f"""SELECT ra, dec, mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, w1sat, w2sat, qual_frame, cc_flags, nb, na
FROM neowiser_p1bs_psd WHERE CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{m.ra},{m.dec},0.0012))=1"""
t=Irsa.query_tap(q).to_table().to_pandas(); t.to_csv(f'{D}/neowise_mira.csv',index=False); print(len(t)); print(t.describe().T)
