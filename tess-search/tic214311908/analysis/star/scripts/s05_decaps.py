"""DECaPS2 (Saydjari et al. 2023) object catalogue around the M3V and the Mira from NOIRLab Data Lab TAP
(Legacy Survey DR10 has no coverage here). Radius 90" around the M3V; plus a 10' x 10' field for star counts.
Outputs: data/decaps2_90.csv, data/decaps2_field10.csv"""
import warnings; warnings.filterwarnings('ignore')
import pyvo
D = 'tic214311908/analysis/star/data'
s = pyvo.dal.TAPService('https://datalab.noirlab.edu/tap')
r = s.search("select * from decaps_dr2.object where 't' = q3c_radial_query(ra, dec, 260.969674, -51.323821, 0.025)").to_table()
r.to_pandas().to_csv(f'{D}/decaps2_90.csv', index=False); print(len(r)); print(r.colnames)
r2 = s.search("select ra, dec, mean_mag_g, mean_mag_r, mean_mag_i, mean_mag_z, mean_mag_y, nmag_ok_g, nmag_ok_r, nmag_ok_i, nmag_ok_z from decaps_dr2.object "
              "where 't' = q3c_radial_query(ra, dec, 260.969674, -51.323821, 0.0833)").to_table()
r2.to_pandas().to_csv(f'{D}/decaps2_field10.csv', index=False); print(len(r2))
