"""v4: Gaia DR3 queries (archive TAP): target full row (raw VOTable kept, masks preserved), 8' field for TESS modelling."""
import pyvo, json, numpy as np, sys
from astropy.table import Table
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/gaia'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
r = svc.run_sync("SELECT * FROM gaiadr3.gaia_source WHERE source_id = 6571343840265900544").to_table()
r.write(f'{O}/target_row.ecsv', overwrite=True, format='ascii.ecsv')
row = {}
for c in r.colnames:
    v = r[c][0]
    masked = bool(np.ma.is_masked(v))
    row[c] = None if masked else (v.item() if hasattr(v, 'item') else str(v))
json.dump(row, open(f'{O}/target_row.json', 'w'), indent=1, default=str)
keys = [k for k in row if k.startswith(('ipd', 'ruwe', 'astrometric', 'scan', 'visib', 'dupl', 'phot_bp_rp', 'phot_bp_n', 'phot_rp_n', 'matched', 'rv_', 'nu_eff', 'pseudo'))]
for k in keys:
    print(f'{k:40s} {row[k]}')
f = svc.run_sync("""SELECT source_id, ra, dec, ref_epoch, pmra, pmdec, parallax, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp,
 ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, astrometric_params_solved, duplicated_source
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 330.887948, -42.835511, 0.1333))""", maxrec=100000).to_table()
f.write(f'{O}/field_8arcmin.ecsv', overwrite=True, format='ascii.ecsv')
print('field rows', len(f))
