import pyvo
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
cols = "source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,phot_g_mean_mag,bp_rp,phot_rp_mean_mag,ruwe,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,astrometric_params_solved"
q = f"SELECT {cols} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE(330.887948,-42.835511,0.105)) AND phot_g_mean_mag < 20.5"
t = svc.run_sync(q, maxrec=100000).to_table()
t.write('field_6arcmin.csv', overwrite=True)
print(len(t))
