"""r06b: consistency check of the per-CCD AL error implied by each of the five catalogue astrometric errors, using
the r06 nominal-scanning-law geometry (best nu0 values). Output: out/sigma_ccd_check.json"""
import sys, json, numpy as np
src = open('tic214311908/analysis/referee2_shape/scripts/r06_ruwe.py').read()
src = src[:src.index('obs = {k:')]
exec(src)
from astropy.time import Time
g = gaia
cat = np.array([g['ra_error'], g['dec_error'], g['pmra_error'], g['pmdec_error'], g['parallax_error']])
res = []
for nu0, sense in [(60, 1), (240, 1), (70, 1), (260, 1), (110, -1)]:
    times, psi, pf = simulate(np.radians(nu0), sense)
    for rep in range(20):
        sel = np.sort(rng.choice(len(times), size=min(37, len(times)), replace=False))
        t = times[sel]; ps = psi[sel]; p_ = pf[sel]
        ty = (t - Time('2016-01-01').jd)/365.25
        A = np.column_stack([np.sin(ps), np.cos(ps), np.sin(ps)*ty, np.cos(ps)*ty, p_])
        Cv = np.linalg.inv(A.T @ A*328/len(t))
        e = np.sqrt(np.diag(Cv))       # errors for sigma_f = 1 mas
        res.append(cat/e)
res = np.array(res)
out = dict(params=['ra', 'dec', 'pmra', 'pmdec', 'parallax'], catalogue=cat.tolist(), sigma_f_implied_median=np.median(res, 0).tolist(),
           sigma_f_implied_p16=np.percentile(res, 16, 0).tolist(), sigma_f_implied_p84=np.percentile(res, 84, 0).tolist())
print(json.dumps(out, indent=1))
json.dump(out, open(f'{RD}/out/sigma_ccd_check.json', 'w'), indent=1)
