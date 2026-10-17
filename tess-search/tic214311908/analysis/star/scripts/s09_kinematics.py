"""Space motion and young-association membership of TIC 214311908 (and TIC 214311910 A at the M3V's RV).
- UVW (heliocentric, U towards the Galactic centre) from Gaia DR3 astrometry + Gaia DR3 RV (16.28 +- 2.22 km/s),
  Monte Carlo over errors; V_tot relative to the LSR (Schoenrich+10 solar motion).
- Thin/thick-disk/halo relative probabilities (Bensby, Feltzing & Oey 2014, Table A.1 kinematic method).
- BANYAN Sigma (Gagne+18; code from github.com/jgagneastro/banyan_sigma, run locally, unmodified) membership in the
  27 young associations within 150 pc, with and without the RV.
Output: out/kinematics.json"""
import sys, json, warnings
warnings.filterwarnings('ignore')
import numpy as np
B = 'tic214311908/analysis/star'
sys.path.insert(0, f'{B}/data/ext/banyan/banyan_sigma/src')
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
star = dict(ra=260.969674, dec=-51.323821, plx=15.5494, eplx=0.036, pmra=-113.671152, epmra=0.031, pmdec=85.54355, epmdec=0.022, rv=16.283682, erv=2.217367)
rng = np.random.default_rng(3); N = 20000
plx = rng.normal(star['plx'], star['eplx'], N); pmra = rng.normal(star['pmra'], star['epmra'], N)
pmde = rng.normal(star['pmdec'], star['epmdec'], N); rv = rng.normal(star['rv'], star['erv'], N)
c = SkyCoord(ra=np.full(N, star['ra']) * u.deg, dec=np.full(N, star['dec']) * u.deg, distance=(1000 / plx) * u.pc,
             pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmde * u.mas / u.yr, radial_velocity=rv * u.km / u.s)
gal = c.galactic
gal.representation_type = 'cartesian'; gal.differential_type = 'cartesian'
U = gal.velocity.d_x.to(u.km / u.s).value; V = gal.velocity.d_y.to(u.km / u.s).value; W = gal.velocity.d_z.to(u.km / u.s).value
q = lambda x: [float(np.median(x)), float(np.std(x))]
Usun, Vsun, Wsun = 11.1, 12.24, 7.25
Ul, Vl, Wl = U + Usun, V + Vsun, W + Wsun
vtot = np.sqrt(Ul ** 2 + Vl ** 2 + Wl ** 2)
# Bensby+14 (Table A.1): sigma_U, sigma_V, sigma_W, V_asym, X (fraction) for thin D, thick TD, halo H; velocities relative to LSR
pops = {'D': (35, 20, 16, -15, 0.85), 'TD': (67, 38, 35, -46, 0.09), 'H': (160, 90, 90, -220, 0.0015)}
def f(Ul, Vl, Wl, p):
    sU, sV, sW, Va, X = p
    k = 1 / ((2 * np.pi) ** 1.5 * sU * sV * sW)
    return X * k * np.exp(-Ul ** 2 / (2 * sU ** 2) - (Vl - Va) ** 2 / (2 * sV ** 2) - Wl ** 2 / (2 * sW ** 2))
fD, fT, fH = (f(Ul, Vl, Wl, pops[k]) for k in ('D', 'TD', 'H'))
out = dict(U=q(U), V=q(V), W=q(W), U_lsr=q(Ul), V_lsr=q(Vl), W_lsr=q(Wl), vtot_lsr=q(vtot),
           TD_over_D=float(np.median(fT / fD)), P_thin=float(np.median(fD / (fD + fT + fH))), P_thick=float(np.median(fT / (fD + fT + fH))))
print(json.dumps(out, indent=0))
try:
    # numpy 2 removed np.chararray; BANYAN only uses np.chararray.strip -> provide that shim (no change to its logic)
    class _CA: strip = staticmethod(np.char.strip)
    np.chararray = _CA
    from banyan_sigma import core as bs
    res = {}
    for lab, use_rv in (('with_rv', True), ('no_rv', False)):
        kw = dict(ra=np.array([star['ra']]), dec=np.array([star['dec']]), pmra=np.array([star['pmra']]), pmdec=np.array([star['pmdec']]),
                  epmra=np.array([star['epmra']]), epmdec=np.array([star['epmdec']]), plx=np.array([star['plx']]), eplx=np.array([star['eplx']]))
        if use_rv:
            kw.update(rv=np.array([star['rv']]), erv=np.array([star['erv']]))
        r = bs.membership_probability(**kw)
        # r is a pandas DataFrame with multi-index columns; extract the ALL probabilities and best hypothesis
        try:
            allp = r['ALL'].iloc[0].to_dict() if 'ALL' in r.columns.get_level_values(0) else {}
        except Exception:
            allp = {}
        best = r['Global'].iloc[0].to_dict() if 'Global' in r.columns.get_level_values(0) else {}
        res[lab] = dict(global_=best, all=allp)
        print(lab, best)
        print('   ', {k: round(float(v), 4) for k, v in allp.items() if np.isfinite(v) and float(v) > 1e-4})
    out['banyan'] = res
except Exception as e:
    import traceback; traceback.print_exc(); out['banyan'] = f'error {e}'
json.dump(out, open(f'{B}/out/kinematics.json', 'w'), indent=1, default=str)
