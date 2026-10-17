"""r06: Gaia RUWE test of the 'white dwarf eclipsed by the M3V' scenario with a nominal-scanning-law (NSL) model.
1. NSL: solar aspect angle 45 deg, spin 6 h, precession rate dnu/dt = dlam_sun/dt (sqrt(S^2 - cos^2 xi) + cos xi cos nu)/sin xi
   with S = 4.22 (Lindegren+2012). The star is observed by both fields of view (106.5 deg apart, 1.775 h) once per
   spin while |z.u| < sin(0.35 deg). DR3 span 2014-07-25 .. 2017-05-28. The unknown precession phase nu0 is chosen
   to reproduce the observed Gaia DR3 scan-direction statistics of this source (scan_direction_strength_k1..k4,
   scan_direction_mean_k1..k4) and the number of matched transits (37); results are also given marginalised over the
   10 best nu0 values.
2. Per-CCD AL formal error sigma_f from the catalogue parallax_error and the NSL design matrix (5-parameter
   solution, 328 good CCD observations); the effective error entering RUWE is u0*sigma_f with u0 = UWE/RUWE.
3. Edge-on (eclipsing) circular orbit, photocentre semi-major axis alpha = a (q/(1+q) - eps), random node position
   angle and orbital phase; RUWE_sim^2 = RUWE0^2 + <res^2>/(u0 sigma_f)^2 + noise; RUWE0 drawn from the empirical
   distribution of comparison M dwarfs (../referee_star/data/mdw_cmd.csv, G 13.5-15.5, ipd_frac_multi_peak <= 2,
   with and without the RUWE < 1.25 cut).
4. Likelihood ratio L(M_WD) = p(RUWE = 1.076 | WD) / p(RUWE = 1.076 | single) by kernel density, and the
   threshold version P(RUWE <= 1.10 | WD)/P(RUWE <= 1.10 | single) used by ../shape s16; 20000 orbits per mass.
usage: r06_ruwe.py [plx | <effective per-CCD sigma in mas>]; 'plx' derives sigma_f from the catalogue parallax_error.
Output: out/ruwe_nsl_{mode}.json"""
import sys, json, numpy as np, pandas as pd
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic, BarycentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
rng = np.random.default_rng(11)
B = 'tic214311908/analysis'
RD = f'{B}/referee2_shape'
gaia = json.load(open(f'{RD}/out/gaia_astrom.json'))
ra, dec = 260.969674, -51.323821
c = SkyCoord(ra*u.deg, dec*u.deg).transform_to(BarycentricTrueEcliptic(equinox='J2016'))
lam, bet = np.radians(c.lon.deg), np.radians(c.lat.deg)
uvec = np.array([np.cos(bet)*np.cos(lam), np.cos(bet)*np.sin(lam), np.sin(bet)])
eps = np.radians(23.4393)
ncp = np.array([0, np.sin(eps), np.cos(eps)])     # north celestial pole in ecliptic coordinates
eE = np.cross(ncp, uvec); eE /= np.linalg.norm(eE); eN = np.cross(uvec, eE)
t_start, t_end = Time('2014-07-25').jd, Time('2017-05-28').jd
dt = 0.01   # d
tt = np.arange(t_start, t_end, dt)
# Sun ecliptic longitude (low-precision formula)
n = tt - 2451545.0
L0 = np.radians((280.460 + 0.9856474*n) % 360); g = np.radians((357.528 + 0.9856003*n) % 360)
lsun = L0 + np.radians(1.915)*np.sin(g) + np.radians(0.020)*np.sin(2*g)
dlsun = np.gradient(lsun if False else np.unwrap(lsun), dt)
xi = np.radians(45.); S = 4.22
ssun = np.array([np.cos(lsun), np.sin(lsun), np.zeros_like(lsun)])      # Gaia -> Sun
Rgaia = -1.01*ssun                                                        # barycentric position of Gaia (au)
ezz = np.array([0., 0., 1.])
def simulate(nu0, sense=1):
    nu = np.empty(len(tt)); nu[0] = nu0
    for i in range(1, len(tt)):
        nu[i] = nu[i-1] + dt*dlsun[i-1]*(np.sqrt(S**2 - np.cos(xi)**2) + np.cos(xi)*np.cos(nu[i-1]))/np.sin(xi)
    # spin axis: z = cos xi * s_sun + sin xi * (cos nu * e_ecl_pole_perp + sin nu * (s_sun x pole))
    k1 = np.cross(ssun.T, ezz).T     # in ecliptic plane, perpendicular to s_sun
    z = np.cos(xi)*ssun + np.sin(xi)*(np.cos(nu)*ezz[:, None] + sense*np.sin(nu)*k1)
    zu = uvec @ z
    inb = np.abs(zu) < np.sin(np.radians(0.35))
    # transits: in each in-band interval one transit per FoV per 6 h spin, random spin phase per interval
    idx = np.where(inb)[0]
    if len(idx) == 0: return None
    br = np.where(np.diff(idx) > 1)[0] + 1
    times = []
    for seg in np.split(idx, br):
        ta, tb = tt[seg[0]], tt[seg[-1]] + dt
        ph = rng.uniform(0, 0.25)
        for t0 in np.arange(ta + ph, tb, 0.25):
            times.append(t0)
            if t0 + 1.775/24 < tb: times.append(t0 + 1.775/24)
    times = np.array(times)
    j = np.clip(((times - t_start)/dt).astype(int), 0, len(tt) - 1)
    zz = z[:, j]
    d = np.cross(zz.T, uvec); d /= np.linalg.norm(d, axis=1)[:, None]        # scan direction (star motion)
    psi = np.arctan2(d @ eE, d @ eN)                                           # position angle of scan direction (N through E)
    pf = -(Rgaia[:, j].T*d).sum(1)                                              # AL parallax factor (au -> parallax units)
    return times, psi, pf
def scan_stats(psi):
    out = {}
    for k in (1, 2, 3, 4):
        z = np.mean(np.exp(1j*k*psi)); out[k] = (float(np.abs(z)), float(np.degrees(np.angle(z))/k))
    return out
obs = {k: (gaia[f'scan_direction_strength_k{k}'], gaia[f'scan_direction_mean_k{k}']) for k in (1, 2, 3, 4)}
print('observed', obs, flush=True)
rows = []
for sense, nu0 in [(sg, v) for sg in (1, -1) for v in np.radians(np.arange(0, 360, 10))]:
    s = simulate(nu0, sense)
    if s is None: continue
    times, psi, pf = s
    st = scan_stats(psi)
    # distance in strength space (k1..k4) plus mean-angle mismatch for k=2 (the relevant harmonic)
    dS = sum((st[k][0] - obs[k][0])**2 for k in (1, 2, 3, 4))
    dpsi2 = np.degrees(np.angle(np.exp(1j*np.radians(2*(st[2][1] - obs[2][1])))))/2
    # Gaia convention for the scan direction may differ by 180 deg (k1) or by sign; k2 mean is used with both signs
    dpsi2b = np.degrees(np.angle(np.exp(1j*np.radians(2*(-st[2][1] - obs[2][1])))))/2
    rows.append(dict(sense=sense, nu0=float(np.degrees(nu0)), ntr=len(times), dS=float(dS), dpsi2=float(dpsi2), dpsi2_flip=float(dpsi2b), **{f'S{k}': st[k][0] for k in st}, **{f'M{k}': st[k][1] for k in st}))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R['score'] = R.dS/0.01 + (np.minimum(np.abs(R.dpsi2), np.abs(R.dpsi2_flip))/10)**2 + ((R.ntr - 40)/10)**2
R = R.sort_values('score')
best = list(zip(R.iloc[:10].nu0.values, R.iloc[:10].sense.values))
print(R.head(10).to_string(), flush=True)
# ---------------- RUWE simulations
m = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
base = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ipd_frac_multi_peak <= 2)]
r0_cut = base[base.ruwe < 1.25].ruwe.values; r0_all = base[base.ruwe < 1.4].ruwe.values
UWE = np.sqrt(gaia['astrometric_chi2_al']/(gaia['astrometric_n_good_obs_al'] - 5)); u0 = UWE/gaia['ruwe']
NCCD = 9; NTR = int(gaia['astrometric_matched_transits']); Nobs = gaia['astrometric_n_good_obs_al']
P = 5.6915087; M1 = 0.370; plx = gaia['parallax']; RUWE_OBS = gaia['ruwe']
def kde(x, x0, h=0.01):
    return float(np.mean(np.exp(-0.5*((x - x0)/h)**2))/(h*np.sqrt(2*np.pi)))
res = dict(UWE=float(UWE), u0=float(u0), scan_fit=R.head(10).to_dict('records'), masses={})
pR0 = dict(cut=float(np.mean(r0_cut <= 1.10)), all=float(np.mean(r0_all <= 1.10)))
dens0 = dict(cut=kde(r0_cut, RUWE_OBS), all=kde(r0_all, RUWE_OBS))
res['null'] = dict(P_le_1p10=pR0, density_at_obs=dens0, n_cut=int(len(r0_cut)), n_all=int(len(r0_all)))
geoms = []
for nu0, sense in best:
    times, psi, pf = simulate(np.radians(nu0), sense)
    for rep in range(10):
        sel = np.sort(rng.choice(len(times), size=min(NTR, len(times)), replace=False))
        t = times[sel]; ps = psi[sel]; p_ = pf[sel]
        ty = (t - Time('2016-01-01').jd)/365.25
        A = np.column_stack([np.sin(ps), np.cos(ps), np.sin(ps)*ty, np.cos(ps)*ty, p_])
        Ainv = np.linalg.inv(A.T @ A*NCCD)
        sig_f = gaia['parallax_error']/np.sqrt(Ainv[4, 4])*np.sqrt(Nobs/(NCCD*len(t)))  # rescale to the actual number of CCD obs
        geoms.append((t, ps, A, sig_f))
sigf = np.array([g[3] for g in geoms])
res['sigma_f_per_ccd_mas'] = dict(median=float(np.median(sigf)), p16=float(np.percentile(sigf, 16)), p84=float(np.percentile(sigf, 84)))
print('sigma_f per CCD', res['sigma_f_per_ccd_mas'], 'u0', u0, flush=True)
NS = 20000
SIGMODE = sys.argv[1] if len(sys.argv) > 1 else 'plx'
if SIGMODE != 'plx':
    geoms = [(g0, g1, g2, float(SIGMODE)/u0) for g0, g1, g2, g3 in geoms]   # fixed effective per-CCD error u0*sigma_f
res['sigma_mode'] = SIGMODE
for M2 in [0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1.0, 1.2]:
    epsG = 0.021    # WD share of G-band light: ~1.9-2.3% for the 5000-6500 K WD the eclipse needs (r08 XP-based estimate)
    q = M2/M1; a_au = (M1 + M2)**(1/3.)*(P/365.25)**(2/3.)
    alpha = a_au*(q/(1 + q) - epsG)*plx
    out = {}
    for which in ('nsl', 'uniform'):
        ru = np.empty(NS)
        for i in range(NS):
            t, ps, A, sf = geoms[rng.integers(len(geoms))]
            if which == 'uniform':
                ps = rng.uniform(0, 2*np.pi, len(t)); A = A.copy(); A[:, 0] = np.sin(ps); A[:, 1] = np.cos(ps)
                ty = (t - Time('2016-01-01').jd)/365.25; A[:, 2] = np.sin(ps)*ty; A[:, 3] = np.cos(ps)*ty
            Om = rng.uniform(0, 2*np.pi); ph0 = rng.uniform(0, 2*np.pi)
            al = alpha*np.cos(2*np.pi*t/P + ph0)*np.cos(ps - Om)       # edge-on orbit along the node PA Om
            rr = al - A @ np.linalg.lstsq(A, al, rcond=None)[0]
            # 9 CCDs per transit see the same displacement; per-CCD noise is in RUWE0
            extra = np.sum(rr**2)*NCCD/(Nobs - 5)/(u0*sf)**2
            ru[i] = np.sqrt(rng.choice(r0_cut)**2 + extra)
        out[which] = dict(median=float(np.median(ru)), p01=float(np.percentile(ru, 1)), p05=float(np.percentile(ru, 5)),
                          P_le_1p10=float(np.mean(ru <= 1.10)), P_le_obs=float(np.mean(ru <= RUWE_OBS)),
                          L_threshold=float(np.mean(ru <= 1.10)/pR0['cut']), L_density=float(kde(ru, RUWE_OBS)/dens0['cut']))
        out[which]['m2lnL_threshold'] = float(-2*np.log(max(out[which]['L_threshold'], 1/NS/pR0['cut'])))
        out[which]['m2lnL_density'] = float(-2*np.log(max(out[which]['L_density'], 1e-12)))
    res['masses'][str(M2)] = dict(alpha_mas=float(alpha), **out)
    print(M2, round(alpha, 3), out, flush=True)
json.dump(res, open(f'{RD}/out/ruwe_nsl_{SIGMODE}.json', 'w'), indent=1)
