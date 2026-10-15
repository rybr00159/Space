# Synthesis s02: (1) separation posterior from TESS localisation map x Gaia IPD window x VHS non-detection,
# (2) position-likelihood factor for the odds (replacing the hard r<=1.24" window, applied once),
# (3) light-fraction / distance likelihood factor from the referee's independent L(f) profiles (r09),
# (4) posterior odds bound : chance-aligned, (5) outer-orbit Monte Carlo predictions.
import numpy as np, json
from scipy.interpolate import RegularGridInterpolator
from astropy.table import Table
rng = np.random.default_rng(20261015)
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
d_pc = 1000 / 3.3622
out = {}
# ---------------- position likelihood pieces ----------------
z = np.load(B + 'localise/v3/results/s05_baseline_map.npz'); g = z['g']; M = z['M']   # M[i=dN, j=dE] = profile dchi2
scale = 6.18 / 6.5417            # map the injection-calibrated 95.45% threshold onto chi2(2 dof)
tess_interp = RegularGridInterpolator((g, g), M, bounds_error=False, fill_value=M.max())
def L_tess(dE, dN):
    return np.exp(-0.5 * scale * tess_interp(np.c_[dN, dE]))
def L_ipd(s, mode='soft'):
    s = np.asarray(s, float)
    if mode == 'hard':
        return ((s >= 0.48) & (s <= 1.2)).astype(float)
    if mode == 'none':
        return np.ones_like(s)
    # soft: ramps 0.35-0.50" and 1.20-1.45" (calibrators: 11 pairs at 0.48-1.43", blended-fraction 45-66% at 0.4-1.2")
    return np.clip((s - 0.35) / 0.15, 0, 1) * np.clip((1.45 - s) / 0.25, 0, 1)
# VHS Ks recovery (3 Ks pawprints; inner h04 + outer h01), mean over frames (frames treated as fully correlated)
h4 = json.load(open(B + 'localise/v3/results/h04_vhs_inner.json'))
seps = sorted({r['sep'] for r in h4[0]['inj']}); dms = sorted({r['dm'] for r in h4[0]['inj']})
R = np.zeros((len(seps), len(dms)))
for fr in h4:
    for r in fr['inj']:
        R[seps.index(r['sep']), dms.index(r['dm'])] += r['frac'] / len(h4)
seps = np.array(seps); dms = np.array(dms)
# beyond 1.2": recovery 1 for dm<=3.0 (h01: 100% at >=1.0" for dm<=3.0)
seps_x = np.r_[seps, 1.5, 3.0]; R_x = np.vstack([R, np.tile(np.where(dms <= 3.0, 1.0, R[-1]), (2, 1))])
# below 0.4": core masked (r~0.7"), treat recovery as falling to 0 at 0.3"
seps_x = np.r_[0.3, seps_x]; R_x = np.vstack([np.zeros(len(dms)), R_x])
rec_interp = RegularGridInterpolator((seps_x, dms), R_x, bounds_error=False, fill_value=None)
def L_vhs(s, dK):
    s = np.atleast_1d(np.clip(np.asarray(s, float), 0.3, 3.0)); dK = np.clip(np.asarray(dK, float), 2.0, 4.0) * np.ones_like(s)
    return np.clip(1 - rec_interp(np.c_[s, dK]), 0, 1)
out['vhs_miss_prob_dK2.6'] = {str(s): round(float(L_vhs(s, 2.6)[0]), 3) for s in (0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2)}
out['vhs_miss_prob_dK2.8'] = {str(s): round(float(L_vhs(s, 2.8)[0]), 3) for s in (0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2)}

def L_pos(dE, dN, dK, ipd='soft', use_vhs=True):
    s = np.hypot(dE, dN)
    L = L_tess(dE, dN) * L_ipd(s, ipd)
    if use_vhs:
        L = L * L_vhs(s, dK)
    return L

# ---------------- bound companions: Monte Carlo of outer orbits ----------------
N = 2_000_000
M_G, M_EB = 1.03, 0.66
Mt = M_G + M_EB
logP = rng.normal(5.03, 2.28, N)                       # Raghavan+2010, days
P_yr = 10**logP / 365.25
a_au = (Mt * P_yr**2)**(1 / 3)
e = rng.uniform(0, 0.8, N)
Mmean = rng.uniform(0, 2 * np.pi, N)
E = Mmean.copy()
for _ in range(30):
    E = E - (E - e * np.sin(E) - Mmean) / (1 - e * np.cos(E))
nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
r = a_au * (1 - e * np.cos(E))
cosi = rng.uniform(-1, 1, N); inc = np.arccos(cosi)
w = rng.uniform(0, 2 * np.pi, N); Om = rng.uniform(0, 2 * np.pi, N)
# position in orbital plane -> sky (X=north?, use standard Thiele-Innes: x toward N, y toward E, z along LOS)
u = w + nu
X = r * (np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(inc))
Y = r * (np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(inc))
Z = r * np.sin(u) * np.sin(inc)
# velocities (au/yr), relative orbit
n = 2 * np.pi / P_yr
h = n * a_au**2 * np.sqrt(1 - e**2)
p = a_au * (1 - e**2)
rdot = h / p * e * np.sin(nu); rf = h / r
vx = rdot * (np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(inc)) + rf * (-np.cos(Om) * np.sin(u) - np.sin(Om) * np.cos(u) * np.cos(inc))
vy = rdot * (np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(inc)) + rf * (-np.sin(Om) * np.sin(u) + np.cos(Om) * np.cos(u) * np.cos(inc))
vz = rdot * np.sin(u) * np.sin(inc) + rf * np.cos(u) * np.sin(inc)
AU_YR_KMS = 4.740470
# relative accelerations (au/yr^2), G M / r^2 along -r
GM = 4 * np.pi**2 * Mt
acc = GM / r**2
az_rel = -acc * Z / r
s_arc = np.hypot(X, Y) / d_pc
dN_ = X / d_pc; dE_ = Y / d_pc
dK = rng.normal(2.62, 0.18, N)       # predicted dKs of the pair (lcmodel 2.62-2.80, referee 2.57-2.77, s01 2.44-2.75)
sel = s_arc < 3.0
Wb = np.zeros(N)
Wb[sel] = L_pos(dE_[sel], dN_[sel], dK[sel])
Wb_hard = np.zeros(N); Wb_hard[sel] = L_pos(dE_[sel], dN_[sel], dK[sel], ipd='hard')
Wb_noipd = np.zeros(N); Wb_noipd[sel] = L_pos(dE_[sel], dN_[sel], dK[sel], ipd='none')
Wb_novhs = np.zeros(N); Wb_novhs[sel] = L_pos(dE_[sel], dN_[sel], dK[sel], use_vhs=False)
EL_b = {k: float(v.mean()) for k, v in dict(soft=Wb, hard=Wb_hard, noipd=Wb_noipd, novhs=Wb_novhs).items()}
f_sep_window = float(((s_arc >= 0.1) & (s_arc <= 1.24)).mean())
out['bound_MC'] = dict(N=N, f_sep_0p1_1p24_MC=f_sep_window, E_Lpos=EL_b)

# ---------------- background: uniform on sky ----------------
xx = np.linspace(-3, 3, 601); XX, YY = np.meshgrid(xx, xx); dA = (xx[1] - xx[0])**2
Ls = {}
for mode in ('soft', 'hard', 'none'):
    acc_ = 0.0
    for dKv, wv in zip(np.linspace(2.26, 2.98, 9), np.exp(-0.5 * ((np.linspace(2.26, 2.98, 9) - 2.62) / 0.18)**2)):
        acc_ += wv * (L_pos(XX.ravel(), YY.ravel(), dKv, ipd=mode).sum() * dA)
    Ls[mode] = acc_ / np.exp(-0.5 * ((np.linspace(2.26, 2.98, 9) - 2.62) / 0.18)**2).sum()
Ls['novhs'] = float(L_pos(XX.ravel(), YY.ravel(), 2.6, use_vhs=False).sum() * dA)
out['background_int_Lpos_dA_arcsec2'] = Ls
circ = np.pi * 1.24**2
# position factor relative to the referee's hard window (0.1-1.24" bound; 1.24" circle background)
fsep_ref = 0.2716531478730537
PF = {k: (EL_b[k] / fsep_ref) / (Ls['none' if k == 'noipd' else k] / circ) for k in ('soft', 'hard', 'noipd', 'novhs')}
out['position_factor_vs_hard_window'] = PF

# ---------------- light-fraction / distance factor from referee L(f) ----------------
r9 = json.load(open(B + 'referee_all/r09_msprior.json'))
G = Table.read(B + 'astro/data/gaia_counts_30arcmin.ecsv')
sep = np.array(G['sep_arcmin']) if 'sep_arcmin' in G.colnames else None
ra0, de0 = 330.887948, -42.835511
if sep is None:
    sep = 60 * np.degrees(np.arccos(np.clip(np.sin(np.radians(de0)) * np.sin(np.radians(np.array(G['dec']))) +
          np.cos(np.radians(de0)) * np.cos(np.radians(np.array(G['dec']))) * np.cos(np.radians(np.array(G['ra']) - ra0)), -1, 1)))
br = np.array(G['bp_rp'].filled(np.nan) if hasattr(G['bp_rp'], 'filled') else G['bp_rp'], float)
Gm = np.array(G['phot_g_mean_mag'], float)
Tm = Gm - 0.00522555 * br**3 + 0.0891337 * br**2 - 0.633923 * br + 0.0324473
area_as2 = np.pi * (30**2 - 1**2) * 3600
host = (sep > 1) & (br >= 2.1) & (br <= 2.5)
f_host = 10**(-0.4 * (Tm[host] - 11.358))
out['n_hosts_BPRP2.1-2.5_f1-15pct'] = int(((f_host > 0.01) & (f_host < 0.15)).sum())
out['n_hosts_window_2.2-6.5'] = int(((f_host >= 0.022) & (f_host <= 0.065)).sum())
res_d = {}
for setname in ('A_Mann', 'B_PM_table'):
    tab = r9[setname]
    ff = np.array([t['f'] for t in tab]); tot = np.array([t['total'] for t in tab]); dist = np.array([t['dist_pc'] for t in tab])
    tmin = tot.min()
    Lf = lambda f: np.exp(-0.5 * (np.interp(np.log(f), np.log(ff), tot, left=1e3, right=1e3) - tmin))
    # background: hosts' f values, effective number relative to hard window count
    Neff = float(Lf(f_host[(f_host > 0.005) & (f_host < 0.3)]).sum())
    # smooth version: host counts per ln f from power-law fit to BP-RP 2.1-2.5 counts in T 13-17
    # bound: f at which MS distance = 297.4 pc, with 7% (0.15 mag) MS-relation scatter
    lnf_297 = np.interp(np.log(d_pc), np.log(dist[::-1]), np.log(ff[::-1]))
    fs = np.exp(rng.normal(lnf_297, 0.15 / 1.0857, 200000))
    Lb = float(Lf(fs).mean())
    res_d[setname] = dict(f_at_297pc=float(np.exp(lnf_297)), L_bound=Lb, Neff_hosts=Neff,
                          n_window=out['n_hosts_window_2.2-6.5'],
                          distance_factor=(Lb / (Neff / out['n_hosts_window_2.2-6.5'])) if Neff > 0 else None,
                          total_min=float(tmin))
out['distance_factor'] = res_d

# ---------------- odds ----------------
r6 = json.load(open(B + 'referee_all/r06_priors.json'))
O_like = r6['odds_like_for_like_Mpair']
out['odds'] = {}
for k in ('soft', 'hard'):
    for setname in res_d:
        out['odds'][f'like_for_like x pos({k}) x dist({setname})'] = O_like * PF[k] * res_d[setname]['distance_factor']
    out['odds'][f'like_for_like x pos({k})'] = O_like * PF[k]
out['odds']['any_host_f0.028 x pos(soft)'] = r6['odds_any_f0.028'] * PF['soft']
# Poisson range on n_match = 7 (68%, 95% intervals of the mean, Garwood)
from scipy.stats import chi2
n = 7
lo68, hi68 = chi2.ppf(0.16, 2 * n) / 2, chi2.ppf(0.84, 2 * (n + 1)) / 2
lo95, hi95 = chi2.ppf(0.025, 2 * n) / 2, chi2.ppf(0.975, 2 * (n + 1)) / 2
out['poisson_n7'] = dict(lo68=lo68, hi68=hi68, lo95=lo95, hi95=hi95, odds_mult68=[n / hi68, n / lo68], odds_mult95=[n / hi95, n / lo95])

# ---------------- separation posterior (bound) ----------------
def wq(x, w, qs):
    o = np.argsort(x); x = x[o]; c = np.cumsum(w[o]); c /= c[-1]
    return [float(np.interp(q, c, x)) for q in qs]
qs = [0.025, 0.16, 0.5, 0.84, 0.975]
m = Wb > 0
out['sep_posterior_arcsec'] = {k: wq(s_arc[m2], W[m2], qs) for k, W, m2 in
                               (('soft', Wb, Wb > 0), ('hard', Wb_hard, Wb_hard > 0), ('noipd', Wb_noipd, Wb_noipd > 0), ('novhs', Wb_novhs, Wb_novhs > 0))}
PA = np.degrees(np.arctan2(dE_, dN_)) % 360
# PA posterior: fraction within NNW sector 300-20 deg vs prior 80/360
inN = (PA >= 300) | (PA <= 20)
out['PA_frac_300_20'] = dict(posterior=float(Wb[inN].sum() / Wb.sum()), prior=80 / 360)
out['proj_sep_au'] = [x * d_pc for x in out['sep_posterior_arcsec']['soft']]
# ---------------- outer orbit predictions (weighted by Wb) ----------------
vrel_los = vz * AU_YR_KMS                         # km/s, EB minus G star
vrel_sky = np.hypot(vx, vy) * AU_YR_KMS
mu_rel = np.hypot(vx, vy) / d_pc * 1000            # mas/yr (au/yr / pc -> arcsec/yr)
a_los_EB = (M_G / Mt) * az_rel                      # au/yr^2, EB's LOS acceleration (relative accel x M_G/Mt)
a_los_G = -(M_EB / Mt) * az_rel
AU_YR2_MS2 = 1.495978707e11 / (3.15576e7)**2
c_light = 2.99792458e8
etv_curv_min_yr2 = np.abs(a_los_EB * AU_YR2_MS2) / (2 * c_light) * (3.15576e7)**2 / 60  # min/yr^2 (O-C = c t^2)
rv_trend_G = np.abs(a_los_G * AU_YR2_MS2) * 3.15576e7                            # m/s/yr
acc_sky_G = (M_EB / Mt) * acc * np.hypot(X, Y) / r                               # au/yr^2 on sky
acc_sky_G_masyr2 = acc_sky_G / d_pc * 1000                                       # mas/yr^2
PA_rate = np.degrees((X * vy - Y * vx) / (X**2 + Y**2))                          # deg/yr
dPdt = np.abs(a_los_EB * AU_YR2_MS2) / c_light * 1.3036723 * 86400 * 3.15576e7  # s/yr
LTTE_semi_amp_min = (M_G / Mt) * a_au * np.sin(inc) * 499.005 / 60               # minutes (semi-amplitude of LTTE over full orbit)
W = Wb
pred = {}
for name, arr in [('a_au', a_au), ('P_out_yr', P_yr), ('e_out', e), ('r_3D_au', r), ('vrel_los_kms', vrel_los), ('|vrel_los|_kms', np.abs(vrel_los)),
                  ('vrel_sky_kms', vrel_sky), ('mu_rel_masyr', mu_rel), ('PA_rate_degyr', np.abs(PA_rate)), ('ETV_curv_min_yr2', etv_curv_min_yr2),
                  ('dPdt_s_per_yr', dPdt), ('RV_trend_G_ms_yr', rv_trend_G), ('acc_sky_G_masyr2', acc_sky_G_masyr2), ('LTTE_semiamp_min', LTTE_semi_amp_min)]:
    pred[name] = wq(arr[m], W[m], qs)
out['outer_orbit_quantiles_2.5_16_50_84_97.5'] = pred
out['P(|vrel_los|<3 km/s)'] = float(W[np.abs(vrel_los) < 3].sum() / W.sum())
out['P(mu_rel<3 mas/yr)'] = float(W[mu_rel < 3].sum() / W.sum())
json.dump(out, open(B + 'synthesis/s02_odds_sep_orbit.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
