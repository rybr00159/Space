"""s10: secondary-eclipse and phase-curve limits (INV full light curve; REF has windows only, its phase-0.5 value
is taken from ../referee_phot/out/phase_scan_clean.json).
1. Secondary eclipse: occultation template = the planet transit shape of the s03 INV free max-likelihood fit with
   uniform disk (u = 0), normalised to unit depth, centred on phase 0.5 (circular) -> linear depth. Empirical error
   from the same fit at 400 random 'phases' (avoiding +-0.06 in phase around transit and secondary); also a scan over
   all phases (eccentric orbits). Depths in the 979 e-/s normalisation.
2. Limits: geometric albedo A_g = depth/(Rp/a)^2, TESS-band brightness temperature (blackbody planet vs
   blackbody 3350 K star and vs a 3350 K star with the true M-dwarf TESS-band flux ratio), maximum Teff of a
   self-luminous occulter.
3. Phase curve / ellipsoidal / beaming: raw (not detrended) INV flux, transits masked, per-TESS-orbit quadratic
   baseline + background and pointing regressors, plus sin/cos at P (beaming/reflection) and P/2 (ellipsoidal).
   Empirical errors from the same regression at 300 random periods in 2-12 d. Predicted amplitudes for companion
   classes (Loeb & Gaudi 2003; Morris 1985 ellipsoidal).
Output: out/secondary_phase.json, out/secondary_scan.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
import batman
rng = np.random.default_rng(7)
z = np.load(f'{SD}/data/inv_full.npz')
fit = json.load(open(f'{SD}/out/fit_inv_free.json')); bst = fit['best']
T0, P, k, b, rho = bst['T0'], bst['P'], bst['k'], bst['b'], 10**bst['log10rho']
aR = aR_from_rho(rho, P)
t, f, sig, cad, sec_ = z['t'], z['f'], z['sig'], z['cad'], z['sector']
ph = ((t - T0)/P) % 1
intr = (np.abs(((t - T0)/P + 0.5) % 1 - 0.5)*P < 0.08)
# occultation template at a given phase offset
def template(phc):
    out = np.zeros(len(t))
    for c in np.unique(cad):
        m = cad == c
        p = batman.TransitParams(); p.t0 = T0 + phc*P; p.per = P; p.rp = k; p.a = aR; p.inc = np.degrees(np.arccos(b/aR)); p.ecc = 0; p.w = 90
        p.limb_dark = 'uniform'; p.u = []
        mm = batman.TransitModel(p, t[m], supersample_factor=SS[int(c)], exp_time=c/DAY)
        out[m] = 1 - mm.light_curve(p)
    return out/(k**2)
def lin_depth(phc):
    s = template(phc); use = (~intr) & np.isfinite(f)
    if s[use].max() <= 0: return np.nan, np.nan
    w = 1/sig[use]**2; r = f[use] - 1
    # per-window offset: fit depth with a local constant (only points within +-0.3 d of the event)
    near = np.abs(((t[use] - T0 - phc*P)/P + 0.5) % 1 - 0.5)*P < 0.3
    if near.sum() < 20 or s[use][near].max() <= 0: return np.nan, np.nan
    ev = np.round((t[use][near] - T0 - phc*P)/P).astype(int)
    _, wi = np.unique(ev, return_inverse=True)
    X = np.zeros((near.sum(), wi.max()+2)); X[np.arange(near.sum()), wi] = 1; X[:, -1] = -s[use][near]
    W = w[near]; A = X.T @ (X*W[:, None]); bb = X.T @ (W*r[near])
    sol = np.linalg.solve(A, bb); cov = np.linalg.inv(A)
    return sol[-1], np.sqrt(cov[-1, -1])
d05, e05 = lin_depth(0.5)
rand = []
while len(rand) < 400:
    phc = rng.uniform(0.06, 0.94)
    if abs(phc - 0.5) < 0.06: continue
    d, e = lin_depth(phc)
    if np.isfinite(d): rand.append((phc, d, e))
rand = np.array(rand)
emp = 1.4826*np.median(np.abs(rand[:, 1] - np.median(rand[:, 1])))
emp_std = np.std(rand[:, 1])
scan = []
for phc in np.arange(0.06, 0.9401, 0.004):
    d, e = lin_depth(phc); scan.append((phc, d, e))
scan = np.array(scan)
pd.DataFrame(scan, columns=['phase', 'depth', 'err_formal']).to_csv(f'{SD}/out/secondary_scan.csv', index=False)
ref = json.load(open('tic214311908/analysis/referee_phot/out/phase_scan_clean.json'))
# limits
lim3_inv = d05 + 3*max(emp, emp_std)
ref05 = (0.29e-3, 0.77e-3)
lim3_ref = ref05[0] + 3*ref05[1]
lim3 = min(lim3_inv, lim3_ref)
RpA = k/aR
tr = np.loadtxt(f'{SD}/data/filt_TESS.Red.dat'); wl = tr[:, 0]*1e-8; T_ = tr[:, 1]
def SB(T): return np.trapezoid(T_*wl/(wl**5*(np.exp(1.4388/(wl*T)) - 1)), wl)
Tg = np.linspace(300, 6000, 5000); sr = np.array([SB(x) for x in Tg])/SB(3350)
dep_tr = 1 - LC(f'{SD}/data/inv_win.npz').transit(T0, P, k, aR, b, bst['u1'], bst['u2']).min()
ratio_max = lim3/(k**2)  # planet/star surface-brightness ratio limit (uniform disks)
Tb_max = float(np.interp(ratio_max, sr, Tg))
# equilibrium temperature and expected secondary for reference
Teq0 = 3350*np.sqrt(1/(2*aR))
exp_reflect_Ag1 = RpA**2
exp_thermal_Teq = k**2*np.interp(Teq0*(4/3)**0.25, Tg, sr)  # dayside, no redistribution
# ---------------- phase curve
raw = z['raw']; rawf = 1 + (raw - 1)*935/979
ok = np.isfinite(rawf) & (~intr) & (np.abs(((t - T0)/P) % 1 - 0.5) > 0.04)
seg = np.zeros(len(t), int); br = np.where(np.diff(t) > 1.0)[0] + 1; 
for i, s0 in enumerate(np.split(np.arange(len(t)), br)): seg[s0] = i
def design(per, idx):
    tt = t[idx]; cols = []
    for sg in np.unique(seg[idx]):
        m = seg[idx] == sg; tm = tt[m].mean()
        for p_ in range(3):
            c = np.zeros(len(tt)); c[m] = ((tt[m] - tm))**p_; cols.append(c)
    ph_ = 2*np.pi*(tt - T0)/per
    cols += [np.sin(ph_), np.cos(ph_), np.sin(2*ph_), np.cos(2*ph_)]
    return np.array(cols).T
def pc_fit(per):
    idx = np.where(ok)[0]
    X = design(per, idx); w = 1/sig[idx]**2; y = rawf[idx] - 1
    A = X.T @ (X*w[:, None]); sol = np.linalg.lstsq(A, X.T @ (w*y), rcond=None)[0]
    a1 = np.hypot(sol[-4], sol[-3]); a2 = np.hypot(sol[-2], sol[-1])
    return sol[-4:], a1, a2
solP, A1, A2 = pc_fit(P)
rp = []
while len(rp) < 300:
    per = rng.uniform(2, 12)
    if min(abs(per - P), abs(per - P/2), abs(per - 2*P), abs(per-P/3)) < 0.15: continue
    s_, a1, a2 = pc_fit(per); rp.append((per, a1, a2, *s_))
rp = np.array(rp)
p95_1 = float(np.percentile(rp[:, 1], 95)); p95_2 = float(np.percentile(rp[:, 2], 95))
p997_1 = float(np.percentile(rp[:, 1], 99.7)); p997_2 = float(np.percentile(rp[:, 2], 99.7))
frac1 = float(np.mean(rp[:, 1] >= A1)); frac2 = float(np.mean(rp[:, 2] >= A2))
# predicted amplitudes (semi-amplitudes) for companions of the M3V (M1 = 0.37, R1 = 0.381, aR ~ 25.3, sin i ~ 1)
M1 = 0.370; R1 = 0.381
def K1(m2_msun):  # km/s
    return (2*np.pi*G/(P*DAY))**(1/3.)*m2_msun*MSUN/((M1 + m2_msun)*MSUN)**(2/3.)/1e5
xb = 1.4388/(8.0e-5*3350)
alpha_beam = xb*np.exp(xb)/(np.exp(xb) - 1)/4.
def A_beam(m2): return 4*alpha_beam*K1(m2)/2.998e5
def A_ell(m2):
    a = (G*(M1+m2)*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN
    u, tau = 0.5, 0.6
    alpha = 0.15*(15+u)*(1+tau)/(3-u)
    return alpha*(m2/M1)*(R1/a)**3
pred = {lab: dict(m2_msun=m, K1_kms=float(K1(m)), A_beam_ppt=float(A_beam(m)*1e3), A_ell_ppt=float(A_ell(m)*1e3))
        for lab, m in [('Neptune_17ME', 17*3.0e-6), ('Jupiter', 9.55e-4), ('BD_30MJ', 0.0286), ('BD_70MJ', 0.0668), ('M_0.1', 0.1), ('M3V_twin_0.37', 0.37), ('WD_0.6', 0.6), ('WD_1.0', 1.0)]}
out = dict(secondary=dict(inv_depth_ppt=float(d05*1e3), inv_err_formal_ppt=float(e05*1e3), inv_err_empirical_ppt=float(max(emp, emp_std)*1e3),
                          inv_err_emp_robust_ppt=float(emp*1e3), inv_err_emp_std_ppt=float(emp_std*1e3), n_random=len(rand),
                          ref_depth_ppt=0.29, ref_err_ppt=0.77, ref_source='../referee_phot/out/phase_scan_clean.json (+0.28 +- 0.75 e-/s)',
                          upper3_inv_ppt=float(lim3_inv*1e3), upper3_ref_ppt=float(lim3_ref*1e3), adopted_upper3_ppt=float(lim3*1e3),
                          scan_max_abs_snr_empirical=float(np.nanmax(np.abs(scan[:, 1]))/max(emp, emp_std)), scan_max_depth_ppt=float(np.nanmax(scan[:, 1])*1e3),
                          scan_phase_of_max=float(scan[np.nanargmax(scan[:, 1]), 0]),
                          scan_max_pos_snr_empirical=float(np.nanmax(scan[:, 1])/max(emp, emp_std)), scan_min_snr_empirical=float(np.nanmin(scan[:, 1])/max(emp, emp_std)),
                          scan_phase_of_min=float(scan[np.nanargmin(scan[:, 1]), 0]),
                          p_random_ge_scanmax=float(np.mean(rand[:, 1] >= np.nanmax(scan[:, 1]))),
                          n_independent_trials=float((0.94-0.06)/((fit['post']['T14_h']['med']/24)/P)),
                          p_global_scanmax=float(1 - (1 - max(np.mean(rand[:, 1] >= np.nanmax(scan[:, 1])), 1/len(rand)))**((0.94-0.06)/((fit['post']['T14_h']['med']/24)/P))),
                          any_phase_upper3_ppt=float((np.nanmax(scan[:, 1]) + 3*max(emp, emp_std))*1e3)),
           limits=dict(Rp_over_a=float(RpA), Ag_upper3=float(lim3/RpA**2), Tb_upper3_K_blackbody=Tb_max, surface_brightness_ratio_upper3=float(ratio_max),
                       transit_depth_model=float(dep_tr), Teq_K_A0_full=float(Teq0), expected_reflection_Ag1_ppm=float(exp_reflect_Ag1*1e6),
                       expected_thermal_dayside_ppm=float(exp_thermal_Teq*1e6)),
           phase_curve=dict(A_P_ppt=float(A1*1e3), A_Phalf_ppt=float(A2*1e3), coeffs=list(map(float, solP)), A_P_95_ppt=p95_1*1e3, A_Phalf_95_ppt=p95_2*1e3,
                            A_P_997_ppt=p997_1*1e3, A_Phalf_997_ppt=p997_2*1e3, frac_random_ge_obs_P=frac1, frac_random_ge_obs_Phalf=frac2,
                            note='raw INV flux, per-orbit quadratic baselines; random-period calibration'),
           predicted_companion_signals=pred, alpha_beam=float(alpha_beam))
json.dump(out, open(f'{SD}/out/secondary_phase.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
