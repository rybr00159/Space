"""v3: combine the six sectors into the sky offset of the eclipsing source, with a lateral-colour nuisance model.

Data: per sector s and eclipse type t the measured offset o_st (dE, dN) of the difference-image source from the
direct-image (G-star) position, bias-corrected with the injection-recovery bias at the target ('emp' injections),
noise = the injection trials' recovered-minus-injected offsets (zero-mean), covariance C_st from those trials.

Chromatic term: a source with SED redder than the G star is displaced radially AWAY from the camera boresight by
   chi_st = c_t * (k_s / k_ref) * (-u_s)        [u_s = unit sky vector towards the boresight in sector s]
with k_s the field-star lateral-colour coefficient of that sector (f02) and c_t >= 0 a free 'colour lever' (arcsec at
k_ref) per eclipse type.  (k_tan ~ 0 and is neglected; variant adds it.)
Models: M0 o = delta;  M1 o = chi (delta = 0);  M2 o = delta + chi;  MT tangential components only (o.t = delta.t).
Calibration: simulated data = model truth + one random injection-trial noise draw per sector/type (non-Gaussian tails
kept); refits give the empirical delta covariance and the distribution of the profile statistic
dchi2(delta_true) = chi2(delta fixed) - chi2_min, used for confidence regions and separation exclusions.
"""
import sys, os, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v3lib import V3, SECTORS, to_json

TYPES = os.environ.get('V3_TYPES', 'deep,shallow').split(',')
DROP = [int(s) for s in os.environ.get('V3_DROP', '').split(',') if s]
KMODE = os.environ.get('V3_KMODE', 'f02')        # f02 | equal | theta
FCOL = os.environ.get('V3_FCOL', 'f02_colour_theta.json')
NOBIAS = os.environ.get('V3_NOBIAS', '0') == '1'
ERRSCALE = float(os.environ.get('V3_ERRSCALE', 1.0))
LABEL = os.environ.get('V3_LABEL', 'baseline')
NS = int(os.environ.get('V3_NS', 3000))
secs = [s for s in SECTORS if s not in DROP]

col = json.load(open(f'{V3}/results/{FCOL}'))
kf = {}
for fr in col['frames']:
    kf[fr['sec']] = dict(k=fr['k_rad'], pa=fr['pa_t'], th=fr['theta_t'])
obs, cov, U, KR, lab, noise = [], [], [], [], [], []
for s in secs:
    real = json.load(open(f"{V3}/results/real_s{s}{os.environ.get('V3_REALTAG', '')}.json"))
    inj = json.load(open(f'{V3}/results/inj_s{s}.json'))
    z = np.load(f'{V3}/results/inj_s{s}.npz')
    pa = np.radians(kf[s]['pa'])
    u = np.array([np.sin(pa), np.cos(pa)])          # towards boresight (E, N)
    if KMODE == 'f02':
        k = kf[s]['k']
    elif KMODE == 'equal':
        k = 1.0
    else:
        k = kf[s]['th']
    for t in TYPES:
        r = real[t]
        nz = z[f'{t}_emp'][0]                         # config 0 = injected at the target
        nz = nz[np.all(np.isfinite(nz), 1)]
        b = np.zeros(2) if NOBIAS else nz.mean(0)
        nz = (nz - nz.mean(0)) * ERRSCALE
        obs.append(np.array([r['dE'], r['dN']]) - b)
        cov.append(np.cov(nz.T)); noise.append(nz); U.append(u); KR.append(k); lab.append(f'S{s}_{t}')
obs = np.array(obs); cov = np.array(cov); U = np.array(U); KR = np.array(KR)
KR = KR / np.mean(KR)
N = len(obs)
tidx = np.array([TYPES.index(l.split('_')[1]) for l in lab])
Wh = np.array([np.linalg.cholesky(np.linalg.inv(c)).T for c in cov])    # r^T C^-1 r = |Wh r|^2
nt = len(TYPES)


def design(kind):
    """Linear model o = A p (stacked over sectors, whitened later). Params: [dE, dN] + [c_t]."""
    cols = []
    if kind in ('M0', 'M2'):
        for j in range(2):
            e = np.zeros((N, 2)); e[:, j] = 1; cols.append(e)
    if kind in ('M1', 'M2'):
        for t in range(nt):
            e = np.zeros((N, 2)); m = tidx == t
            e[m] = -(KR[m][:, None] * U[m])         # away from boresight for c > 0
            cols.append(e)
    return np.stack(cols, -1)                        # (N, 2, P)


def solve(kind, o, fix_delta=None):
    A = design(kind)
    y = o.copy()
    if fix_delta is not None:
        y = y - fix_delta[None, :]
        A = A[:, :, 2:] if kind == 'M2' else A
    Aw = np.einsum('nij,njp->nip', Wh, A).reshape(2 * N, -1)
    yw = np.einsum('nij,nj->ni', Wh, y).reshape(-1)
    if Aw.shape[1] == 0:
        return np.zeros(0), float(yw @ yw)
    p, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    r = yw - Aw @ p
    return p, float(r @ r)


def tangential(o):
    """MT: only the components perpendicular to each sector's radial direction (immune to the radial chromatic term)."""
    T = np.column_stack([U[:, 1], -U[:, 0]])
    st2 = np.array([t @ c @ t for t, c in zip(T, cov)])
    A = T / np.sqrt(st2)[:, None]; y = (o * T).sum(1) / np.sqrt(st2)
    p, *_ = np.linalg.lstsq(A, y, rcond=None)
    return p, float(((y - A @ p) ** 2).sum())


res = {}
for kind in ('M0', 'M1', 'M2'):
    p, c2 = solve(kind, obs)
    res[kind] = dict(p=p, chi2=c2, dof=2 * N - len(p))
pT, c2T = tangential(obs)
res['MT'] = dict(p=pT, chi2=c2T, dof=N - 2)
rng = np.random.default_rng(11)


def sim(truth):
    return truth + np.array([nz[rng.integers(len(nz))] for nz in noise])


A2 = design('M2')
pM2 = res['M2']['p']; truthM2 = np.einsum('nip,p->ni', A2, pM2)
sims = []; chis = {'M2': [], 'M1': [], 'M0': []}; simT = []
for _ in range(NS):
    o = sim(truthM2)
    p, c2 = solve('M2', o); sims.append(p); chis['M2'].append(c2)
    simT.append(tangential(o)[0])
sims = np.array(sims); simT = np.array(simT)
Cd = np.cov(sims[:, :2].T)
# null: delta = 0 (M1 truth) -> distribution of the M2 delta estimate and of dchi2(0)
A1 = design('M1'); pM1 = res['M1']['p']; truth1 = np.einsum('nip,p->ni', A1, pM1)
d0 = []; dchi0 = []; c1n = []
for _ in range(NS):
    o = sim(truth1)
    p, c2 = solve('M2', o); d0.append(p[:2])
    _, c2f = solve('M2', o, fix_delta=np.zeros(2)); dchi0.append(c2f - c2)
    c1n.append(solve('M1', o)[1])
d0 = np.array(d0); dchi0 = np.array(dchi0); c1n = np.array(c1n)
_, c2fix0 = solve('M2', obs, fix_delta=np.zeros(2))
obs_dchi0 = c2fix0 - res['M2']['chi2']
# profile map over delta
g = np.arange(-3.0, 3.0001, 0.05)
M = np.zeros((len(g), len(g)))
for i, dn in enumerate(g):
    for j, de in enumerate(g):
        M[i, j] = solve('M2', obs, fix_delta=np.array([de, dn]))[1] - res['M2']['chi2']
# calibrated thresholds: distribution of dchi2(delta_true) under M2 truth
dtrue = []
for _ in range(min(NS, 2000)):
    o = sim(truthM2)
    p, c2 = solve('M2', o)
    _, c2f = solve('M2', o, fix_delta=pM2[:2])
    dtrue.append(c2f - c2)
dtrue = np.array(dtrue)
thr = {q: float(np.percentile(dtrue, q)) for q in (68.27, 95.45, 99.73)}
EE, NNg = np.meshgrid(g, g); RR = np.hypot(EE, NNg); PA = np.degrees(np.arctan2(EE, NNg)) % 360
prof = []
for r0 in np.arange(0, 3.01, 0.1):
    m = np.abs(RR - r0) <= 0.026 if r0 > 0 else RR < 0.026
    prof.append(dict(sep=float(r0), min_dchi2=float(M[m].min())))
region = {}
for q, tv in thr.items():
    inside = M <= tv
    region[q] = dict(thr=tv, max_sep=float(RR[inside].max()) if inside.any() else None, min_sep=float(RR[inside].min()) if inside.any() else None)
    if inside.any() and RR[inside].min() > 0.06:
        a = np.radians(PA[inside]); mpa = np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())) % 360
        d = (PA[inside] - mpa + 180) % 360 - 180
        region[q]['pa_range'] = [float((mpa + d.min()) % 360), float((mpa + d.max()) % 360)]
# separation exclusion p-values (empirical) for truths at each separation along the best PA (most conservative)
excl = []
for r0 in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
    m = np.abs(RR - r0) <= 0.026
    k = np.argmin(np.where(m, M, np.inf)); dt = np.array([EE.flat[k], NNg.flat[k]])
    obs_d = M.flat[k]
    tr = dt[None, :] + np.einsum('nip,p->ni', A2[:, :, 2:], solve('M2', obs, fix_delta=dt)[0])
    st = []
    for _ in range(1000):
        o = sim(tr)
        _, c2 = solve('M2', o); _, c2f = solve('M2', o, fix_delta=dt); st.append(c2f - c2)
    st = np.array(st)
    excl.append(dict(sep=r0, best_pa=float(np.degrees(np.arctan2(*dt)) % 360), obs_dchi2=float(obs_d), p_emp=float((st >= obs_d).mean()),
                     n=len(st), sim_q99=float(np.percentile(st, 99))))
out = dict(label=LABEL, types=TYPES, sectors=secs, kmode=KMODE, labels=lab, obs=obs, cov=cov, KR_rel=KR,
           u_to_boresight=U, models=res,
           M2_delta=pM2[:2], M2_c=pM2[2:], M2_delta_emp_cov=Cd, M2_delta_emp_sd=np.sqrt(np.diag(Cd)), M2_c_emp_sd=sims[:, 2:].std(0),
           M2_sep=float(np.hypot(*pM2[:2])), M2_pa=float(np.degrees(np.arctan2(*pM2[:2])) % 360),
           M2_chi2_sim_pct=[float(np.percentile(chis['M2'], q)) for q in (5, 50, 95, 99)],
           M2_chi2_p_emp=float((np.array(chis['M2']) >= res['M2']['chi2']).mean()),
           M1_chi2_null_p_emp=float((c1n >= res['M1']['chi2']).mean()),
           delta0_obs_dchi2=obs_dchi0, delta0_p_emp=float((dchi0 >= obs_dchi0).mean()),
           delta_if_on_target_sep_pct=[float(np.percentile(np.hypot(*d0.T), q)) for q in (68.27, 95.45, 99.73)],
           MT_delta=pT, MT_delta_emp_sd=simT.std(0), MT_chi2=c2T,
           thresholds=thr, region=region, profile=prof, exclusion=excl)
to_json(out, f'{V3}/results/s05_{LABEL}.json')
np.savez_compressed(f'{V3}/results/s05_{LABEL}_map.npz', g=g, M=M, sims=sims)
print(f"[{LABEL}] sectors {secs} types {TYPES} kmode {KMODE}")
for kind in ('M0', 'M1', 'M2', 'MT'):
    print(f"  {kind}: p={np.round(res[kind]['p'], 3)} chi2={res[kind]['chi2']:.1f}/{res[kind]['dof']}")
print(f"  M2 delta = ({pM2[0]:+.2f}, {pM2[1]:+.2f})\" emp sd ({np.sqrt(Cd[0, 0]):.2f}, {np.sqrt(Cd[1, 1]):.2f}) corr {Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1]):+.2f}; "
      f"sep {out['M2_sep']:.2f}\" PA {out['M2_pa']:.0f}; c = {np.round(pM2[2:], 2)} +- {np.round(sims[:, 2:].std(0), 2)}")
print(f"  M2 chi2 p_emp = {out['M2_chi2_p_emp']:.3f};  M1 (delta=0) chi2 p_emp = {out['M1_chi2_null_p_emp']:.3f};  "
      f"dchi2(delta=0) = {obs_dchi0:.2f}, p_emp = {out['delta0_p_emp']:.3f}")
print(f"  thresholds (68/95/99.7%): {np.round(list(thr.values()), 2)};  region: " +
      '; '.join(f"{q}%: sep {v['min_sep']:.2f}-{v['max_sep']:.2f}" + (f" PA {v['pa_range'][0]:.0f}-{v['pa_range'][1]:.0f}" if 'pa_range' in v else '') for q, v in region.items()))
print('  exclusion: ' + ', '.join(f"{e['sep']}\":dchi2={e['obs_dchi2']:.1f},p={e['p_emp']:.4f}" for e in excl))
print(f"  MT (tangential only) delta = {np.round(pT, 2)} +- {np.round(simT.std(0), 2)}  chi2 {c2T:.1f}/{N - 2}")
