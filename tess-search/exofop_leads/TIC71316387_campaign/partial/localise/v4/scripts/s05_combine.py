"""v4: combine the six sectors into the sky offset delta of the eclipsing source from the G star, with a lateral-colour
(chromatic) nuisance term, calibrated entirely with the injection-recovery noise.

Per sector s and eclipse type t:  o_st = measured (dE, dN) offset (s02)  minus  b_st = mean recovered-minus-injected
offset of target-position injections (s03, 'emp').  Noise model = the empirical distribution of those injection residuals.
Model:   o_st = delta + c_t * kappa_s * (-u_s) + noise
   u_s     unit sky vector from the target towards the camera boresight in sector s (tess-point; f01 tables);
   kappa_s relative lateral-colour gain at the target in sector s (FFI field-star k_rad, normalised to mean 1; f02),
           variants: equal gains, gains proportional to the field angle theta;
   c_t     colour lever (arcsec) of the eclipsed light of type t (free; c > 0 = redder than the G star).
Models: M0 delta only (fixed sky offset, no chromatic term);  M1 chromatic only (delta = 0);  M2 delta + chromatic.
Calibration: parametric bootstrap with the injection residuals (one random trial per sector/type): empirical covariance
of delta, distribution of the profile statistic dchi2(delta_true) -> calibrated confidence regions, separation
exclusions (worst-case PA), p-values for M0/M1 goodness of fit, and dchi2(delta=0).
Env: V4C_TYPES (deep,shallow), V4C_DROP (sectors), V4C_KAPPA (ffi|equal|theta), V4C_NOBIAS, V4C_ERRX, V4C_LABEL, V4C_NS.
"""
import sys, os, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v4lib import V4, SECTORS, to_json

TYPES = os.environ.get('V4C_TYPES', 'deep,shallow').split(',')
DROP = [int(s) for s in os.environ.get('V4C_DROP', '').split(',') if s]
KAPPA = os.environ.get('V4C_KAPPA', 'ffi')
NOBIAS = os.environ.get('V4C_NOBIAS', '0') == '1'
ERRX = float(os.environ.get('V4C_ERRX', 1.0))
LABEL = os.environ.get('V4C_LABEL', 'baseline')
NS = int(os.environ.get('V4C_NS', 4000))
TAG = os.environ.get('V4_TAG', '')
REALTAG = os.environ.get('V4C_REALTAG', TAG)
FCOL = os.environ.get('V4C_FCOL', 'f02_colour.json')
secs = [s for s in SECTORS if s not in DROP]

col = json.load(open(f'{V4}/results/{FCOL}'))
ksec = {int(k): v for k, v in col['sectors'].items()}
obs, cov, U, KAP, lab, noise, tid = [], [], [], [], [], [], []
for s in secs:
    real = json.load(open(f'{V4}/results/real_s{s}{REALTAG}.json'))
    z = np.load(f'{V4}/results/inj_s{s}{TAG}.npz')
    pa = np.radians(ksec[s]['pa_t'])
    u = np.array([np.sin(pa), np.cos(pa)])
    if KAPPA == 'red':
        rr = json.load(open(f'{V4}/results/f03_persec_red.json'))
        k = ksec[s]['k_rad'] * rr[f'{s}_2.0-3.5']['median_ratio']
    else:
        k = {'ffi': ksec[s]['k_rad'], 'equal': 1.0, 'theta': ksec[s]['theta_t']}[KAPPA]
    for t in TYPES:
        nz = z[f'{t}_emp'][0]
        nz = nz[np.all(np.isfinite(nz), 1)]
        b = np.zeros(2) if NOBIAS else nz.mean(0)
        nz = (nz - nz.mean(0)) * ERRX
        obs.append(np.array([real[t]['dE'], real[t]['dN']]) - b)
        cov.append(np.cov(nz.T)); noise.append(nz); U.append(u); KAP.append(k); lab.append(f'S{s}_{t}'); tid.append(TYPES.index(t))
obs = np.array(obs); cov = np.array(cov); U = np.array(U); KAP = np.array(KAP); tid = np.array(tid)
KAP = KAP / KAP.mean()
N = len(obs); nt = len(TYPES)
Wh = np.array([np.linalg.cholesky(np.linalg.inv(c)).T for c in cov])


def design(kind):
    cols = []
    if kind in ('M0', 'M2'):
        for j in range(2):
            e = np.zeros((N, 2)); e[:, j] = 1; cols.append(e)
    if kind in ('M1', 'M2'):
        for t in range(nt):
            e = np.zeros((N, 2)); m = tid == t
            e[m] = -(KAP[m][:, None] * U[m])
            cols.append(e)
    return np.stack(cols, -1)


DES = {k: design(k) for k in ('M0', 'M1', 'M2')}


def solve(kind, o, fix=None):
    A = DES[kind]; y = o
    if fix is not None:
        y = o - fix[None, :]
        A = A[:, :, 2:]
    Aw = np.einsum('nij,njp->nip', Wh, A).reshape(2 * N, -1)
    yw = np.einsum('nij,nj->ni', Wh, y).reshape(-1)
    if Aw.shape[1] == 0:
        return np.zeros(0), float(yw @ yw)
    p, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    r = yw - Aw @ p
    return p, float(r @ r)


rng = np.random.default_rng(int(os.environ.get('V4C_SEED', 3)))


def draw():
    return np.array([nz[rng.integers(len(nz))] for nz in noise])


res = {k: dict(zip(('p', 'chi2'), solve(k, obs))) for k in ('M0', 'M1', 'M2')}
for k in res:
    res[k]['dof'] = 2 * N - DES[k].shape[2]
pM2 = res['M2']['p']; truth2 = np.einsum('nip,p->ni', DES['M2'], pM2)
pM1 = res['M1']['p']; truth1 = np.einsum('nip,p->ni', DES['M1'], pM1)
pM0 = res['M0']['p']; truth0 = np.einsum('nip,p->ni', DES['M0'], pM0)
sims, chi2sim, dtrue = [], [], []
for _ in range(NS):
    o = truth2 + draw()
    p, c2 = solve('M2', o); sims.append(p); chi2sim.append(c2)
    dtrue.append(solve('M2', o, fix=pM2[:2])[1] - c2)
sims = np.array(sims); dtrue = np.array(dtrue); chi2sim = np.array(chi2sim)
Cd = np.cov(sims[:, :2].T)
thr = {q: float(np.percentile(dtrue, q)) for q in (68.27, 95.45, 99.73)}
# null delta = 0 (M1 truth)
d0, dchi0, c1n = [], [], []
for _ in range(NS):
    o = truth1 + draw()
    p, c2 = solve('M2', o); d0.append(p[:2])
    dchi0.append(solve('M2', o, fix=np.zeros(2))[1] - c2)
    c1n.append(solve('M1', o)[1])
d0 = np.array(d0); dchi0 = np.array(dchi0); c1n = np.array(c1n)
obs_dchi0 = solve('M2', obs, fix=np.zeros(2))[1] - res['M2']['chi2']
# M0 truth (fixed offset, no chromatic term): chi2 distribution of M0 fits
c0n = np.array([solve('M0', truth0 + draw())[1] for _ in range(NS)])
# profile map
g = np.arange(-3.0, 3.0001, 0.05)
M = np.zeros((len(g), len(g)))
for i, dn in enumerate(g):
    for j, de in enumerate(g):
        M[i, j] = solve('M2', obs, fix=np.array([de, dn]))[1] - res['M2']['chi2']
EE, NNg = np.meshgrid(g, g); RR = np.hypot(EE, NNg); PA = np.degrees(np.arctan2(EE, NNg)) % 360
region = {}
for q, tv in thr.items():
    ins = M <= tv
    region[q] = dict(thr=tv, max_sep=float(RR[ins].max()), min_sep=float(RR[ins].min()), area_arcsec2=float(ins.sum() * 0.05 ** 2))
    if RR[ins].min() > 0.07:
        a = np.radians(PA[ins]); mpa = np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())) % 360
        dd = (PA[ins] - mpa + 180) % 360 - 180
        region[q]['pa_range'] = [float((mpa + dd.min()) % 360), float((mpa + dd.max()) % 360)]
# separation exclusion: at each separation take the best (least excluded) PA; calibrate with simulations at that truth
excl = []
for r0 in (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0):
    m = np.abs(RR - r0) <= 0.026
    k = np.argmin(np.where(m, M, np.inf)); dt = np.array([EE.flat[k], NNg.flat[k]])
    od = float(M.flat[k])
    ptr, _ = solve('M2', obs, fix=dt)
    tr = dt[None, :] + np.einsum('nip,p->ni', DES['M2'][:, :, 2:], ptr)
    st = []
    for _ in range(1500):
        o = tr + draw()
        st.append(solve('M2', o, fix=dt)[1] - solve('M2', o)[1])
    st = np.array(st)
    excl.append(dict(sep=r0, best_pa=float(np.degrees(np.arctan2(*dt)) % 360), obs_dchi2=od, p_emp=float((st >= od).mean()), n=len(st)))
# PA marginal for a given separation ring (where is the source if it is at 0.5-1.2"?)
pa_ring = {}
for r_lo, r_hi in ((0.3, 0.7), (0.7, 1.2), (0.5, 1.2)):
    m = (RR >= r_lo) & (RR <= r_hi)
    L = np.exp(-0.5 * M) * m
    pab = np.arange(0, 360, 30)
    hist = [float(L[(PA >= a) & (PA < a + 30)].sum() / L.sum()) for a in pab]
    pa_ring[f'{r_lo}-{r_hi}'] = dict(pa_bins_start=pab.tolist(), prob=hist)
# radial / tangential decomposition per sector (model-free view)
dec = []
for i in range(N):
    ur = -U[i]; ut = np.array([U[i][1], -U[i][0]])
    sr = np.sqrt(ur @ cov[i] @ ur); stn = np.sqrt(ut @ cov[i] @ ut)
    dec.append(dict(label=lab[i], away=float(obs[i] @ ur), away_sd=float(sr), tan=float(obs[i] @ ut), tan_sd=float(stn), kappa=float(KAP[i])))
# tangential-only estimate (immune to any radial term, whatever its amplitude per sector)
T = np.column_stack([U[:, 1], -U[:, 0]])
st2 = np.array([t @ c @ t for t, c in zip(T, cov)])
At = T / np.sqrt(st2)[:, None]; yt = (obs * T).sum(1) / np.sqrt(st2)
pT, *_ = np.linalg.lstsq(At, yt, rcond=None)
simT = []
for _ in range(NS):
    o = truth2 + draw(); y = (o * T).sum(1) / np.sqrt(st2)
    simT.append(np.linalg.lstsq(At, y, rcond=None)[0])
simT = np.array(simT)
# leave-one-sector-out predictions under M2
loo = []
for s in secs:
    keep = np.array([not l.startswith(f'S{s}_') for l in lab])
    if keep.sum() < 4:
        continue
    Aw = np.einsum('nij,njp->nip', Wh[keep], DES['M2'][keep]).reshape(-1, DES['M2'].shape[2])
    yw = np.einsum('nij,nj->ni', Wh[keep], obs[keep]).reshape(-1)
    p, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    pred = np.einsum('nip,p->ni', DES['M2'][~keep], p)
    z2 = float(sum((o - pr) @ np.linalg.inv(c) @ (o - pr) for o, pr, c in zip(obs[~keep], pred, cov[~keep])))
    loo.append(dict(sector=s, delta=p[:2].tolist(), c=p[2:].tolist(), pred=pred.tolist(), obs=obs[~keep].tolist(), chi2_pred=z2, n=int((~keep).sum() * 2)))
out = dict(label=LABEL, types=TYPES, sectors=secs, kappa_mode=KAPPA, labels=lab, obs_bias_corrected=obs, cov=cov, kappa=KAP, u_to_boresight=U,
           models=res, M0_chi2_p_emp=float((c0n >= res['M0']['chi2']).mean()), M0_chi2_sim_q=[float(np.percentile(c0n, q)) for q in (50, 95, 99.9)],
           M1_chi2_p_emp=float((c1n >= res['M1']['chi2']).mean()), M2_chi2_p_emp=float((chi2sim >= res['M2']['chi2']).mean()),
           delta=pM2[:2], delta_cov_emp=Cd, delta_sd_emp=np.sqrt(np.diag(Cd)), delta_corr=float(Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1])),
           c=pM2[2:], c_sd_emp=sims[:, 2:].std(0), sep=float(np.hypot(*pM2[:2])), pa=float(np.degrees(np.arctan2(*pM2[:2])) % 360),
           delta0_dchi2=obs_dchi0, delta0_p_emp=float((dchi0 >= obs_dchi0).mean()),
           sep_if_on_target_q=[float(np.percentile(np.hypot(*d0.T), q)) for q in (50, 68.27, 95.45, 99.73)],
           thresholds=thr, region=region, exclusion=excl, pa_ring=pa_ring, decomposition=dec,
           tangential_only=dict(delta=pT, sd=simT.std(0)), loo=loo)
if nt == 2:
    dc = sims[:, 3] - sims[:, 2]
    out['c_shallow_minus_deep'] = dict(value=float(pM2[3] - pM2[2]), sd_emp=float(dc.std()))
    # null: equal levers (c_s = c_d): refit with a common lever and simulate
    Acom = DES['M2'].copy(); Acom = np.concatenate([Acom[:, :, :2], Acom[:, :, 2:3] + Acom[:, :, 3:4]], -1)
    def solve_com(o):
        Aw = np.einsum('nij,njp->nip', Wh, Acom).reshape(2 * N, -1); yw = np.einsum('nij,nj->ni', Wh, o).reshape(-1)
        p, *_ = np.linalg.lstsq(Aw, yw, rcond=None); r = yw - Aw @ p; return p, float(r @ r)
    pc, c2c = solve_com(obs)
    trc = np.einsum('nip,p->ni', Acom, pc)
    dd = []
    for _ in range(NS):
        o = trc + draw(); dd.append(solve_com(o)[1] - solve('M2', o)[1])
    dd = np.array(dd); od = c2c - res['M2']['chi2']
    out['equal_lever_test'] = dict(dchi2=float(od), p_emp=float((dd >= od).mean()), common_c=float(pc[2]))
to_json(out, f'{V4}/results/s05_{LABEL}.json')
np.savez_compressed(f'{V4}/results/s05_{LABEL}_map.npz', g=g, M=M, sims=sims)
print(f"[{LABEL}] sectors {secs} types {TYPES} kappa {KAPPA}  N={N}")
for k in ('M0', 'M1', 'M2'):
    print(f"  {k}: p={np.round(res[k]['p'], 3)} chi2={res[k]['chi2']:.1f}/{res[k]['dof']}")
print(f"  p_emp: M0 {out['M0_chi2_p_emp']:.4f} (sim chi2 q50/95/99.9 {np.round(out['M0_chi2_sim_q'], 1)}), M1 {out['M1_chi2_p_emp']:.3f}, M2 {out['M2_chi2_p_emp']:.3f}")
print(f"  delta = ({pM2[0]:+.2f}, {pM2[1]:+.2f})\" sd ({np.sqrt(Cd[0, 0]):.2f}, {np.sqrt(Cd[1, 1]):.2f}) corr {out['delta_corr']:+.2f}; sep {out['sep']:.2f}\" PA {out['pa']:.0f}; "
      f"c = {np.round(pM2[2:], 2)} +- {np.round(sims[:, 2:].std(0), 2)}")
print(f"  dchi2(delta=0) = {obs_dchi0:.2f}, p_emp = {out['delta0_p_emp']:.3f};  |delta| if on target q50/68/95/99.7: {np.round(out['sep_if_on_target_q'], 2)}")
print(f"  thresholds {np.round(list(thr.values()), 2)}; " + '; '.join(f"{q}%: sep {v['min_sep']:.2f}-{v['max_sep']:.2f}\" area {v['area_arcsec2']:.2f}" +
      (f" PA {v['pa_range'][0]:.0f}-{v['pa_range'][1]:.0f}" if 'pa_range' in v else '') for q, v in region.items()))
print('  exclusion (worst PA): ' + ', '.join(f"{e['sep']}\":{e['obs_dchi2']:.1f}/p={e['p_emp']:.4f}" for e in excl))
print(f"  tangential-only delta = {np.round(pT, 2)} +- {np.round(simT.std(0), 2)}")
if nt == 2:
    print(f"  c_shallow - c_deep = {out['c_shallow_minus_deep']['value']:+.2f} +- {out['c_shallow_minus_deep']['sd_emp']:.2f};  equal-lever test dchi2 {out['equal_lever_test']['dchi2']:.1f} p_emp {out['equal_lever_test']['p_emp']:.4f}")
for d in dec:
    print(f"   {d['label']:12s} away-from-boresight {d['away']:+.2f} +- {d['away_sd']:.2f}   tangential {d['tan']:+.2f} +- {d['tan_sd']:.2f}   kappa {d['kappa']:.2f}")
for l in loo:
    print(f"   LOO S{l['sector']}: delta {np.round(l['delta'], 2)}  chi2_pred {l['chi2_pred']:.1f}/{l['n']}")
