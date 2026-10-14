"""Robustness variants of the sky-offset fit (model M2: delta + chromatic levers) and a profile chi2 map of delta.

Variants: lateral-colour scaling from narrow box (baseline, x11) / wide box (x12) / equal in all sectors / proportional
to field angle; k_tan = 0; leave-one-sector-out; deep only; shallow only; PRF-injection covariance; no bias correction;
covariance inflated x1.5.
Profile chi2 map over delta (c_deep, c_shallow profiled) for the baseline -> exclusion radii and PA constraint.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
from loclib import R, SECTORS, to_json

TYPES = ['deep', 'shallow']
real = json.load(open(f'{R}/results/real_fast.json'))
bs = json.load(open(f'{R}/results/x02_boresight.json'))
x11 = json.load(open(f'{R}/results/x11_colour_terms.json'))
x12 = json.load(open(f'{R}/results/x12_red_end.json'))


def kset(name):
    kr, kt = {}, {}
    for s in SECTORS:
        if name == 'x12':
            rr = [c for c in x12['frames'] if c['sec'] == s]
        else:
            rr = [c for c in x11 if c['sec'] == s]
        kr[s] = np.mean([c['k_rad'] for c in rr]); kt[s] = np.mean([c['k_tan'] for c in rr])
    if name == 'equal':
        m = np.mean(list(kr.values())); kr = {s: m for s in SECTORS}; kt = {s: 0.0 for s in SECTORS}
    if name == 'theta':
        th = {s: bs[str(s)]['field_angle_deg'] for s in SECTORS}; mt = np.mean(list(th.values())); m = np.mean(list(kr.values()))
        kr = {s: m * th[s] / mt for s in SECTORS}; kt = {s: 0.0 for s in SECTORS}
    return kr, kt


def build(ks='x11', ktan=True, sectors=SECTORS, types=TYPES, covmode='emp', bias=True, infl=1.0):
    kr, kt = kset(ks)
    kref = np.mean([kr[s] for s in SECTORS])
    obs, cov, U, Tn, KR, KT, ti = [], [], [], [], [], [], []
    for s in sectors:
        inj = json.load(open(f'{R}/results/inj_s{s}.json'))
        pa = np.radians(bs[str(s)]['pa_to_boresight'])
        u = np.array([np.sin(pa), np.cos(pa)]); tn = np.array([u[1], -u[0]])
        for t in types:
            r = real[str(s)][t]; s0 = inj[t]['summary'][0]
            b = np.array(s0[covmode]['bias']) if bias else np.zeros(2)
            obs.append(np.array([r['dE'], r['dN']]) - b); cov.append(np.array(s0[covmode]['cov']) * infl ** 2)
            U.append(u); Tn.append(tn); KR.append(kr[s] / kref); KT.append((kt[s] / kref) if ktan else 0.0); ti.append(TYPES.index(t))
    return dict(obs=np.array(obs), cov=np.array(cov), U=np.array(U), Tn=np.array(Tn), KR=np.array(KR), KT=np.array(KT),
                ti=np.array(ti), types=types)


def fit(D, fix_delta=None):
    W = np.array([np.linalg.cholesky(np.linalg.inv(c)) for c in D['cov']])
    nt = len(D['types'])
    tmap = np.array([D['types'].index(TYPES[i]) for i in D['ti']])

    def res(p):
        d = np.zeros_like(D['obs'])
        delta = p[:2] if fix_delta is None else fix_delta
        c = np.asarray(p[-nt:])[tmap]
        d += delta
        d += (c * D['KR'])[:, None] * D['U'] + (c * D['KT'])[:, None] * D['Tn']
        r = D['obs'] - d
        return np.concatenate([W[i].T @ r[i] for i in range(len(r))])
    p0 = ([0, 0] if fix_delta is None else []) + [-3.0] * nt
    f = least_squares(res, p0)
    C = np.linalg.inv(f.jac.T @ f.jac)
    return f.x, float((f.fun ** 2).sum()), C


if __name__ == '__main__':
    variants = {'baseline (x11 k)': {}, 'k from wide box (x12)': dict(ks='x12'), 'k equal all sectors': dict(ks='equal'),
                'k proportional to field angle': dict(ks='theta'), 'k_tan = 0': dict(ktan=False),
                'deep only': dict(types=['deep']), 'shallow only': dict(types=['shallow']),
                'PRF-injection covariance': dict(covmode='prf'), 'no bias correction': dict(bias=False),
                'errors x1.5': dict(infl=1.5)}
    for s in SECTORS:
        variants[f'drop S{s}'] = dict(sectors=[x for x in SECTORS if x != s])
    out = {}
    for name, kw in variants.items():
        D = build(**kw)
        p, c2, C = fit(D)
        dof = 2 * len(D['obs']) - len(p)
        sd = np.sqrt(np.diag(C))
        out[name] = dict(delta=p[:2].tolist(), delta_sd=sd[:2].tolist(), c=p[2:].tolist(), c_sd=sd[2:].tolist(), chi2=c2, dof=dof,
                         sep=float(np.hypot(*p[:2])), pa=float(np.degrees(np.arctan2(p[0], p[1])) % 360))
        print(f"{name:32s} delta=({p[0]:+.2f}+-{sd[0]:.2f}, {p[1]:+.2f}+-{sd[1]:.2f})\"  |d|={np.hypot(*p[:2]):.2f}\" PA={out[name]['pa']:5.0f}  "
              f"c={np.round(p[2:],2)}  chi2={c2:.1f}/{dof}")
    # profile chi2 map (baseline)
    D = build()
    pb, c2b, _ = fit(D)
    g = np.arange(-4.0, 4.01, 0.05)
    M = np.zeros((len(g), len(g)))
    for i, e in enumerate(g):
        for j, n in enumerate(g):
            _, c2, _ = fit(D, fix_delta=np.array([e, n]))
            M[j, i] = c2 - c2b
    EE, NN = np.meshgrid(g, g)
    rr = np.hypot(EE, NN)
    lim = {}
    for lvl, name in [(2.30, '1sigma'), (6.18, '2sigma'), (11.83, '3sigma'), (28.74, '5sigma')]:
        inside = M <= lvl
        lim[name] = dict(dchi2=lvl, max_sep=float(rr[inside].max()), min_sep=float(rr[inside].min()))
        # PA range of allowed region (for points with r > 0.3")
    print('\nprofile chi2 (baseline): allowed region radii', {k: (round(v['min_sep'], 2), round(v['max_sep'], 2)) for k, v in lim.items()})
    print('dchi2 at delta=0:', round(float(M[np.argmin(np.abs(g)), np.argmin(np.abs(g))]), 2))
    # exclusion by separation: minimum dchi2 over all PAs at each radius
    rad_bins = np.arange(0, 4.01, 0.25)
    excl = []
    for r0 in rad_bins:
        s = np.abs(rr - r0) < 0.03
        excl.append(dict(sep=float(r0), min_dchi2=float(M[s].min()) if s.any() else None))
    print('min dchi2 vs separation (any PA):', [(e['sep'], round(e['min_dchi2'], 1)) for e in excl])
    out['profile'] = dict(levels=lim, dchi2_at_zero=float(M[np.argmin(np.abs(g)), np.argmin(np.abs(g))]), min_dchi2_vs_sep=excl)
    np.savez_compressed(f'{R}/results/s06_chi2map.npz', g=g, M=M)
    to_json(out, f'{R}/results/s06_variants.json')
