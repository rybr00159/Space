"""Multi-source fits of the real difference/amplitude images with free amplitudes (item 2):
 (a) unconstrained least squares with the main candidate sources {M3V, pair (A+B), Mira (measured apparent position),
     the four faint Gaia stars within 5" of the M3V or of the pair} + constant, per sector and jointly (shared
     amplitudes across sectors = inverse-covariance-weighted mean);
 (b) non-negative least squares with ALL 93 catalogued stars within 80" (T < 17.6) + constant (+/-), per sector and
     joint (stacked design over sectors with shared amplitudes);
 (c) the Mira test from flux scaling: per-sector amplitude under the Mira hypothesis versus the Mira's TESS flux in
     that sector (a dip on the Mira scales with its flux; a dip on a constant star does not).
Output: out/joint_fits.json"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from scipy.optimize import nnls
out = {}
mp = json.load(open(f'{BASE}/out/mira_pos.json'))
SETS = {'3src': ['M3V', 'pair', 'Mira'], '4src': ['M3V', 'pair', 'Mira', '5925373272827260544'],
        '7src': ['M3V', 'pair', 'Mira', '5925373272827260544', '5925373272827262208', '5925373272827256576', '5925373272827256448']}
import itertools
for (meth, (setname, NAMES)) in itertools.product(('DIFF', 'AMP'), SETS.items()):
    per = {}; Ms = {}
    for sec in SECTORS:
        M = SectorModel(sec, meth); Ms[sec] = M
        h = M.hyp; idx = lambda n: int(np.where(h.name == n)[0][0])
        names = NAMES
        cols = []
        for n in names:
            if n == 'pair':
                # pair A + B with catalogue flux ratio (the dip could be on either; they are 0.83" apart)
                fa, fb = h.F[idx('pairA')], h.F[idx('pairB')]
                cols.append((fa * M.H[idx('pairA')] + fb * M.H[idx('pairB')]) / (fa + fb))
            else:
                cols.append(M.H[idx(n)])
        X = np.column_stack(cols + [np.ones(M.npix)])
        w = np.sqrt(M.w)
        c, *_ = np.linalg.lstsq(X * w[:, None], M.d_real * w, rcond=None)
        cov = np.linalg.inv((X * w[:, None]).T @ (X * w[:, None]))
        rr_ = (X @ c - M.d_real) * w; chi2 = float(rr_ @ rr_)
        per[sec] = dict(a=c[:-1], cov=cov[:-1, :-1], chi2=chi2, dof=M.npix - X.shape[1])
        # NNLS with all hypotheses
        XA = np.column_stack([M.H.T, np.ones(M.npix), -np.ones(M.npix)]) * w[:, None]
        cn, rn = nnls(XA, M.d_real * w)
        nz = np.where(cn[:-2] > 0.5)[0]
        per[sec]['nnls'] = {h.name[i]: float(cn[i]) for i in nz}
        per[sec]['nnls_chi2'] = float(rn ** 2)
    # joint (a): inverse-covariance weighted
    Ci = sum(np.linalg.inv(per[s]['cov']) for s in SECTORS)
    aj = np.linalg.solve(Ci, sum(np.linalg.inv(per[s]['cov']) @ per[s]['a'] for s in SECTORS))
    cj = np.linalg.inv(Ci)
    # consistency of per-sector amplitude vectors with the joint
    chi_cons = float(sum((per[s]['a'] - aj) @ np.linalg.inv(per[s]['cov']) @ (per[s]['a'] - aj) for s in SECTORS))
    corr = cj / np.sqrt(np.outer(np.diag(cj), np.diag(cj)))
    # joint (b): stacked NNLS with shared amplitudes across sectors (constant per sector)
    blocks = []; rhs = []
    nh = len(Ms[SECTORS[0]].hyp)
    for k, s in enumerate(SECTORS):
        M = Ms[s]; w = np.sqrt(M.w)
        cst = np.zeros((M.npix, 2 * len(SECTORS))); cst[:, 2 * k] = 1; cst[:, 2 * k + 1] = -1
        blocks.append(np.column_stack([M.H.T, cst]) * w[:, None]); rhs.append(M.d_real * w)
    XS = np.vstack(blocks); yS = np.concatenate(rhs)
    cs, rs = nnls(XS, yS)
    hn = Ms[SECTORS[0]].hyp.name.values
    nnls_joint = {hn[i]: float(cs[i]) for i in np.where(cs[:nh] > 0.3)[0]}
    # (c) Mira flux scaling
    am = []; sm = []; fm = []
    for s in SECTORS:
        M = Ms[s]; c2, a, b, sa = M.fit(M.d_real, 'H'); i = int(np.where(M.hyp.name == 'Mira')[0][0])
        am.append(a[i]); sm.append(sa[i]); fm.append(mp[str(s)]['base']['flux'])
        i3 = int(np.where(M.hyp.name == 'M3V')[0][0])
    am, sm, fm = map(np.array, (am, sm, fm))
    wv = 1 / sm ** 2
    cprop = (wv * am * fm).sum() / (wv * fm * fm).sum(); chi_prop = float((wv * (am - cprop * fm) ** 2).sum())
    cconst = (wv * am).sum() / wv.sum(); chi_const = float((wv * (am - cconst) ** 2).sum())
    out[f'{meth}_{setname}'] = dict(names=names, per_sector={str(s): dict(a=per[s]['a'].tolist(), err=np.sqrt(np.diag(per[s]['cov'])).tolist(),
                                                            chi2=per[s]['chi2'], dof=per[s]['dof'], nnls=per[s]['nnls'],
                                                            nnls_chi2=per[s]['nnls_chi2']) for s in SECTORS},
                     joint=dict(a=aj.tolist(), err=np.sqrt(np.diag(cj)).tolist(), corr=corr.tolist(), chi2_consistency=chi_cons,
                                dof_consistency=len(names) * (len(SECTORS) - 1)),
                     nnls_joint=nnls_joint,
                     mira_scaling=dict(amp=am.tolist(), err=sm.tolist(), mira_flux=fm.tolist(), chi2_proportional=chi_prop,
                                       chi2_constant=chi_const, dof=len(SECTORS) - 1, c_prop=cprop, depth_if_mira_ppm=cprop * 1e6))
    print(f'== {meth} {setname}')
    for s in SECTORS:
        p = per[s]
        print(f"S{s}: " + ", ".join(f"{n[:8]} {a:+5.1f}+-{e:.1f}" for n, a, e in zip(names, p['a'], np.sqrt(np.diag(p['cov'])))) +
              f"  chi2 {p['chi2']:.0f}/{p['dof']}; NNLS: " + ", ".join(f"{k[:10]} {v:.1f}" for k, v in sorted(p['nnls'].items(), key=lambda x: -x[1])))
    print('JOINT: ' + ", ".join(f"{n[:8]} {a:+5.1f}+-{e:.1f}" for n, a, e in zip(names, aj, np.sqrt(np.diag(cj)))) +
          f"; per-sector consistency chi2 {chi_cons:.1f}/{len(names)*(len(SECTORS)-1)}")
    print('   correlation M3V-pair %.2f, M3V-Mira %.2f, pair-Mira %.2f' % (corr[0, 1], corr[0, 2], corr[1, 2]) + (', M3V-260544 %.2f' % corr[0, 3] if len(names) > 3 else ''))
    print('JOINT NNLS (all 93 stars, shared amplitudes): ' + ", ".join(f"{k} {v:.1f}" for k, v in sorted(nnls_joint.items(), key=lambda x: -x[1])))
    print(f'Mira scaling: amplitudes {np.round(am,1).tolist()} +- {np.round(sm,1).tolist()} vs Mira flux {np.round(fm).tolist()}: '
          f'chi2 proportional {chi_prop:.1f}, constant {chi_const:.1f} (dof {len(SECTORS)-1}); depth if Mira {cprop*1e6:.0f} ppm')
to_json(out, f'{BASE}/out/joint_fits.json')
