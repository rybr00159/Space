"""Free-amplitude joint fits on the combined LQ images (W=0.30): per sector  d = a1 PRF(M3V) + a2 [PRF(pairA)+0.4 PRF(pairB)]/1.4
+ a3 PRF(Mira) + plane, weights from formal variance x k_med; amplitudes combined over sectors by inverse covariance.
Also the 4-source version adding the G17.5 star.  Output: out/joint.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
IM = np.load(f'{ROOT}/out/img_W0.30.npz'); WM = json.load(open(f'{ROOT}/out/wcs_mira.json')); FB = json.load(open(f'{ROOT}/out/fit_base.json'))
g = gaia(); out = {}
for nsrc in (3, 4):
    Fi = 0; Fa = 0; per = {}
    for sec in SECTORS:
        S = Sec(sec, cube=False); S.off = np.array(WM[str(sec)]['offset_px']); psf = S.prf()
        xy = S.xy(g); jm = int(np.where(g.source_id == ID['Mira'])[0][0])
        xy[jm] = S.xy(g.iloc[[jm]], use_colour=False)[0] + S.Ji @ np.array(WM[str(sec)]['mira']['R3.0_blur0']['dEN'])
        P_ = {k: psf(*xy[int(np.where(g.source_id == ID[k])[0][0])]).ravel() for k in ('M3V', 'pairA', 'pairB', 'Mira', 'G175')}
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]; p3 = xy[int(np.where(g.source_id == ID['M3V'])[0][0])]
        m = (np.hypot(xx - p3[0], yy - p3[1]) <= 5.0).ravel()
        cols = [P_['M3V'], (P_['pairA'] + 0.4 * P_['pairB']) / 1.4, P_['Mira']] + ([P_['G175']] if nsrc == 4 else [])
        X = np.column_stack([c[m] for c in cols] + [np.ones(m.sum()), (xx.ravel() - p3[0])[m], (yy.ravel() - p3[1])[m]])
        D = IM[f's{sec}_real_D'].reshape(len(IM[f's{sec}_real_D']), -1)[:, m]; V = IM[f's{sec}_real_V'].reshape(len(D), -1)[:, m]
        w = 1 / V; d = (D * w).sum(0) / w.sum(0); wc = w.sum(0) / FB['real']['per_sector'][str(sec)]['k_med'] ** 2
        Xw = X * np.sqrt(wc)[:, None]; C = np.linalg.inv(Xw.T @ Xw); a = C @ Xw.T @ (d * np.sqrt(wc))
        k = len(cols); Ci = np.linalg.inv(C[:k, :k])
        Fi = Fi + Ci; Fa = Fa + Ci @ a[:k]
        per[sec] = dict(a=a[:k].tolist(), err=np.sqrt(np.diag(C[:k, :k])).tolist())
    Cc = np.linalg.inv(Fi); ac = Cc @ Fa
    names = ['M3V', 'pair', 'Mira'] + (['G175'] if nsrc == 4 else [])
    out[f'{nsrc}src'] = dict(names=names, a=ac.tolist(), err=np.sqrt(np.diag(Cc)).tolist(),
                             corr=(Cc / np.sqrt(np.outer(np.diag(Cc), np.diag(Cc)))).round(3).tolist(), per_sector=per)
    print(nsrc, 'sources:', ', '.join(f'{n} {a:.1f}+-{e:.1f}' for n, a, e in zip(names, ac, np.sqrt(np.diag(Cc)))))
    print('   corr', np.round(Cc / np.sqrt(np.outer(np.diag(Cc), np.diag(Cc))), 2).tolist())
jdump(out, f'{ROOT}/out/joint.json')
