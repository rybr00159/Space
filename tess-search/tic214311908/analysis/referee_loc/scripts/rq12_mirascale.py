"""Mira-hypothesis flux-scaling test: per-sector amplitude of a single PRF at the Mira's measured position (+plane),
with errors from formal variance x k_med, against the Mira's fitted TESS flux (rq02).  chi2 for 'constant' and
'proportional to Mira flux'; same for the M3V hypothesis (constant expected).  Output: out/mirascale.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
IM = np.load(f'{ROOT}/out/img_W0.30.npz'); WM = json.load(open(f'{ROOT}/out/wcs_mira.json')); FB = json.load(open(f'{ROOT}/out/fit_base.json'))
g = gaia(); res = {k: [] for k in ('Mira', 'M3V')}; F = []
for sec in SECTORS:
    S = Sec(sec, cube=False); S.off = np.array(WM[str(sec)]['offset_px']); psf = S.prf()
    xy = S.xy(g); jm = int(np.where(g.source_id == ID['Mira'])[0][0])
    xy[jm] = S.xy(g.iloc[[jm]], use_colour=False)[0] + S.Ji @ np.array(WM[str(sec)]['mira']['R3.0_blur0']['dEN'])
    p3 = xy[int(np.where(g.source_id == ID['M3V'])[0][0])]
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]; m = (np.hypot(xx - p3[0], yy - p3[1]) <= 5.0).ravel()
    D = IM[f's{sec}_real_D'].reshape(len(IM[f's{sec}_real_D']), -1)[:, m]; V = IM[f's{sec}_real_V'].reshape(len(D), -1)[:, m]
    w = 1 / V; d = (D * w).sum(0) / w.sum(0); wc = w.sum(0) / FB['real']['per_sector'][str(sec)]['k_med'] ** 2
    for k in res:
        X = np.column_stack([psf(*xy[int(np.where(g.source_id == ID[k])[0][0])]).ravel()[m], np.ones(m.sum()), (xx.ravel() - p3[0])[m], (yy.ravel() - p3[1])[m]])
        Xw = X * np.sqrt(wc)[:, None]; C = np.linalg.inv(Xw.T @ Xw); a = C @ Xw.T @ (d * np.sqrt(wc))
        res[k].append((a[0], np.sqrt(C[0, 0])))
    F.append(WM[str(sec)]['mira']['R3.0_blur0']['flux'])
F = np.array(F); out = dict(mira_flux=F.tolist())
for k, v in res.items():
    a, e = np.array(v).T; w = 1 / e ** 2
    c0 = (a * w).sum() / w.sum(); cp = (a * F * w).sum() / (F * F * w).sum()
    out[k] = dict(amp=a.tolist(), err=e.tolist(), chi2_const=float(((a - c0) ** 2 * w).sum()), chi2_prop=float(((a - cp * F) ** 2 * w).sum()))
    print(k, np.round(a, 1), np.round(e, 1), 'chi2 const %.1f prop %.1f (5 dof)' % (out[k]['chi2_const'], out[k]['chi2_prop']))
print('Mira flux', F.round())
jdump(out, f'{ROOT}/out/mirascale.json')
