"""Per-transit (single-window) fits and the Mira-hypothesis test with the Mira's EMPIRICAL image as its PSF.

Single windows (LQ images, W=0.30 d): amplitude at the M3V, chi2(pairA)-chi2(M3V), chi2(Mira)-chi2(M3V), using the
plane model and per-window formal errors x the sector inflation factor (rq04 'k_med').
Empirical Mira PSF: median direct image minus all other catalogued stars (catalogue fluxes, apparent positions) minus
a constant background (median of an annulus 6-9 px from the Mira), within 4.5 px of the Mira; normalised by its sum
within 4.5 px times the PRF's encircled fraction there.  The real combined image is fitted with a*E + plane, and
1000 injections of A*E on fake-window noise give the null distribution of chi2(Mira_emp)-chi2(M3V) for a true Mira source.
Output: out/perwindow.csv, out/mira_emp.json
"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *

IM = np.load(f'{ROOT}/out/img_W0.30.npz')
WM = json.load(open(f'{ROOT}/out/wcs_mira.json'))
FB = json.load(open(f'{ROOT}/out/fit_base.json'))
g = gaia()
rows = []; emp = {}
rng = np.random.default_rng(7)
A_ref = FB['real']['A_M3V_comb']
tot_real = 0.0; tot_inj = np.zeros(1000)
for sec in SECTORS:
    S = Sec(sec); S.off = np.array(WM[str(sec)]['offset_px']); psf = S.prf()
    xy = S.xy(g)
    jm = int(np.where(g.source_id == ID['Mira'])[0][0])
    mg = S.xy(g.iloc[[jm]], use_colour=False)[0] + S.Ji @ np.array(WM[str(sec)]['mira']['R3.0_blur0']['dEN'])
    xy[jm] = mg
    pos = {k: xy[int(np.where(g.source_id == ID[k])[0][0])] for k in ('M3V', 'pairA', 'Mira')}
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    p3 = pos['M3V']
    mask = (np.hypot(xx - p3[0], yy - p3[1]) <= 5.0).ravel()
    Pl = np.column_stack([np.ones(S.ny * S.nx), xx.ravel() - p3[0], yy.ravel() - p3[1]])[mask]
    kmed = FB['real']['per_sector'][str(sec)]['k_med']
    D = IM[f's{sec}_real_D'].reshape(len(IM[f's{sec}_real_D']), -1)[:, mask].astype(float)
    V = IM[f's{sec}_real_V'].reshape(len(D), -1)[:, mask].astype(float) * kmed ** 2
    tcs = IM[f's{sec}_real_tc']
    def fitchi(d, w, img):
        X = np.column_stack([img, Pl]) * np.sqrt(w)[:, None]
        c, *_ = np.linalg.lstsq(X, d * np.sqrt(w), rcond=None)
        r = d * np.sqrt(w) - X @ c
        return float(r @ r), float(c[0]), float(np.sqrt(np.linalg.inv(X.T @ X)[0, 0]))
    imgs = {k: psf(*p).ravel()[mask] for k, p in pos.items()}
    for i in range(len(D)):
        w = 1 / V[i]
        c3, a3, e3 = fitchi(D[i], w, imgs['M3V']); cp, ap, _ = fitchi(D[i], w, imgs['pairA']); cm, am, _ = fitchi(D[i], w, imgs['Mira'])
        n = int(round((tcs[i] - T0) / P))
        rows.append(dict(sec=sec, epoch=n, tc=tcs[i], A_M3V=a3, A_err=e3, dchi2_pairA=cp - c3, dchi2_Mira=cm - c3))
    # empirical Mira PSF
    O = np.median(S.flux[S.good], 0)
    others = sum(g.F.values[j] * psf(*xy[j]) for j in range(len(g)) if j != jm and g['T'].values[j] < 17.5 and np.hypot(*(xy[j] - mg)) < 10)
    rr = np.hypot(xx - mg[0], yy - mg[1])
    bk = np.median((O - others)[(rr > 6) & (rr < 9)])
    E = np.where(rr <= 4.5, O - others - bk, 0.0)
    P0 = psf(*mg)
    E = E / E[rr <= 4.5].sum() * P0[rr <= 4.5].sum()
    E = np.where(rr <= 4.5, E, P0)                    # PRF wings outside 4.5 px
    e_img = E.ravel()[mask]
    # combined real image (inverse-variance) and noise draws
    w_all = 1 / V; dcomb = (D * w_all).sum(0) / w_all.sum(0); wc = w_all.sum(0)
    c3, a3, _ = fitchi(dcomb, wc, imgs['M3V']); ce, ae, _ = fitchi(dcomb, wc, e_img); cmp_, amp_, _ = fitchi(dcomb, wc, imgs['Mira'])
    FD = IM[f's{sec}_fake_D'].reshape(len(IM[f's{sec}_fake_D']), -1)[:, mask].astype(float)
    FV = IM[f's{sec}_fake_V'].reshape(len(FD), -1)[:, mask].astype(float) * kmed ** 2
    ftc = IM[f's{sec}_fake_tc']
    inj = []
    for k in range(1000):
        sel = []
        for ii in rng.permutation(len(ftc)):
            if all(abs(ftc[ii] - ftc[j]) >= 0.7 for j in sel): sel.append(ii)
            if len(sel) == len(D): break
        wf = 1 / FV[sel]; dn = (FD[sel] * wf).sum(0) / wf.sum(0)
        dd = dn + A_ref * e_img
        inj.append(fitchi(dd, wc, e_img)[0] - fitchi(dd, wc, imgs['M3V'])[0])
    inj = np.array(inj)
    tot_real += ce - c3; tot_inj += inj
    emp[sec] = dict(dchi2_MiraEmp_minus_M3V=ce - c3, dchi2_MiraPRF_minus_M3V=cmp_ - c3, amp_emp=ae,
                    inj_true_mira_emp_q=np.percentile(inj, [0.15, 2.5, 50, 97.5, 99.85]).tolist(),
                    emp_vs_prf_rms_frac=float(np.sqrt(np.mean(((E - P0)[rr <= 3]) ** 2)) / P0.max()))
    print(sec, {k: (np.round(v, 2) if not isinstance(v, list) else np.round(v, 1).tolist()) for k, v in emp[sec].items()}, flush=True)
emp['combined'] = dict(real=tot_real, inj_q=np.percentile(tot_inj, [0.15, 2.5, 50, 97.5, 99.85]).tolist(),
                       frac_inj_ge_real=float(np.mean(tot_inj >= tot_real)))
print('combined', emp['combined'])
pd.DataFrame(rows).to_csv(f'{ROOT}/out/perwindow.csv', index=False)
print(pd.DataFrame(rows).round(2).to_string())
jdump(emp, f'{ROOT}/out/mira_emp.json')
