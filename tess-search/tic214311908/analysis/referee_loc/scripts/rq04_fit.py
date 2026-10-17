"""Referee hypothesis fits, free-position fits and injection-recovery on the LQ window images (rq03_img.py).

Image model per sector:  d = a * PRF(apparent position) + plane (c0 + c1 x + c2 y), weighted LSQ, pixels within RFIT
px of the M3V's apparent position.  Plane-projected closed form for all hypotheses / grid points at once.
Noise: weights = 1/(formal combined variance * k^2), k = per-sector median of |plane-only residual| / formal sigma
over fake combos (variant: per-pixel inflation).  Noise realisations = combinations of fake windows (same number as
real, >= 0.7 d apart).
Hypotheses: every Gaia star within 80" of the Mira with T < 17.6, apparent position = Gaia + PM + header WCS +
referee WCS offset + lateral-colour shift; the Mira at the referee's measured apparent position (rq02).
Free fit: 0.5" sky grid co-moving with the M3V (source colour CFREE), quadratic refinement.
Injections: true source j in {M3V, pairA, pairB, Mira, G175, G177, G170}: noise combo + A*PRF(pos_j + perturbation);
perturbations: WCS N(0, SWCS px) per axis per sector; colour-curve scale common to all sectors N(1, 0.2) plus a
per-sector radial N(0, sc) with sc 0.15/0.3/0.5" for C <1.8/<2.3/>=2.3; Mira N(0, 1.0") per axis; pairB N(0,1").
PSF modes for injection: prf (same as fit) or blur (PRF convolved with Gaussian sigma 0.2 px + random 0.05 px shift).
Usage: rq04_fit.py IMGTAG CMODEL NINJ PSFMODE OUTTAG [extra options as key=value: RFIT, CFREE, INFL(pix|med), AMP(ref|own),
       MIRA(meas|loc|curve), DROP=<sector>]
"""
import sys, os, time
sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
from scipy.ndimage import gaussian_filter

imgtag, CMODEL, NINJ, PSFMODE, OUTTAG = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
opt = dict(a.split('=') for a in sys.argv[6:])
RFIT = float(opt.get('RFIT', 5.0)); CFREE = float(opt.get('CFREE', 2.76)); INFL = opt.get('INFL', 'med')
AMPMODE = opt.get('AMP', 'ref'); MIRAMODE = opt.get('MIRA', 'meas'); SWCS = float(opt.get('SWCS', 0.03))
DROP = int(opt.get('DROP', 0)); NNOISE = int(opt.get('NNOISE', 4000)); CSCALE = float(opt.get('CSCALE', 1.0))
OFFMODE = opt.get('OFF', 'ref'); OFFKEY = opt.get('OFFKEY', '')
C175 = opt.get('C175', '')
secs = [s for s in SECTORS if s != DROP]
rng = np.random.default_rng(2026)
IM = np.load(f'{ROOT}/out/img_{imgtag}.npz')
WM = json.load(open(f'{ROOT}/out/wcs_mira.json')) if os.path.exists(f'{ROOT}/out/wcs_mira.json') else {}
g = gaia()
mira_row = g[g.source_id == ID['Mira']].iloc[0]
sepm = np.hypot(g.dE - mira_row.dE, g.dN - mira_row.dN)
HYP = g[(sepm < 80) & (g['T'] < 17.6)].reset_index(drop=True)
if C175:
    HYP.loc[HYP.source_id == ID['G175'], 'C'] = float(C175)
hnames = HYP.name.values
MAIN = ['M3V', 'pairA', 'pairB', 'Mira', 'G175', 'G177', 'G170']
imain = [int(np.where(hnames == k)[0][0]) for k in MAIN]
i3 = imain[0]
m3 = g[g.source_id == ID['M3V']].iloc[0]
GE = np.arange(-30, 30.01, 0.5); GN = np.arange(-15, 15.01, 0.5)
GG = np.array([(e, n) for n in GN for e in GE])
print(f'{len(HYP)} hypotheses; grid {len(GG)} points', flush=True)


class SM:
    def __init__(self, sec):
        S = Sec(sec, cube=False)
        w = WM.get(str(sec), {})
        S.off = (np.array(w['global_variants'][OFFKEY]['off']) if OFFKEY else np.array(w['offset_px'])) if OFFMODE == 'ref' else np.array(json.load(open(
            'tic214311908/analysis/loc/out/wcs_cal.json'))[str(sec)]['offset_px'])
        self.S = S; psf = S.prf(); self.psf = psf
        xy = S.xy(HYP, CMODEL, CSCALE)
        jm = int(np.where(hnames == 'Mira')[0][0])
        mg = S.xy(HYP.iloc[[jm]], use_colour=False)[0]
        if MIRAMODE == 'meas':
            dEN = np.array(w['mira']['R3.0_blur0']['dEN'])
            xy[jm] = mg + S.Ji @ dEN
        elif MIRAMODE == 'loc':
            dEN = np.array(json.load(open('tic214311908/analysis/loc/out/mira_pos.json'))[str(sec)]['base']['dEN'])
            xy[jm] = mg + S.Ji @ dEN
        self.xy = xy
        p3 = xy[i3]
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]
        self.mask = (np.hypot(xx - p3[0], yy - p3[1]) <= RFIT).ravel()
        self.npix = int(self.mask.sum())
        Pl = np.column_stack([np.ones(S.ny * S.nx), xx.ravel() - p3[0], yy.ravel() - p3[1]])[self.mask]
        D = IM[f's{sec}_real_D'].reshape(len(IM[f's{sec}_real_D']), -1)[:, self.mask].astype(float)
        V = IM[f's{sec}_real_V'].reshape(len(D), -1)[:, self.mask].astype(float)
        FD = IM[f's{sec}_fake_D'].reshape(len(IM[f's{sec}_fake_D']), -1)[:, self.mask].astype(float)
        FV = IM[f's{sec}_fake_V'].reshape(len(FD), -1)[:, self.mask].astype(float)
        self.ftc = IM[f's{sec}_fake_tc']; self.nreal = len(D)
        self.FD, self.FV = FD, FV
        wr = 1 / V
        self.d_real = (D * wr).sum(0) / wr.sum(0); self.v_real = 1 / wr.sum(0)
        self.d_single = D; self.v_single = V
        # noise combos
        combos = [self.draw() for _ in range(NNOISE)]
        self.noise = np.array([self.comb(c)[0] for c in combos])
        vform = np.array([self.comb(c)[1] for c in combos])
        # inflation: residual of plane-only fit / formal sigma
        w0 = 1 / np.median(vform, 0)
        Q = self._qr(Pl, w0)
        r = self.noise * np.sqrt(w0) - (self.noise * np.sqrt(w0)) @ Q @ Q.T
        z = r / np.sqrt(np.median(vform, 0) * w0)
        rr = 1.4826 * np.median(np.abs(z), 0)
        self.k_med = float(np.median(rr)); self.k_pix = rr
        kk = np.sqrt(0.5 * rr ** 2 + 0.5 * self.k_med ** 2) if INFL == 'pix' else np.full(self.npix, self.k_med)
        self.sig = np.sqrt(self.v_real) * kk
        # noise realisations need their own per-combo formal variance ~ real; use the real weights (same n windows)
        self.w = 1 / self.sig ** 2
        self.sw = np.sqrt(self.w)
        self.Q = self._qr(Pl, self.w)
        # hypothesis matrix
        self.H = np.array([psf(*p).ravel()[self.mask] for p in xy])
        self.Hp = self.proj(self.H)
        # grid: sky offsets from the M3V (co-moving), colour CFREE
        dtp = S.epoch - 2016.0
        m3_now = np.array([m3.dE + m3.pmra * dtp / 1000 * 1.0, m3.dN + m3.pmdec * dtp / 1000])
        # pixel position of sky offset (relative to M3V now): use M3V apparent position without colour + Jacobian
        p3nc = S.xy(HYP.iloc[[i3]], use_colour=False)[0]
        cs = CSCALE * colour_shift(CFREE, CMODEL) * S.u
        self.gxy = p3nc[None, :] + (S.Ji @ (GG + cs[None, :]).T).T
        self.G = np.array([psf(*p).ravel()[self.mask] for p in self.gxy])
        self.Gp = self.proj(self.G)
        self.p3nc = p3nc

    def _qr(self, Pl, w):
        q, _ = np.linalg.qr(Pl * np.sqrt(w)[:, None]); return q

    def proj(self, M):
        Mw = M * self.sw[None, :]
        return Mw - (Mw @ self.Q) @ self.Q.T

    def draw(self):
        for _ in range(500):
            idx = rng.permutation(len(self.ftc)); sel = []
            for i in idx:
                if all(abs(self.ftc[i] - self.ftc[j]) >= 0.7 for j in sel): sel.append(i)
                if len(sel) == self.nreal: return sel
        raise RuntimeError

    def comb(self, sel):
        w = 1 / self.FV[sel]
        return (self.FD[sel] * w).sum(0) / w.sum(0), 1 / w.sum(0)

    def chi2(self, D, which='H'):
        """D: (K, npix) data -> chi2 (K, nhyp), amplitude (K, nhyp)."""
        Dp = self.proj(D)
        Mp = self.Hp if which == 'H' else self.Gp
        nn = (Mp ** 2).sum(1)
        pr = Dp @ Mp.T
        c0 = (Dp ** 2).sum(1)
        return c0[:, None] - pr ** 2 / nn[None, :], pr / nn[None, :], c0


def refine(c2map):
    n_e, n_n = len(GE), len(GN)
    i = int(np.argmin(c2map)); iy, ix = divmod(i, n_e)
    ix = min(max(ix, 1), n_e - 2); iy = min(max(iy, 1), n_n - 2)
    Z = c2map.reshape(n_n, n_e)[iy - 1:iy + 2, ix - 1:ix + 2]
    gx = (Z[1, 2] - Z[1, 0]) / 2; gy = (Z[2, 1] - Z[0, 1]) / 2
    hxx = Z[1, 2] - 2 * Z[1, 1] + Z[1, 0]; hyy = Z[2, 1] - 2 * Z[1, 1] + Z[0, 1]
    hxy = (Z[2, 2] - Z[2, 0] - Z[0, 2] + Z[0, 0]) / 4
    Hm = np.array([[hxx, hxy], [hxy, hyy]])
    try:
        dl = -np.linalg.solve(Hm, [gx, gy]) if np.all(np.linalg.eigvalsh(Hm) > 0) else np.zeros(2)
    except np.linalg.LinAlgError:
        dl = np.zeros(2)
    dl = np.clip(dl, -1, 1)
    cmin = Z[1, 1] - 0.5 * (np.array([gx, gy]) @ dl) if np.any(dl) else Z[1, 1]
    return np.array([GE[ix] + dl[0] * 0.5, GN[iy] + dl[1] * 0.5]), float(cmin)


t0 = time.time()
M = {s: SM(s) for s in secs}
print(f'models built {time.time()-t0:.0f}s; k_med ' + ' '.join(f'S{s}:{M[s].k_med:.2f}' for s in secs), flush=True)
res = dict(secs=secs, hyp_names=hnames.tolist(), main=MAIN, opts=dict(imgtag=imgtag, CMODEL=CMODEL, PSFMODE=PSFMODE, **opt))

# ---------------- real data
realc = {}; reala = {}; realg = {}; per = {}
for s in secs:
    c, a, c0 = M[s].chi2(M[s].d_real[None, :], 'H'); cg, ag, _ = M[s].chi2(M[s].d_real[None, :], 'G')
    realc[s] = c[0]; reala[s] = a[0]; realg[s] = cg[0]
    pos, cmin = refine(cg[0])
    dof = M[s].npix - 4
    per[s] = dict(chi2_M3V=float(c[0, i3]), dof=dof, chi2_null=float(c0[0]), amp=dict(zip(MAIN, a[0, imain].tolist())),
                  dchi2=dict(zip(MAIN, (c[0, imain] - c[0, i3]).tolist())), best=str(hnames[np.argmin(c[0])]),
                  dchi2_best_minus_M3V=float(c[0].min() - c[0, i3]), free_pos=pos, free_dchi2_vs_M3V=float(cmin - c[0, i3]),
                  amp_err_M3V=float(1 / np.sqrt((M[s].Hp[i3] ** 2).sum())), k_med=M[s].k_med)
    print(f"S{s}: chi2(M3V) {c[0,i3]:.1f}/{dof}  null {c0[0]:.1f}  A(M3V) {a[0,i3]:.1f}+-{per[s]['amp_err_M3V']:.1f}  "
          + ' '.join(f"{k} {c[0,i]-c[0,i3]:+.1f}" for k, i in zip(MAIN[1:], imain[1:])) + f"  best {per[s]['best']}  free {pos.round(1)} ({cmin-c[0,i3]:+.1f})", flush=True)
csum = sum(realc[s] for s in secs); gsum = sum(realg[s] for s in secs)
posc, cminc = refine(gsum)
res['real'] = dict(per_sector=per, comb_dchi2=dict(zip(MAIN, (csum[imain] - csum[i3]).tolist())),
                   comb_best=str(hnames[np.argmin(csum)]), comb_free_pos=posc, comb_free_dchi2_vs_M3V=float(cminc - csum[i3]),
                   comb_dchi2_all=dict(zip(hnames.tolist(), (csum - csum[i3]).tolist())))
# common amplitude (weighted) for the M3V
wA = np.array([1 / (1 / (M[s].Hp[i3] ** 2).sum()) for s in secs]); aA = np.array([reala[s][i3] for s in secs])
A_ref = float((wA * aA).sum() / wA.sum())
res['real']['A_M3V_comb'] = A_ref; res['real']['A_M3V_comb_err'] = float(1 / np.sqrt(wA.sum()))
A_own = {}
for k, i in zip(MAIN, imain):
    wk = np.array([(M[s].Hp[i] ** 2).sum() for s in secs]); ak = np.array([reala[s][i] for s in secs])
    A_own[k] = float((wk * ak).sum() / wk.sum())
res['real']['A_own'] = A_own
res['real']['gsum_map'] = (gsum - csum[i3]).astype(np.float32)
print('COMBINED real: ' + ' '.join(f"{k} {csum[i]-csum[i3]:+.1f}" for k, i in zip(MAIN, imain)) +
      f" | best {res['real']['comb_best']} | free {posc.round(2)} dchi2 {cminc-csum[i3]:+.2f} | A_M3V {A_ref:.2f}+-{res['real']['A_M3V_comb_err']:.2f}", flush=True)

# ---------------- injections
def inj_pos(s, j, r):
    S = M[s].S; p = M[s].xy[j].copy()
    p += r.normal(0, SWCS, 2)
    nm = hnames[j]
    if nm == 'Mira' or nm == 'pairB':
        p += S.Ji @ r.normal(0, 1.0, 2)
    return p


def psf_img(s, p, r):
    psf = M[s].psf
    if PSFMODE == 'prf':
        return psf(*p).ravel()[M[s].mask]
    if PSFMODE == 'blur':
        q = p + r.normal(0, 0.05, 2)
        return gaussian_filter(psf(*q), 0.2, mode='constant').ravel()[M[s].mask]
    raise ValueError


if NINJ > 0:
    inj = {}
    for k, j in zip(MAIN, imain):
        A = A_ref if AMPMODE == 'ref' else A_own[k]
        r = np.random.default_rng(abs(hash(k)) % 2 ** 32)
        cs_common = r.normal(1.0, float(opt.get('CSSD', 0.2)), NINJ)                 # common colour-curve scale per realisation
        Ck = HYP.C.values[j]
        sc = 0.5 if Ck >= 2.3 else 0.3 if Ck >= 1.8 else 0.15
        chi = np.zeros((NINJ, len(secs), len(MAIN)), np.float32); amin_other = np.zeros((NINJ, len(secs)), np.float32)
        free = np.zeros((NINJ, len(secs), 2), np.float32); freec = np.zeros((NINJ, len(secs)), np.float32)
        gsumI = np.zeros((NINJ, len(GG)), np.float32); csumI = np.zeros((NINJ, len(HYP)), np.float32)
        for si, s in enumerate(secs):
            S = M[s].S
            idx = r.integers(0, NNOISE, NINJ)
            sig = []
            base_cs = colour_shift(Ck, CMODEL) if hnames[j] != 'Mira' else 0.0
            for n in range(NINJ):
                p = inj_pos(s, j, r)
                if hnames[j] != 'Mira':
                    drad = (cs_common[n] - 1.0) * base_cs + r.normal(0, sc)
                    p = p + S.Ji @ (drad * S.u)
                sig.append(psf_img(s, p, r))
            Dn = M[s].noise[idx] + A * np.array(sig)
            c, a, _ = M[s].chi2(Dn, 'H'); cg, _, _ = M[s].chi2(Dn, 'G')
            chi[:, si] = c[:, imain]
            csumI += c; gsumI += cg
            for n in range(NINJ):
                free[n, si], freec[n, si] = refine(cg[n])
        cpos = np.array([refine(gsumI[n]) for n in range(NINJ)], dtype=object)
        inj[k] = dict(chi=chi, free=free, freec=freec, comb_free=np.array([x[0] for x in cpos]).astype(np.float32),
                      comb_freec=np.array([x[1] for x in cpos]).astype(np.float32),
                      comb_best=np.argmin(csumI, 1), comb_minother=np.array([np.min(np.delete(csumI[n], imain)) for n in range(NINJ)]),
                      comb_c=csumI[:, imain],
                      gsum_at_true=None)
        # delta-chi2 at the true star's grid point is not needed: store chi2(true hypothesis) - min over grid
        inj[k]['comb_true_minus_free'] = csumI[:, j] - inj[k]['comb_freec']
        dc = csumI[:, imain] - csumI[:, [i3]]
        q = np.percentile(dc, [0.15, 2.5, 50, 97.5, 99.85], axis=0)
        print(f'INJ {k} (A={A:.1f}) combined chi2(h)-chi2(M3V) 0.15/2.5/50/97.5/99.85%: ' +
              ' '.join(f'{h}[{q[0,ii]:.0f},{q[1,ii]:.0f},{q[2,ii]:.0f},{q[3,ii]:.0f},{q[4,ii]:.0f}]' for ii, h in enumerate(MAIN[1:], 1)) +
              f'  best=true {np.mean(inj[k]["comb_best"] == j):.3f}  ({time.time()-t0:.0f}s)', flush=True)
    np.savez_compressed(f'{ROOT}/out/inj_{OUTTAG}.npz', **{f'{k}__{kk}': v for k, d in inj.items() for kk, v in d.items() if v is not None},
                        real_comb_dchi2=np.array([res['real']['comb_dchi2'][k] for k in MAIN]),
                        real_per_sector=np.array([[per[s]['dchi2'][k] for k in MAIN] for s in secs]), GE=GE, GN=GN,
                        real_gsum=res['real']['gsum_map'])
    # summary stats
    st = {}
    for k in MAIN:
        dc = inj[k]['comb_c'] - inj[k]['comb_c'][:, [0]]
        st[k] = dict(q=np.percentile(dc, [0.15, 2.5, 16, 50, 84, 97.5, 99.85], axis=0).tolist(),
                     frac_ge_real=[float(np.mean(dc[:, ii] >= res['real']['comb_dchi2'][h])) for ii, h in enumerate(MAIN)],
                     frac_le_real=[float(np.mean(dc[:, ii] <= res['real']['comb_dchi2'][h])) for ii, h in enumerate(MAIN)],
                     best_is_true=float(np.mean(inj[k]['comb_best'] == imain[MAIN.index(k)])),
                     free_pos_mean=inj[k]['comb_free'].mean(0).tolist(), free_pos_sd=inj[k]['comb_free'].std(0).tolist(),
                     true_minus_free_q=np.percentile(inj[k]['comb_true_minus_free'], [50, 68.27, 95.45, 99.73]).tolist())
    res['inj'] = st
if int(opt.get('RING', 0)) > 0:
    NR = int(opt['RING']); ring = {}
    for rr_ in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        stat = []
        for pa in range(0, 360, 45):
            off = rr_ * np.array([np.sin(np.radians(pa)), np.cos(np.radians(pa))])
            r = np.random.default_rng(int(rr_ * 1000 + pa))
            csum_ = np.zeros((NR, len(HYP))); gsum_ = np.zeros((NR, len(GG)))
            for s in secs:
                S = M[s].S
                cs = colour_shift(CFREE, CMODEL) * S.u
                p0 = M[s].p3nc + S.Ji @ (off + cs)
                sig = np.array([psf_img(s, p0 + r.normal(0, SWCS, 2) + S.Ji @ (r.normal(0, 0.5) * S.u), r) for _ in range(NR)])
                Dn = M[s].noise[r.integers(0, NNOISE, NR)] + A_ref * sig
                c, _, _ = M[s].chi2(Dn, 'H'); cg, _, _ = M[s].chi2(Dn, 'G')
                csum_ += c; gsum_ += cg
            fmin = np.array([refine(gsum_[n])[1] for n in range(NR)])
            stat.append(csum_[:, i3] - fmin)
        ring[str(rr_)] = np.concatenate(stat)
    np.savez_compressed(f'{ROOT}/out/ring_{OUTTAG}.npz', **ring)
    m3 = None
    if NINJ > 0:
        thr = np.percentile(inj['M3V']['comb_true_minus_free'], [95.45, 99.73])
        res['ring'] = {k: dict(frac_excl_954=float(np.mean(v > thr[0])), frac_excl_9973=float(np.mean(v > thr[1])),
                               median_dchi2=float(np.median(v))) for k, v in ring.items()}
        res['ring_thresholds'] = thr.tolist()
        for k, v in res['ring'].items(): print(f"RING r={k}\": M3V excluded at 95.45% in {v['frac_excl_954']:.2f}, at 99.73% in {v['frac_excl_9973']:.2f}; median dchi2 {v['median_dchi2']:.1f}")
res['real'].pop('gsum_map')
jdump(res, f'{ROOT}/out/fit_{OUTTAG}.json')
print('done', time.time() - t0)
