"""Out-of-eclipse variability at P_EB and P_EB/2.
Data: eclipses masked (|dt| < 0.06 d of either conjunction; T14/2 = 0.032 d), 5-sigma clip vs 0.5-d running median,
segments split at gaps > 0.5 d, per-segment Legendre trend of degree ceil(span/2.5 d) (>=2) fitted SIMULTANEOUSLY with
    a1 cos(2 pi phi) + b1 sin(2 pi phi) + a2 cos(4 pi phi) + b2 sin(4 pi phi),   phi = 0 at primary eclipse.
Physical meaning: ellipsoidal  -> a2 = -A_ell (minima at both conjunctions);  reflection -> a1 < 0 (max at phase 0.5 for the
heated/brighter face of star 2... sign depends on which star dominates);  beaming -> b1;  spots on a synchronised star -> a1,b1
(and a2,b2) with phases that drift between sectors.
Errors: (i) fake-period null (400 periods 0.8-3.0 d avoiding P, P/2, 2P, 3P/2), (ii) sector-to-sector scatter.
Products: default (PDCSAP S1/S28/S68 + own TESScut S95/S105/S106), SAP S1/28/68 + TC, TC all six."""
import numpy as np, sys, json, glob
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import lclib as L
from astropy.io import fits
from numpy.polynomial import legendre as LG
import pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
P, T0 = L.P0, L.T00
def spoc(sec, col):
    f = glob.glob(f'{L.RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
    d = fits.open(f)[1].data
    t = np.array(d['TIME'], float); fl = np.array(d[col], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0)
    return dict(t=t[m], f=fl[m] / np.median(fl[m]), sector=sec)
def prep(lc, mask_hw=0.06):
    t, f = lc['t'], lc['f']
    ph = ((t - T0) / P) % 1
    dd = np.minimum(np.minimum(ph, 1 - ph), np.abs(ph - 0.5)) * P
    m = dd > mask_hw; t, f = t[m], f[m]
    rm = pd.Series(f, index=pd.to_datetime(t, unit='D')).rolling('0.5D', center=True, min_periods=5).median().values
    r = f - rm; s = 1.4826 * np.nanmedian(np.abs(r)); k = ~(np.abs(r) > 5 * s); t, f = t[k], f[k]
    br = np.r_[0, np.where(np.diff(t) > 0.5)[0] + 1, t.size]
    seg = np.zeros(t.size, int)
    for i in range(len(br) - 1): seg[br[i]:br[i + 1]] = i
    sig = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
    return dict(t=t, f=f, seg=seg, sig=sig)
def trend_cols(t, seg, dper=2.5):
    cols = []
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < 20: continue
        ts = t[m]; x = 2 * (ts - ts.min()) / (ts.max() - ts.min()) - 1
        deg = int(max(2, np.ceil((ts.max() - ts.min()) / dper)))
        for j in range(deg + 1):
            c = np.zeros(t.size); c[m] = LG.legval(x, np.eye(deg + 1)[j]); cols.append(c)
    return np.array(cols).T
def harm(t, per):
    ph = 2 * np.pi * (t - T0) / per
    return np.vstack([np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)]).T
class Joint:
    """trend columns projected out once; harmonic fit for any period"""
    def __init__(self, ds):
        A = []; y = []; w = []
        for d in ds:
            A.append(trend_cols(d['t'], d['seg'])); y.append(d['f'] - 1); w.append(np.full(d['t'].size, 1 / d['sig']))
        n = sum(a.shape[0] for a in A); m = sum(a.shape[1] for a in A)
        M = np.zeros((n, m)); i0 = j0 = 0
        for a in A:
            M[i0:i0 + a.shape[0], j0:j0 + a.shape[1]] = a; i0 += a.shape[0]; j0 += a.shape[1]
        self.w = np.concatenate(w); self.t = np.concatenate([d['t'] for d in ds]); yv = np.concatenate(y)
        self.Q, _ = np.linalg.qr(M * self.w[:, None])
        yw = yv * self.w; self.yp = yw - self.Q @ (self.Q.T @ yw); self.y = yv; self.M = M
    def fit(self, per):
        H = harm(self.t, per) * self.w[:, None]; Hp = H - self.Q @ (self.Q.T @ H)
        c, *_ = np.linalg.lstsq(Hp, self.yp, rcond=None)
        return c, np.sqrt(np.diag(np.linalg.inv(Hp.T @ Hp)))
rng = np.random.default_rng(42)
fake = []
while len(fake) < 400:
    pf = rng.uniform(0.8, 3.0)
    if min(abs(pf / x - 1) for x in (P, P / 2, 2 * P, 1.5 * P, P / 3)) > 0.03: fake.append(pf)
fake = np.array(fake)
prod = {}
for s in (1, 28, 68):
    prod[(s, 'PDCSAP')] = prep(spoc(s, 'PDCSAP_FLUX')); prod[(s, 'SAP')] = prep(spoc(s, 'SAP_FLUX'))
for s in (1, 28, 68, 95, 105, 106):
    lc = L.load_tc(s); prod[(s, 'TC')] = prep(lc)
COMBOS = {'default': [(1, 'PDCSAP'), (28, 'PDCSAP'), (68, 'PDCSAP'), (95, 'TC'), (105, 'TC'), (106, 'TC')],
          'sap+tc': [(1, 'SAP'), (28, 'SAP'), (68, 'SAP'), (95, 'TC'), (105, 'TC'), (106, 'TC')],
          'tc_all': [(s, 'TC') for s in (1, 28, 68, 95, 105, 106)]}
out = {}
for cname, keys in COMBOS.items():
    J = Joint([prod[k] for k in keys]); c, e = J.fit(P)
    nul = np.array([J.fit(pf)[0] for pf in fake]); esc = 1.4826 * np.median(np.abs(nul - np.median(nul, 0)), 0)
    A1 = np.hypot(c[0], c[1]); A2 = np.hypot(c[2], c[3])
    n1 = np.hypot(nul[:, 0], nul[:, 1]); n2 = np.hypot(nul[:, 2], nul[:, 3])
    per = {}
    for k in keys:
        Jk = Joint([prod[k]]); ck, ek = Jk.fit(P)
        nk = np.array([Jk.fit(pf)[0] for pf in fake[:200]]); ek_emp = 1.4826 * np.median(np.abs(nk - np.median(nk, 0)), 0)
        per[f'S{k[0]}_{k[1]}'] = dict(a1=ck[0] * 1e6, b1=ck[1] * 1e6, a2=ck[2] * 1e6, b2=ck[3] * 1e6, emp_err=(ek_emp * 1e6).tolist(),
                                     formal_err=(ek * 1e6).tolist(), A1=float(np.hypot(ck[0], ck[1]) * 1e6),
                                     phase_max_deg=float(np.degrees(np.arctan2(ck[1], ck[0])) % 360),
                                     A1_snr_emp=float(np.hypot(ck[0], ck[1]) / np.mean(ek_emp[:2])),
                                     p_A1=float(np.mean(np.hypot(nk[:, 0], nk[:, 1]) >= np.hypot(ck[0], ck[1]))))
    # ellipsoidal: -a2.  (i) joint fit with fake-period empirical error; (ii) inverse-variance mean of per-sector -a2 with
    # error = max(propagated empirical, sector scatter / sqrt(N))
    a2s = np.array([-per[f'S{k[0]}_{k[1]}']['a2'] for k in keys]); ea2 = np.array([per[f'S{k[0]}_{k[1]}']['emp_err'][2] for k in keys])
    wts = 1 / ea2 ** 2; mu = np.sum(a2s * wts) / wts.sum(); e_prop = 1 / np.sqrt(wts.sum())
    e_scat = np.std(a2s, ddof=1) / np.sqrt(len(a2s)); chi_c = float(np.sum((a2s - mu) ** 2 * wts))
    e_sec = max(e_prop, e_scat)
    res = dict(a1=c[0] * 1e6, b1=c[1] * 1e6, a2=c[2] * 1e6, b2=c[3] * 1e6, formal_err=(e * 1e6).tolist(), emp_err=(esc * 1e6).tolist(),
               A1=A1 * 1e6, p_A1=float(np.mean(n1 >= A1)), A2=A2 * 1e6, p_A2=float(np.mean(n2 >= A2)),
               Aell_joint=-c[2] * 1e6, Aell_joint_err=esc[2] * 1e6, Aell_joint_UL95=(max(-c[2], 0) + 1.645 * esc[2]) * 1e6,
               Aell_sectors=float(mu), Aell_sectors_err=float(e_sec), Aell_sectors_chi2=chi_c, Aell_sectors_UL95=float(max(mu, 0) + 1.645 * e_sec),
               beaming_b1=c[1] * 1e6, per_sector=per)
    out[cname] = res
    print(f"[{cname}] joint: a1={c[0]*1e6:+.1f} b1={c[1]*1e6:+.1f} a2={c[2]*1e6:+.1f} b2={c[3]*1e6:+.1f} ppm (emp err {np.round(esc*1e6,1)}) "
          f"A1={A1*1e6:.0f} (p={res['p_A1']:.3f}) A2={A2*1e6:.0f} (p={res['p_A2']:.3f})")
    print(f"   A_ell = -a2: joint {-c[2]*1e6:+.1f} +- {esc[2]*1e6:.1f} (UL95 {res['Aell_joint_UL95']:.0f});  sector mean {mu:+.1f} +- {e_sec:.1f} "
          f"(chi2 {chi_c:.1f}/{len(a2s)-1}; UL95 {res['Aell_sectors_UL95']:.0f}) ppm")
    for k in keys:
        q = per[f'S{k[0]}_{k[1]}']
        print(f"   S{k[0]:3d} {k[1]:6s}: A1={q['A1']:5.0f} ppm (emp sd {np.mean(q['emp_err'][:2]):4.0f}, p={q['p_A1']:.3f}) phase of max {q['phase_max_deg']:5.0f} deg | "
              f"a2={q['a2']:+6.0f} b2={q['b2']:+6.0f} (emp sd {q['emp_err'][2]:.0f})")
json.dump(out, open(f'{FD}/results/ooe.json', 'w'), indent=1, default=float)
# Lomb-Scargle per sector (default products) to show where the strongest periodicity lies
from astropy.timeseries import LombScargle
freq = np.linspace(1 / 8.0, 1 / 0.3, 30000)
fig, axs = plt.subplots(7, 1, figsize=(8, 13))
ls_out = {}
for ax, k in zip(axs, COMBOS['default']):
    d = prod[k]; A = trend_cols(d['t'], d['seg'], dper=4.0)
    c, *_ = np.linalg.lstsq(A, d['f'] - 1, rcond=None); r = d['f'] - 1 - A @ c
    pw = LombScargle(d['t'], r).power(freq); j = np.argmax(pw)
    ls_out[f'S{k[0]}'] = dict(peak_period=float(1 / freq[j]), peak_power=float(pw[j]), frac_offset_from_PEB=float(1 / freq[j] / P - 1))
    print(f'LS S{k[0]}: highest peak 0.3-8 d at {1/freq[j]:.4f} d (P_EB = {P:.4f}); power {pw[j]:.3f}')
    ax.plot(1 / freq, pw, 'k-', lw=0.6); ax.axvline(P, color='r', ls=':'); ax.axvline(P / 2, color='b', ls=':'); ax.set_xscale('log')
    ax.set_ylabel(f'S{k[0]}'); ax.set_xlim(0.3, 8)
# spot-modulation phase and amplitude per sector (default products): phase drift => spots, not ellipsoidal/reflection
pp = out['default']['per_sector']; ks = list(pp.keys())
tmid = {1: 1338, 28: 2073, 68: 3164, 95: 3890, 105: 4220, 106: 4248}
xs = [tmid[int(k.split('_')[0][1:])] for k in ks]
axs[-1].errorbar(xs, [pp[k]['phase_max_deg'] for k in ks], [np.degrees(np.mean(pp[k]['emp_err'][:2]) / pp[k]['A1']) for k in ks], fmt='o', color='k')
for x_, k in zip(xs, ks): axs[-1].annotate(f"{pp[k]['A1']:.0f} ppm", (x_, pp[k]['phase_max_deg']), fontsize=7, xytext=(4, 4), textcoords='offset points')
axs[-1].set_xlabel('BTJD'); axs[-1].set_ylabel('phase of max (deg)')
axs[-1].set_title('P_EB modulation: phase of maximum and amplitude per sector', fontsize=9)
plt.tight_layout(); plt.savefig(f'{FD}/plots/ooe.png', dpi=100)
json.dump(dict(ls=ls_out), open(f'{FD}/results/ooe_ls.json', 'w'), indent=1)
