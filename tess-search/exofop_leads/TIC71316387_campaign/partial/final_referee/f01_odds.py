# Final referee f01: independent recomputation of bound : chance-aligned odds for the TIC 71316387 EB.
# Own implementation (does not import or reuse synthesis/referee code); uses only data products:
#   localise TESS map (s05_baseline_map.npz), VHS injections (h04), Gaia IPD calibration (v03),
#   field Gaia tables (synthesis field_Mdwarfs_1p5deg.ecsv, astro gaia_counts_30arcmin.ecsv),
#   referee r09 L(f) profiles.
# Structure: O = [CSF * q_w] / [n_C per arcsec^2] * [E_b(Lpos)/ (int Lpos dA)] * F_LC * E
# (the separation prior enters only through E_b(Lpos), so no separate f_sep factor).
import numpy as np, json
from astropy.table import Table
from scipy.stats import fisher_exact, chi2 as chi2d
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
O = B + 'final_referee/'
rng = np.random.default_rng(1234)
d_pc = 1000 / 3.3622
out = {}

# ---------- 1. position likelihood pieces (own implementation) ----------
z = np.load(B + 'localise/v3/results/s05_baseline_map.npz'); g = z['g']; Mmap = z['M']
scale = 6.18 / 6.5417
Lt_map = np.exp(-0.5 * scale * (Mmap - Mmap.min()))     # [dN, dE]
GN, GE = np.meshgrid(g, g, indexing='ij'); R = np.hypot(GN, GE)
dg = g[1] - g[0]
# azimuthal average of TESS likelihood vs radius (bound orbits are isotropic in PA)
rb = np.arange(0, 3.0001, 0.05); rc = 0.5 * (rb[1:] + rb[:-1])
Lt_r = np.array([Lt_map[(R >= a) & (R < b)].mean() for a, b in zip(rb[:-1], rb[1:])])
# VHS miss probability: mean over 3 frames; beyond 1.2" -> recovered for dm<=3; inside 0.4" -> unknown (miss=1 at 0.3")
h4 = json.load(open(B + 'localise/v3/results/h04_vhs_inner.json'))
S = sorted({r['sep'] for r in h4[0]['inj']}); D = sorted({r['dm'] for r in h4[0]['inj']})
Rec = np.zeros((len(S), len(D)))
for fr in h4:
    for r in fr['inj']:
        Rec[S.index(r['sep']), D.index(r['dm'])] += r['frac'] / len(h4)
S = np.array(S); D = np.array(D)
def miss(s, dk):
    # bilinear in (s, dk) with extrapolation rules
    s = np.atleast_1d(s).astype(float)
    col = np.array([np.interp(dk, D, Rec[i]) for i in range(len(S))])
    rec = np.interp(s, np.r_[0.3, S, 1.5], np.r_[0.0, col, 1.0 if dk <= 3.0 else col[-1]])
    return 1 - rec
dKs = np.linspace(2.26, 2.98, 13); wK = np.exp(-0.5 * ((dKs - 2.62) / 0.18)**2); wK /= wK.sum()
Vmiss_r = sum(w * miss(rc, k) for k, w in zip(dKs, wK))
# Gaia IPD separation likelihoods
cal = json.load(open(B + 'localise/v3/results/v03_gaia_like.json'))
def emp_rates(dmset):
    edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 1.6, 2.5]; kk = []; nn = []
    for lo in edges[:-1]:
        k = n = 0
        for c in cal['all_blend']['cells']:
            if abs(c['sep'][0] - lo) < 1e-6 and tuple(c['dm']) in dmset:
                k += c['k']; n += c['n']
        kk.append(k); nn.append(n)
    return np.array(edges), np.array(kk), np.array(nn)
def ipd_fun(kind):
    if kind == 'soft':
        return lambda s: np.clip((s - 0.35) / 0.15, 0, 1) * np.clip((1.45 - s) / 0.25, 0, 1)
    if kind == 'hard':
        return lambda s: ((s >= 0.48) & (s <= 1.2)).astype(float)
    if kind == 'none':
        return lambda s: np.ones_like(np.asarray(s, float))
    dmset = {(3, 4), (4, 5)} if kind == 'emp35' else {(2, 3), (3, 4), (4, 5)}
    e, k, n = emp_rates(dmset)
    rate = (k + 0.5) / (n + 1)
    single = cal['all_blend']['single']['L']
    def f(s):
        s = np.asarray(s, float); o = np.full(s.shape, single)
        i = np.searchsorted(e, s, side='right') - 1; ok = (i >= 0) & (i < len(rate)); o[ok] = rate[i[ok]]
        o[s >= e[-1]] = rate[-1]; return o
    return f
# calibration Fisher tests (own)
for lab, dmset in (('dI3-5', {(3, 4), (4, 5)}), ('dI2-5', {(2, 3), (3, 4), (4, 5)})):
    e, k, n = emp_rates(dmset)
    kin, nin = k[:2].sum(), n[:2].sum(); kout, nout = k[2:5].sum(), n[2:5].sum()
    p = fisher_exact([[kin, nin - kin], [kout, nout - kout]], alternative='less')[1]
    out[f'ipd_calib_{lab}'] = dict(k_n_0p1_0p5=[int(kin), int(nin)], k_n_0p5_1p2=[int(kout), int(nout)], fisher_one_sided_p=float(p),
                                   per_bin=[[float(a), int(b), int(c)] for a, b, c in zip(e[:-1], k, n)])

# bound projected-separation distribution (Raghavan), eccentric orbits, random phase/orientation
N = 3_000_000
Mt = 1.03 + 0.66
logP = rng.normal(5.03, 2.28, N); P_yr = 10**logP / 365.25; a = (Mt * P_yr**2)**(1 / 3)
e_ = rng.uniform(0, 0.8, N); Mm = rng.uniform(0, 2 * np.pi, N); E = Mm.copy()
for _ in range(25):
    E -= (E - e_ * np.sin(E) - Mm) / (1 - e_ * np.cos(E))
r3 = a * (1 - e_ * np.cos(E))
nu = 2 * np.arctan2(np.sqrt(1 + e_) * np.sin(E / 2), np.sqrt(1 - e_) * np.cos(E / 2))
w = rng.uniform(0, 2 * np.pi, N); cosi = rng.uniform(-1, 1, N)
# projected separation of a point at angle (w+nu) in an orbit of inclination i
proj = r3 * np.sqrt(np.cos(w + nu)**2 + (np.sin(w + nu) * cosi)**2)
s = proj / d_pc
a_tab = np.array([1, 5, 10, 20, 32, 50, 100, 160, 1e4]); ex_tab = np.array([1.0, 0.98, 0.70, 0.55, 0.25, 0.05, 0.011, 0.001, 0.0])
Wpma = 1 - np.interp(np.log(a), np.log(a_tab), ex_tab)
# combine: ETV exclusion (lcmodel e04, dchi2>9) for a<=50 au -- applied in the 'pma+etv' variant
a_etv = np.array([1, 2, 5, 10, 14, 20, 25, 30, 40, 50, 80, 160, 1e4]); ex_etv = np.array([.955, .98, .99, .94, .84, .66, .46, .30, .10, .05, .007, 0, 0])
Wetv = 1 - np.interp(np.log(a), np.log(a_etv), ex_etv)
out['bound_fsep_0.1_1.24'] = float(((s >= 0.1) & (s <= 1.24)).mean())
Lt_s = np.interp(s, rc, Lt_r, right=0.0)
Vm_s = np.interp(s, rc, Vmiss_r, right=0.0)
# background integral over area (uniform), same azimuthally averaged ingredients -> use full 2D map for TESS
Rflat = R.ravel(); Ltflat = Lt_map.ravel()
Vm_flat = np.interp(Rflat, rc, Vmiss_r)
PF = {}
for kind in ('soft', 'hard', 'emp35', 'emp25', 'none'):
    f = ipd_fun(kind)
    for wname, W in (('noPMa', np.ones(N)), ('PMa', Wpma), ('PMa+ETV', Wpma * Wetv)):
        Eb = float((W * Lt_s * Vm_s * f(s)).mean())
        Ibg = float((Ltflat * Vm_flat * f(Rflat)).sum() * dg * dg)
        PF[f'{kind}|{wname}'] = dict(E_bound=Eb, int_bg_arcsec2=Ibg, ratio_per_arcsec2=Eb / Ibg)
out['position_terms'] = PF
# check vs synthesis convention: PF_syn = (Eb/fsep)/(Ibg/circ) -> ratio * circ / fsep
circ = np.pi * 1.24**2
out['position_factor_synthesis_convention'] = {k: v['ratio_per_arcsec2'] * circ / 0.27165 for k, v in PF.items()}

# ---------- 2. host densities (own counts) ----------
fm = Table.read(B + 'synthesis/field_Mdwarfs_1p5deg.ecsv')
gc = Table.read(B + 'astro/data/gaia_counts_30arcmin.ecsv')
def sepdeg(t):
    ra0, de0 = np.radians(330.887948), np.radians(-42.835511)
    ra = np.radians(np.array(t['ra'])); de = np.radians(np.array(t['dec']))
    return np.degrees(np.arccos(np.clip(np.sin(de0) * np.sin(de) + np.cos(de0) * np.cos(de) * np.cos(ra - ra0), -1, 1)))
def Tmag(G, br):
    return G - 0.00522555 * br**3 + 0.0891337 * br**2 - 0.633923 * br + 0.0324473
def hosts(t, rmin_deg, rmax_deg):
    sp = sepdeg(t); br = np.array(t['bp_rp'].filled(np.nan) if hasattr(t['bp_rp'], 'filled') else t['bp_rp'], float)
    G = np.array(t['phot_g_mean_mag'], float); f = 10**(-0.4 * (Tmag(G, br) - 11.358))
    m = (sp > rmin_deg) & (sp < rmax_deg) & np.isfinite(br)
    area = np.pi * (rmax_deg**2 - rmin_deg**2) * 3600**2   # arcsec^2
    return br[m], f[m], area
br1, f1, A1 = hosts(fm, 1 / 60, 1.5)
br2, f2, A2 = hosts(gc, 1 / 60, 0.5)
# r09 L(f)
r9 = json.load(open(B + 'referee_all/r09_msprior.json'))
def Lfun(setname, col):
    tab = r9[setname]; ff = np.array([x['f'] for x in tab]); c = np.array([x[col] for x in tab])
    return (lambda x: np.exp(-0.5 * (np.interp(np.log(np.clip(x, 1e-4, 1)), np.log(ff), c, left=1e3, right=1e3) - c.min()))), ff, np.array([x['dist_pc'] for x in tab])
# bound L: f at which the MS pair is at 297.4 pc, 0.15 mag MS scatter + 0.05 mag depth-scale systematic
def Lbound(setname, col):
    L, ff, dist = Lfun(setname, col)
    lnf = np.interp(np.log(d_pc), np.log(dist[::-1]), np.log(ff[::-1]))
    x = np.exp(rng.normal(lnf, np.hypot(0.15, 0.05) / 1.0857, 200000))
    return float(L(x).mean()), float(np.exp(lnf))
# colour windows: BP-RP range of host, and q windows (MS mapping BP-RP -> M1 [approx. Pecaut & Mamajek], M2 0.17-0.23)
bp_tab = np.array([1.43, 1.64, 1.84, 1.97, 2.10, 2.22, 2.32, 2.45, 2.65, 2.85, 3.05])
m_tab = np.array([0.70, 0.63, 0.57, 0.50, 0.47, 0.44, 0.40, 0.37, 0.27, 0.23, 0.20])
M_G = 1.03
def qwin(c1, c2):
    m1a, m1b = np.interp(c2, bp_tab, m_tab), np.interp(c1, bp_tab, m_tab)
    tot = ((m1a + 0.17) / M_G, min((m1b + 0.23) / M_G, 1.0))
    prim = (m1a / M_G, m1b / M_G)
    return dict(q_tot=tot, w_tot=(tot[1] - tot[0]) / 0.9, q_prim=prim, w_prim=(prim[1] - prim[0]) / 0.9)
windows = {'2.1-2.5': (2.1, 2.5, 'fm'), '1.8-2.7': (1.8, 2.7, 'fm'), '2.0-2.6': (2.0, 2.6, 'fm'), '1.5-2.8 (gc)': (1.5, 2.8, 'gc'), '2.1-2.5 (gc)': (2.1, 2.5, 'gc')}
CSF = 0.58
res = {}
for lab, (c1, c2, src) in windows.items():
    br, f, A = (br1, f1, A1) if src == 'fm' else (br2, f2, A2)
    inc = (br >= c1) & (br <= c2)
    win = inc & (f >= 0.022) & (f <= 0.065)
    n_win = int(win.sum())
    qw = qwin(c1, c2)
    o = dict(n_window=n_win, area_deg2=A / 3600**2, dens_per_arcsec2=n_win / A, q=qw)
    for setname in ('A_Mann', 'B_PM_table'):
        for col in ('total', 'dchi2_lc'):
            L, ff, dist = Lfun(setname, col)
            # H_C: sum of L over ALL colour-matched hosts with f in 0.5-30% (shape-only variant: cap at SED limit 12%)
            fmax = 0.30 if col == 'total' else 0.12
            sel = inc & (f > 0.005) & (f < fmax)
            Neff = float(L(f[sel]).sum())
            Lb, f297 = Lbound(setname, col)
            o[f'{setname}|{col}'] = dict(Neff=Neff, Lbound=Lb, f297=f297,
                                         dens_eff_per_arcsec2=Neff / A,
                                         F_LC_equiv=Lb / (Neff / max(n_win, 1)))
    res[lab] = o
out['hosts'] = res

# ---------- 3. odds grid ----------
grid = []
for lab, o in res.items():
    if o['n_window'] == 0:
        continue
    for qconv in ('w_tot', 'w_prim'):
        P_T = CSF * o['q'][qconv]
        for setname in ('A_Mann', 'B_PM_table'):
            for col in ('total', 'dchi2_lc'):
                h = o[f'{setname}|{col}']
                for pk in ('soft|PMa', 'hard|PMa', 'emp35|PMa+ETV', 'emp25|PMa+ETV'):
                    pr = PF[pk]['ratio_per_arcsec2']
                    # odds = P_T * Lbound * E_b(Lpos) / (dens_eff * int Lpos dA)
                    odds = P_T * h['Lbound'] * pr / h['dens_eff_per_arcsec2']
                    grid.append(dict(colour=lab, qconv=qconv, rel=setname, L=col, pos=pk, odds=odds))
out['grid'] = grid
og = np.array([x['odds'] for x in grid])
def summ(mask):
    v = og[mask]
    return dict(n=int(v.size), min=float(v.min()), p16=float(np.percentile(v, 16)), median=float(np.median(v)), p84=float(np.percentile(v, 84)), max=float(v.max()))
cols = np.array([x['colour'] for x in grid]); Ls = np.array([x['L'] for x in grid]); pos = np.array([x['pos'] for x in grid]); qc = np.array([x['qconv'] for x in grid])
fmwin = np.isin(cols, ['2.1-2.5', '1.8-2.7', '2.0-2.6'])
out['summary'] = {
    'all': summ(np.ones(og.size, bool)),
    'MS_L, fm windows': summ(fmwin & (Ls == 'total')),
    'MS_L, fm windows, soft/hard pos': summ(fmwin & (Ls == 'total') & np.isin(pos, ['soft|PMa', 'hard|PMa'])),
    'shape-only L, fm windows': summ(fmwin & (Ls == 'dchi2_lc')),
    'gc 1.5-2.8 window (K5-M4 hosts)': summ(cols == '1.5-2.8 (gc)'),
}
# Poisson on the smallest relevant count (2.1-2.5 window, 1.5 deg field)
n = res['2.1-2.5']['n_window']
out['poisson_2.1-2.5'] = dict(n=n, mult68=[n / (chi2d.ppf(0.84, 2 * (n + 1)) / 2), n / (chi2d.ppf(0.16, 2 * n) / 2)])
json.dump(out, open(O + 'f01_odds.json', 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ('bound_fsep_0.1_1.24', 'position_factor_synthesis_convention', 'summary', 'poisson_2.1-2.5', 'ipd_calib_dI3-5', 'ipd_calib_dI2-5')}, indent=1, default=float))
for lab, o in res.items():
    print(lab, o['n_window'], round(o['area_deg2'], 3), {k: round(v, 3) for k, v in o['q'].items() if k.startswith('w')},
          {k: (round(v['Neff'], 1), round(v['Lbound'], 3), round(v['F_LC_equiv'], 2)) for k, v in o.items() if '|' in k})
