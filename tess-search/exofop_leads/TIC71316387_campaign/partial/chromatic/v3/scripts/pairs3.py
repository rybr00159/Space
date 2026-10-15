"""Joint consistency map for main-sequence K/M pairs: TESS eclipse SHAPE (both eclipses, 200-s sectors),
TESS shallow/deep ratio, and the archival CHROMATIC depth ratios (WASP, W1, W2).

For each Mamajek pair (primary K5V..M2V, cooler secondary), with MS radii/masses, circular orbit and TESS LD:
  * scan the inclination; for each, fit free deep & shallow amplitudes to the TESS data (chunk baselines,
    fixed white-noise weights) -> best inclination, chi2, and the light fractions f1, f2 the amplitudes imply;
  * implied MS distance from each star: d_i = 297.4 pc * sqrt(f_i(MS at 297.4 pc) / f_i(fit));
  * chromatic chi2 vs the measured ratios (combine3 covariances), with predictions from the generic
    small-light-fraction model at Teff1 (deep) / Teff2 (shallow), mid geometry.
A bound pair at the G star's distance needs d1 ~ d2 ~ 297 pc (Mamajek radii scatter ~10% -> ~10% in d)."""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import tess3 as TS
import predict3 as PR
import sed3 as S

data = TS.load_all()
SECS = (68, 95, 105, 106)
T3 = json.load(open(L.V3 + '/results/tess3.json'))
C3 = json.load(open(L.V3 + '/results/combine3.json'))
P3 = json.load(open(L.V3 + '/results/predict3.json'))
shp = T3['shapes']
pre = {}
for s in SECS:
    d = data[s]
    oot = (np.abs(L.dphase(d['t'], 0)) > 0.08) & (np.abs(L.dphase(d['t'], 0.5)) > 0.08)
    sig = L.mad_std(np.diff(d['f'][oot])) / np.sqrt(2)
    a, e, r = TS.fit_sector(d, shp)
    k = np.abs(r) < 5 * L.mad_std(r)
    t, f = d['t'][k], d['f'][k]
    w = np.ones(k.sum()) / sig ** 2
    pr = L.Projector(t, w, TS.chunks(t), 1)
    pre[s] = dict(t=t, w=w, pr=pr, fp=pr(f), ex=d['ex'])
tg = np.linspace(-0.15, 0.15, 451)


def cols(t, ex, pd_, ps_):
    out = []
    for prof, phc in ((pd_, shp['d']['phc']), (ps_, shp['s']['phc'])):
        acc = np.zeros(len(t)); nsub = 4
        for o in ((np.arange(nsub) + 0.5) / nsub - 0.5) * ex:
            acc += np.interp(L.dphase(t + o, phc), tg, prof, left=0, right=0)
        out.append(acc / nsub)
    return out


def chi2_free(pd_, ps_):
    tot = 0; Ad = []; As = []; Wt = []
    for s in SECS:
        P = pre[s]
        cd, cs = cols(P['t'], P['ex'], pd_, ps_)
        Sd, Ss = P['pr'](-cd), P['pr'](-cs)
        G = np.array([[np.sum(P['w'] * Sd * Sd), np.sum(P['w'] * Sd * Ss)], [np.sum(P['w'] * Sd * Ss), np.sum(P['w'] * Ss * Ss)]])
        b = np.array([np.sum(P['w'] * Sd * P['fp']), np.sum(P['w'] * Ss * P['fp'])])
        a = np.linalg.solve(G, b)
        tot += np.sum(P['w'] * P['fp'] ** 2) - b @ a
        Ad.append(a[0]); As.append(a[1]); Wt.append(G[0, 0])
    Wt = np.array(Wt)
    return tot, float(np.sum(Wt * Ad) / Wt.sum()), float(np.sum(Wt * As) / Wt.sum())


c_tmpl = chi2_free(*[L.template(tg, 0.0, shp[k]['T'], shp[k]['p'], shp[k]['q'], t0=0.0, per=L.P) for k in 'ds'])[0]
print('empirical-template chi2 %.1f' % c_tmpl, flush=True)
tab = S.mamajek()
V_SYS, DIST, AV = 11.857, 297.4, 0.02
Mv_sys = V_SYS - 5 * np.log10(DIST / 10) - AV
FV_sys = S.ZP['V'] * 10 ** (-0.4 * Mv_sys)
gm = S.star_model(teff=5883.0, tab=tab)
rel = lambda m, b: S.band_flux(m['spec'], b) / S.band_flux(m['spec'], 'V')
meas = C3['measurements']
bgt = P3['background']
TEg = np.array(sorted(int(t) for t in bgt))


def generic(te, b, g='mid'):
    y = np.array([bgt[str(t)][g][b] for t in TEg])
    return float(np.exp(np.interp(te, TEg, np.log(y))))


def chrom_chi2(te1, te2):
    tot = 0
    for nm in ('WASP', 'NEOWISE_W1', 'NEOWISE_W2'):
        m = meas[nm]; b = m['band']
        pv = np.array([generic(te1, b), generic(te2, b)])
        dl = np.array(m['ratio']) - pv
        tot += dl @ np.linalg.solve(np.array(m['cov']), dl)
    return float(tot)


cache = {}
def comp(spt):
    if spt not in cache:
        cache[spt] = S.star_model(spt, tab=tab)
    return cache[spt]


prims = ['K5V', 'K6V', 'K7V', 'K8V', 'K9V', 'M0V', 'M0.5V', 'M1V', 'M1.5V', 'M2V']
order = [k for k in tab if k[0] in 'KM' and np.isfinite(tab[k]['R_Rsun']) and tab[k]['Teff'] < 4500]
order = sorted(order, key=lambda k: -tab[k]['Teff'])
res = {}
sod_obs = T3['s_over_d']
print('%-12s %5s %5s | inc  dchi2_shape | f1fit  f1MS  d1(pc) | f2fit  f2MS  d2(pc) | s/d obs-fit %.3f | chrom chi2' % ('pair', 'T1', 'T2', sod_obs))
for p1 in prims:
    r1 = tab[p1]
    secs = [k for k in order if 250 <= r1['Teff'] - tab[k]['Teff'] <= 1000]
    for p2 in secs:
        r2 = tab[p2]
        c1, c2 = comp(p1), comp(p2)
        F1V = S.ZP['V'] * 10 ** (-0.4 * r1['Mv']); F2V = S.ZP['V'] * 10 ** (-0.4 * r2['Mv'])
        FGV = FV_sys - F1V - F2V
        Ft = (FGV * rel(gm, 'TESS'), F1V * rel(c1, 'TESS'), F2V * rel(c2, 'TESS'))
        f1MS, f2MS = Ft[1] / sum(Ft), Ft[2] / sum(Ft)
        a = (PR.GSUN * (r1['Msun'] + r2['Msun']) * PR.MSUN * (L.P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / PR.RSUN
        R1, R2 = r1['R_Rsun'], r2['R_Rsun']
        ld1 = S.ld_quadratic(r1['Teff'], 4.7, 'TESS'); ld2 = S.ld_quadratic(r2['Teff'], 4.8, 'TESS')
        imin = np.degrees(np.arccos(min((R1 + R2) / a, 1))) + 0.1
        best = None
        for inc in np.arange(imin, 90.01, 0.5):
            k1, k2 = PR.lightcurves(R1, R2, a, inc, ld1, ld2, tg)
            if k1.max() <= 0 or k2.max() <= 0:
                continue
            c, Ad, As = chi2_free(k1, k2)          # amplitudes = light fractions directly (profiles in units of each star's flux)
            if best is None or c < best[1]:
                best = (inc, c, Ad, As)
        inc, c, f1fit, f2fit = best
        d1 = DIST * np.sqrt(f1MS / f1fit) if f1fit > 0 else np.nan
        d2 = DIST * np.sqrt(f2MS / f2fit) if f2fit > 0 else np.nan
        cc = chrom_chi2(r1['Teff'], r2['Teff'])
        key = p1 + '+' + p2
        res[key] = dict(Teff1=r1['Teff'], Teff2=r2['Teff'], inc=float(inc), dchi2_shape=float(c - c_tmpl), f1_fit=f1fit, f2_fit=f2fit,
                        f1_MS=float(f1MS), f2_MS=float(f2MS), d1=float(d1), d2=float(d2), chrom_chi2=cc,
                        lumratio_fit_over_MS=float((f2fit / f1fit) / (f2MS / f1MS)))
        print('%-12s %5.0f %5.0f | %4.1f %+7.1f | %.4f %.4f %5.0f | %.4f %.4f %5.0f | f2/f1 fit/MS %.2f | %.2f' % (
            key, r1['Teff'], r2['Teff'], inc, c - c_tmpl, f1fit, f1MS, d1, f2fit, f2MS, d2, res[key]['lumratio_fit_over_MS'], cc), flush=True)
json.dump(dict(template_chi2=c_tmpl, pairs=res), open(L.V3 + '/results/pairs3.json', 'w'), indent=1)
