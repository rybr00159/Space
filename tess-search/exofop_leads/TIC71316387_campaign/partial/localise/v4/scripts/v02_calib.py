"""v4 empirical calibration of Gaia DR3 image/astrometric diagnostics with imaging-resolved pairs and imaging singles.

Target: G 11.78, ipd_frac_multi_peak (fmp) 41, RUWE 1.026, ipd_gof_harmonic_amplitude (gofha) 0.029, phase 18.2 deg,
ipd_frac_odd_win 0, BP/RP blended in 47/47 and 48/48 transits, corrected BP/RP excess C* = +0.018, nearest Gaia source 18.2".

Definitions
  primary      brightest Gaia DR3 source within 3" (KOI/TOI) or 4" (WDS-based SOAR rows) of the catalogue position;
  companion    the brightest imaging companion at 0.08-3" (median over repeat measurements within 30% in separation);
               contrasts: SOAR/Ziegler I band only, KOIs i/880nm/692nm/LP600 (in that order of preference);
  unresolved   no other Gaia DR3 source within 3" of the primary (target: nearest at 18.2");
  single       observed (non-detection table, or KOI with Nc=0) and no companion row in any catalogue.
Diagnostic set D: 30 <= fmp <= 55, RUWE < 1.4, unresolved.  D+: D and BP-blended fraction >= 0.9.
Rates are kernel-smoothed in (log10 sep, dm): Gaussian kernels (0.12 dex, 0.6 mag) -> rate, effective n.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy.stats import beta
from v4lib import V4, to_json

O = V4 + '/speckle'
R = pd.read_csv(f'{O}/rows.csv')
X = pd.read_csv(f'{O}/xm.csv').rename(columns={'BP-RP': 'bprp', 'E(BP/RP)': 'E'})
X['bpb'] = X.NBPblend / X.o_BPmag
X['rpb'] = X.NRPblend / X.o_RPmag


def cstar(E, c):
    """Riello et al. 2021 corrected BP/RP flux excess C* and its 1-sigma scatter as a function of G."""
    c = np.asarray(c, float); E = np.asarray(E, float)
    f = np.where(c < 0.5, 1.154360 + 0.033772 * c + 0.032277 * c ** 2,
                 np.where(c < 4.0, 1.162004 + 0.011464 * c + 0.049255 * c ** 2 - 0.005879 * c ** 3, 1.057572 + 0.140537 * c))
    return E - f


def sig_cstar(G):
    return 0.0059898 + 8.817481e-12 * np.asarray(G, float) ** 7.618399


X['Cstar'] = cstar(X.E, X.bprp); X['Cstar_sig'] = X.Cstar / sig_cstar(X.Gmag)
# primary per position
wds = R.groupby('pos').cat.first().str.contains('155_235|158_48|160_7|162_41|164_58|166_139|168_28')
rmax = np.where(wds.reindex(X.pos).values, 4.0, 3.0)
Xc = X[X.angDist < rmax].sort_values(['pos', 'Gmag'])
prim = Xc.groupby('pos').first().reset_index()
# nearest other Gaia source around each primary (from the 5" list) -> 'unresolved' if none within 3"
nn = []
for p, g in X.groupby('pos'):
    pr = prim[prim.pos == p]
    if not len(pr):
        continue
    pr = pr.iloc[0]
    d = np.hypot((g.RAdeg - pr.RAdeg) * np.cos(np.radians(pr.DEdeg)), g.DEdeg - pr.DEdeg) * 3600
    o = g.Source.values != pr.Source
    nn.append((p, float(d[o].min()) if o.any() else 99.0, float((g.Gmag.values[o] - pr.Gmag)[np.argmin(d[o])]) if o.any() else np.nan))
nn = pd.DataFrame(nn, columns=['pos', 'dnn', 'dG_nn'])
prim = prim.merge(nn, on='pos', how='left')
prim['unres'] = prim.dnn > 3.0
keep = ['pos', 'Source', 'Gmag', 'bprp', 'Plx', 'RUWE', 'IPDfmp', 'IPDgofha', 'IPDgofhp', 'IPDfow', 'epsi', 'sepsi', 'bpb', 'rpb', 'Cstar',
        'Cstar_sig', 'dnn', 'dG_nn', 'unres', 'Dup', 'Nper']
R = R.merge(prim[keep], on='pos', how='inner')
soar = R.cat.str.contains('155_235|158_48|160_7|162_41|164_58|166_139|168_28|159_19|162_192')
R['dmu'] = np.where(soar & ~R.band.fillna('').str.strip().isin(['I']), np.nan, R.dm)
res = R[(R.kind == 'res') & (R.sep > 0.03) & (R.sep < 6)]
pairs = []
for sid, g in res.groupby('Source'):
    g2 = g[(g.sep >= 0.08) & (g.sep <= 3.0) & np.isfinite(g.dmu)].sort_values('dmu')
    if not len(g2):
        continue
    b = g2.iloc[0]
    same = g2[np.abs(g2.sep - b.sep) < 0.3 * b.sep + 0.05]
    pa2 = np.radians(2 * same.pa.dropna().values)
    pairs.append(dict(Source=sid, sep=float(same.sep.median()), dm=float(same.dmu.median()), nmeas=len(same),
                      pa=float((np.degrees(np.arctan2(np.sin(pa2).mean(), np.cos(pa2).mean())) / 2) % 180) if len(pa2) else np.nan,
                      cats=';'.join(sorted(set(same.cat))), n_comp=int(g.sep.round(1).nunique()), min_sep_any=float(g.sep.min())))
pairs = pd.DataFrame(pairs).merge(prim.drop_duplicates('Source')[keep[1:]], on='Source', how='left')
resolved_ids = set(res.Source)
sing = R[(R.kind == 'unres') & ~R.Source.isin(resolved_ids)].drop_duplicates('Source')
print(f'pairs (0.08-3", contrast known): {len(pairs)}; singles: {len(sing)}')
pairs.to_csv(f'{O}/pairs.csv', index=False); sing.to_csv(f'{O}/singles.csv', index=False)


def jeff(k, n):
    return [float(beta.ppf(0.16, k + 0.5, n - k + 0.5)), float(beta.ppf(0.84, k + 0.5, n - k + 0.5))] if n else [None, None]


def D(d, plus=False, fmp=(30, 55)):
    m = (d.IPDfmp >= fmp[0]) & (d.IPDfmp <= fmp[1]) & (d.RUWE < 1.4) & d.unres.astype(bool)
    if plus:
        m &= d.bpb >= 0.9
    return m


out = dict(n_pairs=len(pairs), n_singles=len(sing))
# ---------------- singles
for lab, sel in [('G10.5-13', (sing.Gmag > 10.5) & (sing.Gmag < 13)), ('G9-14', (sing.Gmag > 9) & (sing.Gmag < 14))]:
    s = sing[sel]
    r = dict(n=int(len(s)))
    for nm, m in [('fmp30-55', (s.IPDfmp >= 30) & (s.IPDfmp <= 55)), ('fmp>=20', s.IPDfmp >= 20), ('D', D(s)), ('D+', D(s, True)),
                  ('bpb>=0.9', s.bpb >= 0.9), ('bpb>=0.9&unres', (s.bpb >= 0.9) & s.unres), ('Cstar_sig>2', s.Cstar_sig > 2),
                  ('unres', s.unres)]:
        k = int(m.sum()); r[nm] = dict(k=k, frac=k / max(len(s), 1), ci68=jeff(k, len(s)))
    r['fmp_q50_95_99'] = [float(np.percentile(s.IPDfmp, q)) for q in (50, 95, 99)]
    r['gofha_q50_84_95'] = [float(np.nanpercentile(s.IPDgofha, q)) for q in (50, 84, 95)]
    r['D_list'] = s[D(s)][['Source', 'Gmag', 'IPDfmp', 'RUWE', 'bpb', 'Cstar_sig', 'cat', 'lim015', 'lim1']].round(3).to_dict('records')
    out[f'singles_{lab}'] = r
    print(f"\nsingles {lab}: n={len(s)}  fmp30-55 {r['fmp30-55']['k']}  D {r['D']['k']}  D+ {r['D+']['k']}  bpb>=0.9 {r['bpb>=0.9']['k']} "
          f"(unres {r['bpb>=0.9&unres']['k']})  C*>2sig {r['Cstar_sig>2']['k']}  fmp q50/95/99 {r['fmp_q50_95_99']}")
    for x in r['D_list']:
        print('    D single:', x)
# ---------------- pairs: like the target
G = (pairs.Gmag > 9.5) & (pairs.Gmag < 13.5)
for lab, m in [('G9.5-13.5 & D', G & D(pairs)), ('G9.5-13.5 & D+', G & D(pairs, True)), ('allG & D', D(pairs)), ('allG & D+', D(pairs, True))]:
    p = pairs[m]
    out[f'pairs_{lab}'] = dict(n=int(len(p)), sep_q=[float(np.percentile(p.sep, q)) for q in (5, 16, 50, 84, 95)] if len(p) else None,
                               dm_q=[float(np.percentile(p.dm, q)) for q in (5, 16, 50, 84, 95)] if len(p) else None,
                               rows=p[['Source', 'Gmag', 'sep', 'dm', 'pa', 'IPDfmp', 'RUWE', 'IPDgofha', 'IPDgofhp', 'bpb', 'Cstar_sig', 'cats']].round(3).to_dict('records'))
    print(f"\npairs {lab}: n={len(p)} sep q5/16/50/84/95 {np.round(out[f'pairs_{lab}']['sep_q'], 2) if len(p) else None}  "
          f"dm {np.round(out[f'pairs_{lab}']['dm_q'], 2) if len(p) else None}")
    for x in out[f'pairs_{lab}']['rows'][:40]:
        print('   ', x)
# ---------------- kernel-smoothed rate maps
P9 = pairs[G].copy()
ls = np.log10(P9.sep.values); dmv = P9.dm.values
ind = {'D': D(P9).values, 'D+': D(P9, True).values}
sg = np.arange(np.log10(0.1), np.log10(3.0) + 1e-9, 0.05); dg = np.arange(0.0, 7.01, 0.25)
maps = {}
for nm, y in ind.items():
    Rm = np.zeros((len(dg), len(sg))); Nm = np.zeros_like(Rm)
    for i, d0 in enumerate(dg):
        for j, s0 in enumerate(sg):
            w = np.exp(-0.5 * ((ls - s0) / 0.12) ** 2 - 0.5 * ((dmv - d0) / 0.6) ** 2)
            Nm[i, j] = w.sum(); Rm[i, j] = ((w * y).sum() + 0.5 * 0.01) / (w.sum() + 0.01)
    maps[nm] = (Rm, Nm)
# posterior for flat priors in log sep (0.1-3") and dm (0-7), cells with n_eff >= 3 only
post = {}
for nm, (Rm, Nm) in maps.items():
    L = np.where(Nm >= 3, Rm, 0.0)
    pp = L / L.sum()
    ms = pp.sum(0); md = pp.sum(1)
    cs = np.cumsum(ms); cd = np.cumsum(md)
    q = lambda c, x, v: float(np.interp(v, c, x))
    post[nm] = dict(sep_q16_50_84=[10 ** q(cs, sg, v) for v in (0.16, 0.5, 0.84)], sep_q025_975=[10 ** q(cs, sg, v) for v in (0.025, 0.975)],
                    dm_q16_50_84=[q(cd, dg, v) for v in (0.16, 0.5, 0.84)], dm_q025_975=[q(cd, dg, v) for v in (0.025, 0.975)],
                    max_rate=float(Rm[Nm >= 3].max()), argmax_sep=float(10 ** sg[np.unravel_index(np.argmax(np.where(Nm >= 3, Rm, 0)), Rm.shape)[1]]),
                    argmax_dm=float(dg[np.unravel_index(np.argmax(np.where(Nm >= 3, Rm, 0)), Rm.shape)[0]]))
    print(f"\nkernel posterior ({nm}, flat in log sep & dm, n_eff>=3): sep 16/50/84% = {np.round(post[nm]['sep_q16_50_84'], 2)} "
          f"(95%: {np.round(post[nm]['sep_q025_975'], 2)}); dm 16/50/84% = {np.round(post[nm]['dm_q16_50_84'], 2)} "
          f"(95%: {np.round(post[nm]['dm_q025_975'], 2)}); max rate {post[nm]['max_rate']:.3f} at sep {post[nm]['argmax_sep']:.2f} dm {post[nm]['argmax_dm']:.2f}")
out['kernel_posterior'] = post
# coarse table of k/n (for transparency)
tab = []
SE = [0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0]; DM = [0, 1.5, 2.5, 3.5, 4.5, 7]
print('\ncounts k(D)/n  [k(D+)]  for G 9.5-13.5 pairs:')
print('            ' + ''.join(f'dm {a}-{b}'.ljust(16) for a, b in zip(DM[:-1], DM[1:])))
for s0, s1 in zip(SE[:-1], SE[1:]):
    line = f'{s0:4.2f}-{s1:4.2f}  '
    for d0, d1 in zip(DM[:-1], DM[1:]):
        m = (P9.sep >= s0) & (P9.sep < s1) & (P9.dm >= d0) & (P9.dm < d1)
        k = int((m & D(P9)).sum()); kp = int((m & D(P9, True)).sum()); n = int(m.sum())
        f30 = int((m & (P9.IPDfmp >= 30)).sum()); fr = float((m & ~P9.unres).sum() / max(n, 1))
        tab.append(dict(sep=[s0, s1], dm=[d0, d1], n=n, kD=k, kDplus=kp, k_fmp30=f30, frac_gaia_resolved=fr,
                        fmp_med=float(P9.IPDfmp[m].median()) if n else None, ruwe_med=float(P9.RUWE[m].median()) if n else None,
                        bpb_med=float(P9.bpb[m & P9.unres].median()) if (m & P9.unres).sum() else None))
        line += f'{k:2d}/{n:<3d}[{kp:2d}]'.ljust(16)
    print(line)
out['count_table'] = tab
# ---------------- BP blending vs separation for Gaia-unresolved pairs; C*
print('\nBP-blended fraction >= 0.9 among Gaia-UNRESOLVED G 9.5-13.5 pairs, by separation (any dm):')
bl = []
for s0, s1 in zip(SE[:-1], SE[1:]):
    m = (P9.sep >= s0) & (P9.sep < s1) & P9.unres
    k = int((m & (P9.bpb >= 0.9)).sum()); n = int(m.sum())
    kc = int((m & (P9.Cstar_sig > 2)).sum())
    bl.append(dict(sep=[s0, s1], n=n, k_bpb09=k, frac=k / max(n, 1), k_cstar2=kc))
    print(f'   {s0:.2f}-{s1:.2f}": {k}/{n}  (C*>2 sigma: {kc}/{n})')
out['bp_blend_vs_sep'] = bl
# ---------------- harmonic phase vs companion PA
hp = pairs[(pairs.sep > 0.1) & (pairs.sep < 2.0) & pairs.unres & np.isfinite(pairs.pa) & np.isfinite(pairs.IPDgofhp)]
for lab, m in [('gofha>0.05', hp.IPDgofha > 0.05), ('0.02<gofha<0.04 & fmp>=20', (hp.IPDgofha > 0.02) & (hp.IPDgofha < 0.04) & (hp.IPDfmp >= 20)),
               ('0.02<gofha<0.04 & fmp 30-55 & ruwe<1.4', (hp.IPDgofha > 0.02) & (hp.IPDgofha < 0.04) & (hp.IPDfmp >= 30) & (hp.IPDfmp <= 55) & (hp.RUWE < 1.4))]:
    h = hp[m]
    d = ((h.pa - h.IPDgofhp) % 180).values
    a = np.radians(2 * d)
    n = len(h)
    Rb = float(np.hypot(np.cos(a).mean(), np.sin(a).mean())) if n else np.nan
    mu = float((np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())) / 2) % 180) if n else np.nan
    z = n * Rb ** 2
    pR = float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n) - (24 * z - 132 * z ** 2 + 76 * z ** 3 - 9 * z ** 4) / (288 * n ** 2))) if n > 5 else None
    within = float(np.mean(np.abs(((d - mu + 90) % 180) - 90) < 30)) if n else None
    out[f'phase_{lab}'] = dict(n=n, mean_PA_minus_phase=mu, Rbar=Rb, rayleigh_p=pR, frac_within_30=within,
                               hist_30deg=np.histogram(d, bins=np.arange(0, 181, 30))[0].tolist())
    print(f'phase test [{lab}]: n={n} <(PA - phase) mod 180> = {mu:.1f} deg  Rbar={Rb:.2f}  p={pR}  frac within +-30 = {within}')
to_json(out, f'{V4}/results/v02_calib.json')
