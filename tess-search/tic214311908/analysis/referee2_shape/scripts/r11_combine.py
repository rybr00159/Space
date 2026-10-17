"""r11: combined evidence against 'white dwarf totally occulted by the M3V' as a function of M_WD:
  LC + secondary + M-dwarf priors (r10, cool branch; secondary = own INV value or REF r18 value)
  + Gaia RUWE (r06, NSL scan pattern, edge-on orbit; sigma_eff 0.18 (plx), 0.22, 0.25, 0.35 mas)
  + Gaia XP blue index (r08, red-only regression; WD Teff solved from R_WD and the 1.57% TESS share using the
    target's XP spectrum as the M3V flux and a blackbody WD)
  + Gaia RVS combined-CCF error (r13, most favourable 5th-percentile smearing, w = 16 km/s, C0 = 0.3) - listed
    separately because its error model is approximate.
Output: out/wd_combined.json"""
import json, numpy as np, pandas as pd
RD = 'tic214311908/analysis/referee2_shape'
B = RD.rsplit('/', 1)[0]
xp = json.load(open(f'{RD}/out/xp_test.json'))
rv = json.load(open(f'{RD}/out/rv_smear.json'))
X = pd.read_csv(f'{RD}/data/xp_sampled.csv')
wl = np.arange(336, 1021, 2.0)
ft = np.array(json.loads(X[X.source_id == 5925373272850784768].flux.iloc[0]))
filt = np.loadtxt(f'{B}/shape/data/filt_TESS.Red.dat'); tw = np.interp(wl*10, filt[:, 0], filt[:, 1], left=0, right=0)
d_m = 64.31*3.0857e16; RSUN_m = 6.957e8
def planck_lam(T):   # W m^-2 m^-1 sr^-1 -> per nm
    lam = wl*1e-9; return 2*6.626e-34*2.998e8**2/lam**5/(np.exp(1.4388e-2/(lam*T)) - 1)*1e-9
def teff_needed(RWD, fA):
    target = fA/(1 - fA)*np.sum(ft*wl*tw)
    Ts = np.linspace(3000, 30000, 2701)
    vals = np.array([np.sum(np.pi*planck_lam(T)*(RWD*RSUN_m/d_m)**2*wl*tw) for T in Ts])
    return float(np.interp(target, vals, Ts))
xpT = np.array([float(k) for k in xp['consistent']['red_only']['pred']]); xpD = np.array([v['dchi2_wd_minus_planet'] for v in xp['consistent']['red_only']['pred'].values()])
rvm = {float(k): v['p05']['m2lnL'] for k, v in rv['m2lnL_wd_vs_planet']['w16_C0.3'].items()}
out = {}
for ds in ('inv', 'ref'):
    W = json.load(open(f'{RD}/out/wd_total_{ds}.json'))
    rus = {}
    for sm in ('plx', '0.22', '0.25', '0.35'):
        try: rus[sm] = json.load(open(f'{RD}/out/ruwe_nsl_{sm}.json'))['masses']
        except FileNotFoundError: pass
    rows = []
    for r in W['rows']:
        if r['branch'] != 'cool': continue
        M = r['M_WD']
        T = teff_needed(r['R_WD'], r['fA'])
        xpd = float(np.interp(T, xpT, xpD)) if T <= xpT.max() else float(xpD[-1])
        row = dict(M_WD=M, sec_meas=r['sec_meas'], lc_sec=r['dchi2_total_vs_planet'], lc_only=r['dchi2_lc_vs_planet'], secondary_ppt=r['secondary_ppt'], R_WD=r['R_WD'], Teff_needed=T, xp=xpd,
                   rv=float(np.interp(M, sorted(rvm), [rvm[k] for k in sorted(rvm)])))
        for sm, R in rus.items():
            mm = min(R, key=lambda q: abs(float(q) - M))
            row[f'ruwe_{sm}'] = R[mm]['nsl']['m2lnL_density']
            row[f'total_{sm}'] = row['lc_sec'] + row[f'ruwe_{sm}']
            row[f'total_{sm}_xp'] = row[f'total_{sm}'] + max(xpd, 0)    # XP only counted when it disfavours the WD
        rows.append(row)
    D = pd.DataFrame(rows)
    res = dict(table=D.round(3).to_dict('records'))
    for sec in ('own', 'ref_r18'):
        q = D[D.sec_meas == sec]
        res[f'min_{sec}'] = {c: dict(value=float(q[c].min()), at_M=float(q.loc[q[c].idxmin(), 'M_WD'])) for c in q.columns if c.startswith('total') or c == 'lc_sec'}
    out[ds] = res
    print(ds)
    print(D[D.sec_meas == 'own'].round(2).to_string())
    print({k: (round(v['value'], 1), v['at_M']) for k, v in res['min_own'].items()})
    print({k: (round(v['value'], 1), v['at_M']) for k, v in res['min_ref_r18'].items()})
json.dump(out, open(f'{RD}/out/wd_combined.json', 'w'), indent=1)
