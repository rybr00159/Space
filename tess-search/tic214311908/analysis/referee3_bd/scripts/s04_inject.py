"""s04: fake-epoch injection calibration of the grazing-BD shape statistic Lambda_G = chi2(grazBD, k >= 0.19) - chi2(planet).
Truths (from s03, beta noise, per data set): 'planet' = best planet fit of the real cleaned data; 'grazBD' = best
grazing-BD fit of the real data (k 0.19, b ~0.96, rho_c ~1.8). Each trial injects the truth (D = 1, nominal LD, dT0 = 0)
into a set of fake windows with the real sector composition and refits both scenarios with lc3 (grid + Nelder-Mead,
dT0 and LD free), giving Lambda_G under each truth.
Fake windows: REF = referee_phot fake windows (data/win_s*.npz, fake_*; flux 1 + aM/979, per-window error = robust OOT
std/979 x clean's REF beta_transit per sector), drawn without replacement within a trial and cycling over shuffled
pools across trials, composition S12:3, S39:4, S66:3, S93:4, S103:4, S104:3 (as the real REF set). INV = the cleaned
INV curve at fake ephemerides T0 + phi P, phi on a 0.02 grid in 0.13-0.39 and 0.61-0.87 (in-transit parts
independent; baselines partly shared), windows avoiding |dt_real| < 0.3 d; errors flux_err x beta_transit.
Usage: s04_inject.py <ref|inv> <planet|grazBD> <ntrials> <seed>  -> out/s04_inject_<ds>_<truth>.csv"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from lc3 import *
import time, json
ds, truth, ntr, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
fr = json.load(open(f'{R3}/out/s03_fit_real.json'))['fits'][f'{ds}_beta'][truth]
TK, TB, TR = fr['k'], fr['b'], fr['rho_c']
truth_curves = binned(TK, TB, float(aR_of_rho(TR)))
print('truth', ds, truth, TK, TB, TR, flush=True)
bank = TemplateBank(np.r_[[0.10, 0.11, 0.12, 0.13, 0.14], [0.19, 0.21, 0.24, 0.28]], nb=26, rhos=np.exp(np.linspace(np.log(0.7), np.log(40.), 26)))
print('bank', len(bank.items), flush=True)
# ---------------- fake data sets
sets = []
if ds == 'ref':
    import pandas as pd
    beta = pd.read_csv(f'{B}/clean/lc_M3V_REF_clean.csv').groupby('sector').beta_transit.first().to_dict()
    comp = {12: 3, 39: 4, 66: 3, 93: 4, 103: 4, 104: 3}
    pools = {}
    for s in comp:
        z = np.load(f'{B}/referee_phot/data/win_s{s}.npz', allow_pickle=True)
        t, aM, w, oot, tc = z['fake_t'], z['fake_aM'], z['fake_win'], z['fake_oot'], z['fake_tc']
        cad = float(z['exptime'])*DAY
        P_ = []
        for i in range(len(tc)):
            m = w == i
            if m.sum() < 20: continue
            sg = rstd(aM[m & oot])/979.
            P_.append(dict(t=t[m], f=1 + aM[m]/979., e=np.full(m.sum(), sg*beta[s]), cad=np.full(m.sum(), round(cad)), tc=np.full(m.sum(), tc[i]), s=s))
        order = list(rng.permutation(len(P_)))
        pools[s] = (P_, order)
    ptr = {s: 0 for s in comp}
    for it in range(ntr):
        parts = []
        for s, n in comp.items():
            P_, order = pools[s]
            for j in range(n):
                if ptr[s] >= len(order): order[:] = list(rng.permutation(len(P_))); ptr[s] = 0
                parts.append((P_[order[ptr[s]]], s*1000 + j)); ptr[s] += 1
        cat = lambda key: np.concatenate([p[0][key] for p in parts])
        win = np.concatenate([np.full(len(p[0]['t']), p[1]) for p in parts])
        sets.append(LC(cat('t'), cat('f'), cat('e'), np.concatenate([np.full(len(p[0]['t']), p[0]['s']) for p in parts]), cat('cad'), win, cat('tc')))
else:
    d = inv_frame()
    dtr = d.time_btjd_tdb.values - (T0C + np.round((d.time_btjd_tdb.values - T0C)/PC)*PC)
    phis = np.r_[np.arange(0.13, 0.391, 0.02), np.arange(0.61, 0.871, 0.02)]
    phis = phis[rng.permutation(len(phis))][:ntr]
    for ph in phis:
        sets.append(load_inv(T0=T0C + ph*PC, avoid=np.abs(dtr) > 0.3))
rows = []
t0 = time.time()
for it, lc in enumerate(sets):
    s_true = lc.template(truth_curves)
    lc.f = lc.f + s_true; lc.y = lc.f - 1; lc.yp = lc.perp(lc.y)
    rp = fit_scen(lc, bank, 'planet', nstart=2)
    rg = fit_scen(lc, bank, 'grazBD', nstart=2)
    row = dict(trial=it, N=len(lc.t), nw=int(lc.nw), Lambda_G=rg['chi2'] - rp['chi2'],
               **{f'p_{k}': rp[k] for k in ('chi2', 'k', 'b', 'rho_c', 'T14_h', 'flat', 'dT0_min')},
               **{f'g_{k}': rg[k] for k in ('chi2', 'k', 'b', 'rho_c', 'T14_h', 'dT0_min')})
    rows.append(row)
    print(it, round(row['Lambda_G'], 2), 'planet', round(rp['k'], 3), round(rp['b'], 2), round(rp['rho_c'], 1), ' graz', round(rg['k'], 3), round(rg['b'], 3), round(rg['rho_c'], 2), round(time.time() - t0), flush=True)
    import pandas as pd
    pd.DataFrame(rows).to_csv(f'{R3}/out/s04_inject_{ds}_{truth}.csv', index=False)
