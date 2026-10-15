"""Per-event eclipse times for every primary and secondary eclipse, with empirical errors.
Template: global best-fit physical model (results/best_<tag>.json), with the product-specific depth scale.
Empirical errors: the template (primary or secondary shape) is injected into the real light curve at fake
conjunctions at orbital phases 0.25 and 0.75 (eclipse-free) with a random offset in +-5 min, and recovered with the
identical procedure; per product and eclipse type the robust sd of (recovered - injected) is the timing error, and
its ratio to the formal error is the error-scale factor.
Usage: e02_etv_measure.py <best-fit tag> [timewarp_c_min_per_yr2]   (the optional second argument applies the
positive-control time warp t -> t + c (t - 2022.5 yr)^2 to all time stamps before measuring)."""
import numpy as np, sys, json, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import data as D, etv
V = D.V
tag = sys.argv[1] if len(sys.argv) > 1 else 'default'
warp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
best = json.load(open(f'{V}/results/best_{tag}.json'))
p = best['p']
scale = {k[2:]: v for k, v in p.items() if k.startswith('s_')}
tm = etv.Template(p, scale)
T0, P = p['T0'], p['P']
TMID = 3165.6697 - 365.25 * 0.6   # BTJD of ~2022.5, warp pivot
rng = np.random.default_rng(1 if warp == 0 else 2)
lcs = [D.load_spoc(1), D.load_spoc(28), D.load_spoc(68), D.load_tc(95), D.load_tc(105), D.load_tc(106)]
if len(sys.argv) > 3 and sys.argv[3] == 'qlp95':
    lcs[3] = D.load_qlp(95)
summ = json.load(open(f'{V}/results/data_summary.json'))['default']['info']
events = []; inj = []
for lc in lcs:
    lc = D.clip(lc)
    key = f"{lc['prod']}{lc['sector']}"
    if warp != 0.0:
        lc['t'] = lc['t'] + warp / 1440.0 * ((lc['t'] - TMID) / 365.25) ** 2
    sig = summ[key]['sigma'] if key in summ else lc['sig_pt']
    for w in D.windows(lc, P=P, T0=T0, phases=(0.0, 0.5)):
        i = w['idx']; ph = w['ph']
        tc = T0 + w['n'] * P + (tm.dsec if ph == 0.5 else 0.0)
        tt = lc['t'][i]
        if np.abs(tt - tc).max() > 0.215: i = i[np.abs(lc['t'][i] - tc) < 0.2]; tt = lc['t'][i]
        dl, sf, chi, n = etv.measure(tt, lc['f'][i], np.full(i.size, sig), np.full(i.size, lc['texp']), tc, ph, key, tm)
        events.append(dict(key=key, sector=lc['sector'], n=w['n'], ph=ph, tc_pred=tc, delta_min=dl * 1440, sig_formal_min=sf * 1440, chi2=chi, npts=n,
                           t_obs=tc + dl if np.isfinite(dl) else np.nan))
    if warp != 0.0:
        continue
    for w in D.windows(lc, P=P, T0=T0, phases=(0.25, 0.75)):
        i = w['idx']; tcf = w['tc']
        for ph in (0.0, 0.5):
            for rep in range(2):
                di = rng.uniform(-5, 5) / 1440
                y = lc['f'][i] * (1 - tm.depth(lc['t'][i] - tcf - di, lc['texp'], ph, key))
                dl, sf, chi, n = etv.measure(lc['t'][i], y, np.full(i.size, sig), np.full(i.size, lc['texp']), tcf, ph, key, tm)
                inj.append(dict(key=key, ph=ph, n=w['n'], fake_phase=w['ph'], inj_min=di * 1440, rec_min=dl * 1440, sig_formal_min=sf * 1440))
    print(key, 'done', flush=True)
out = dict(tag=tag, warp_min_per_yr2=warp, T0=T0, P=P, dsec_d=tm.dsec, events=events, injections=inj)
stats = {}
if warp == 0.0:
    for key in sorted(set(e['key'] for e in events)):
        for ph in (0.0, 0.5):
            a = np.array([(x['rec_min'] - x['inj_min'], x['sig_formal_min']) for x in inj if x['key'] == key and x['ph'] == ph and np.isfinite(x['rec_min'])])
            ev = np.array([(e['delta_min'], e['sig_formal_min']) for e in events if e['key'] == key and e['ph'] == ph and np.isfinite(e['delta_min'])])
            err = a[:, 0]
            rsd = 1.4826 * np.median(np.abs(err - np.median(err)))
            z = err / a[:, 1]
            st = dict(n_inj=int(len(a)), bias_min=float(np.median(err)), robust_sd_min=float(rsd), sd_min=float(np.std(err)),
                      frac_out_3rsd=float(np.mean(np.abs(err - np.median(err)) > 3 * rsd)),
                      zscale=float(1.4826 * np.median(np.abs(z - np.median(z)))), median_formal_inj=float(np.median(a[:, 1])),
                      n_ev=int(len(ev)), median_formal_ev=float(np.median(ev[:, 1])),
                      ev_robust_rms=float(1.4826 * np.median(np.abs(ev[:, 0] - np.median(ev[:, 0])))), ev_median=float(np.median(ev[:, 0])))
            stats[f'{key}_ph{ph}'] = st
            print(f"{key} ph{ph}: inj N={st['n_inj']} bias {st['bias_min']:+.2f} rsd {rsd:.2f} min (formal {st['median_formal_inj']:.2f}, z-scale {st['zscale']:.2f}, out>3rsd {st['frac_out_3rsd']:.3f}) | "
                  f"real N={st['n_ev']} median O-C {st['ev_median']:+.2f} robust rms {st['ev_robust_rms']:.2f} min (formal {st['median_formal_ev']:.2f})", flush=True)
out['stats'] = stats
sfx = f'_warp{warp:g}' if warp != 0 else ''
if len(sys.argv) > 3: sfx += '_' + sys.argv[3]
json.dump(out, open(f'{V}/results/etv_{tag}{sfx}.json', 'w'), indent=1)
