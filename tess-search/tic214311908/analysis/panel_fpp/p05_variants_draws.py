"""p05: sign-fixed re-run of the referee2_fpp prior-dimension variants and joint prior draws (which referee3_bd did not
re-run: its 64/192 variants vary only light curve x noise x planet e prior x BD e law x f_BD).
Adopted settings: BD e law Ferreira+26 short-period, f_BD 1.3e-4 (referee3_bd). Outputs: out/p05_variants.csv,
out/p05_draws.csv, out/p05_summary.json"""
import json, time, numpy as np, pandas as pd
from fpplib import *
os.makedirs('out', exist_ok=True)
t0 = time.time()
V = []
def var(name, ds, **kw):
    rvf = kw.pop('rvf', 1.0)
    p = dict(ADOPT, ds=ds); p.update(kw)
    r = run(p, rvf=rvf); r.update(variant=name, ds=ds); V.append(r)
    print(f"{ds:9s} {name:45s} FPP {r['FPP']:.5f} strict {r['FPP_strict']:.4f} NFPP {r['NFPP']:.2e} PTP {r['PTP']:.2e} BD {r['BD']:.1e}  t={time.time()-t0:.0f}s", flush=True)
LR = R.LR
pess = dict(beta_N=-1.0, f_mult=0.40, F_twin=0.0, S_close=0.5, ipd='L_ipd84', ast='soft', ecc='vaneylen+bdsp', f_BD=3.5e-4,
            lr_pairA=10**-2.8, dchi2_wd=6.0, r_wd=2e-4)
for ds in ('ref_beta', 'inv_beta'):
    var('adopted (kipping, bdsp, f_BD 1.3e-4)', ds)
    var('adopted + Gaia RV factor 0.49', ds, rvf=0.49)
    var('planet e Rayleigh 0.19 (Sagear & Ballard)', ds, ecc='rayleigh19+bdsp')
    for v in (0.001, 0.01): var(f'f_N={v}', ds, f_N=v)
    for v in (-1.0, 2.0): var(f'beta_N={v}', ds, beta_N=v)
    var('radius prior R^-3', ds, alpha=3.0)
    for v in (0.2, 0.45): var(f'f_mult={v}', ds, f_mult=v)
    for v in (0.05, 1.0): var(f'S_close={v}', ds, S_close=v)
    for v in (0.0, 0.25): var(f'F_twin={v}', ds, F_twin=v)
    var('CMD x only', ds, cmd='xonly'); var('CMD flat (ignored)', ds, cmd='flat')
    var('CMD analyst Gaussian 0.30', ds, cmd='gauss', sig_s=0.30)
    for ip in ('L_ipd84', 'L_ipd_red', 'L_ipd_far', 'off'): var(f'IPD {ip}', ds, ipd=ip)
    var('astrometry off', ds, ast='off')
    var('all archival constraints off', ds, ipd='off', ast='off', dec='off', wide=False)
    for v in (0.5, 1.0): var(f'dT_twin={v}', ds, dT_twin=v)
    var('pair A LR sys-inflated worst (10^-2.8)', ds, lr_pairA=10**LR['sysinfl']['log10LR_sysinfl_at_min_variant'])
    var('pair A LR t median (typical)', ds, lr_pairA=10**LR['log10LR_t']['median'])
    var('WD dchi2 6, r_wd 2e-4', ds, dchi2_wd=6.0, r_wd=2e-4)
    var('WD dchi2 22 (RV-smearing term added)', ds, dchi2_wd=22.0)
    var('PESSIMISTIC realistic (sign fixed, f_BD 3.5e-4)', ds, **pess)
    var('PESSIMISTIC + IPD only beyond 0.25"', ds, **dict(pess, ipd='L_ipd_far'))
dfv = pd.DataFrame(V); dfv.to_csv('out/p05_variants.csv', index=False)
# joint draws (r05 ranges, sign fixed, BD inputs updated to referee3_bd)
rng2 = np.random.default_rng(2026)
def lu(lo, hi): return float(np.exp(rng2.uniform(np.log(lo), np.log(hi))))
Rr = []
NDR = int(os.environ.get('PF_NDRAW', 1500))
for i in range(NDR):
    pe = str(rng2.choice(['kipping', 'vaneylen', 'circ', 'rayleigh19'])); be = 'bdsp' if rng2.random() < 0.8 else 'bdlp'
    p = dict(ds=str(rng2.choice(['ref_beta', 'inv_beta', 'ref_white', 'inv_white'], p=[0.4, 0.4, 0.1, 0.1])),
             ecc=f'{pe}+{be}', alpha=float(rng2.choice([0.0, 2.0, 3.0])),
             f_N=lu(0.001, 0.01), f_S=lu(0.0001, 0.003), f_J=lu(0.0005, 0.005), f_BD=lu(5e-5, 3.5e-4),
             beta_N=float(rng2.uniform(-1, 2)), beta_S=float(rng2.uniform(0, 2)), beta_J=float(rng2.uniform(0, 3)),
             f_mult=float(rng2.uniform(0.22, 0.42)), mu_loga=float(rng2.uniform(1.0, 1.5)), sig_loga=float(rng2.uniform(0.9, 1.4)),
             gamma_q=float(rng2.uniform(-0.5, 1.0)), F_twin=float(rng2.uniform(0, 0.3)), S_close=lu(0.05, 1.0),
             cmd=str(rng2.choice(['colour', 'colour', 'xonly'])), cmd_Ft=float(rng2.choice([0.0, 0.1, 0.25])),
             cmd_sb=['marg', 'marg', 0.30, 0.35, 0.40, 0.45][int(rng2.integers(6))], xobs=float(rng2.uniform(-0.85, -0.78)),
             ipd=str(rng2.choice(['L_ipd', 'L_ipd84', 'L_ipd_red', 'L_ipd_far'])), ast=str(rng2.choice(['base', 'soft'])),
             dT_twin=float(rng2.uniform(0.5, 1.0)), r_wd=lu(2e-6, 2e-4), dchi2_wd=float(rng2.uniform(6, 13)),
             f_occ=lu(0.03, 0.3), dchi2_heb=float(rng2.uniform(20, 33.7)),
             lr_pairA=10**float(rng2.uniform(-9.0, -2.8)), S_BD=bool(rng2.random() < 0.5))
    rvf = float(rng2.uniform(0.25, 0.80))
    r = run(p, rvf=rvf, seed=1000 + i)
    r.update({k: (v if not isinstance(v, (list, tuple)) else str(v)) for k, v in p.items()}); r['rvf'] = rvf
    Rr.append(r)
    if i % 100 == 0: print('draw', i, f"t={time.time()-t0:.0f}s", flush=True)
dfr = pd.DataFrame(Rr); dfr.to_csv('out/p05_draws.csv', index=False)
q = lambda x, p_: float(np.percentile(x, p_))
summ = dict(n=len(dfr), FPP_pct_5_16_50_84_95=[q(dfr.FPP, p_) for p_ in (5, 16, 50, 84, 95)],
            frac_FPP_lt_0p015=float(np.mean(dfr.FPP < 0.015)), frac_FPP_lt_0p01=float(np.mean(dfr.FPP < 0.01)),
            NFPP_pct_5_50_95=[q(dfr.NFPP, p_) for p_ in (5, 50, 95)], frac_NFPP_lt_1em3=float(np.mean(dfr.NFPP < 1e-3)),
            FPP_strict_pct_5_50_95=[q(dfr.FPP_strict, p_) for p_ in (5, 50, 95)], frac_strict_lt_0p015=float(np.mean(dfr.FPP_strict < 0.015)),
            FPP_system_pct_5_50_95=[q(dfr.FPP_system, p_) for p_ in (5, 50, 95)],
            FPP_median_by_ds=dfr.groupby('ds').FPP.median().to_dict(), FPP_median_by_ipd=dfr.groupby('ipd').FPP.median().to_dict(),
            variants_FPP_range=[float(dfv.FPP.min()), float(dfv.FPP.max())], variants_n_ge_0p015=int((dfv.FPP >= 0.015).sum()),
            method=__doc__)
json.dump(summ, open('out/p05_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
