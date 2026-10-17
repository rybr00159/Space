"""r11: how robust is the brown-dwarf (BD) channel added in r05? Re-evaluates the r05 scenario model for the BD terms
under (i) the four light-curve/noise choices (REF/INV x beta-inflated/white), (ii) the three eccentricity priors, and
(iii) a calibrated grazing penalty: referee2_shape r07b found that a V-shaped (grazing) truth reproduces the observed
flat-bottom preference in 0/25 REF and 1/25 INV injections, i.e. p(data | grazing) / p(data | flat) <~ 0.04; the grid
likelihood with beta noise gives exp(-dchi2/2) = 0.05 (REF) / 0.11 (INV) for the best grazing solution, so the beta grid
is at least as lenient as the calibrated value.
Output: out/r11_bd_noise.json"""
import numpy as np, json, sys
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
import r05_fpp as R
out = {}
for ds in ('ref_beta', 'inv_beta', 'ref_white', 'inv_white'):
    for ecc in ('kipping', 'vaneylen', 'circ'):
        p = dict(R.BASE); p.update(ds=ds, ecc=ecc)
        r = R.evaluate(p, detail=True)
        bd = r['post']['BD_on_M3V_single'] + r['post']['BD_on_primary_of_binary'] + r['post']['BDHEB_on_companion']
        out[f'{ds}|{ecc}'] = dict(FPP=r['FPP'], FPP_without_BD_terms=None, BD_total=bd, BD_on_M3V=r['post']['BD_on_M3V_single'],
                                  BD_primary=r['post']['BD_on_primary_of_binary'], BDHEB=r['post']['BDHEB_on_companion'],
                                  PTP=r['post']['PTP_companion'], dominant=r['dominant_FP'])
        p0 = dict(p); p0['f_BD'] = 0.0
        out[f'{ds}|{ecc}']['FPP_without_BD_terms'] = R.evaluate(p0)['FPP']
        print(f"{ds:9s} {ecc:8s} FPP {r['FPP']:.4f} (no BD {out[f'{ds}|{ecc}']['FPP_without_BD_terms']:.4f})  BD total {bd:.2e} "
              f"[M3V {r['post']['BD_on_M3V_single']:.1e}, primary {r['post']['BD_on_primary_of_binary']:.1e}, companion {r['post']['BDHEB_on_companion']:.1e}]  PTP {r['post']['PTP_companion']:.1e}", flush=True)
# best grazing solution on the M3V (D = 1) in each grid: dchi2 relative to the global planet optimum
gr = {}
for ds in ('ref_beta', 'inv_beta', 'ref_white', 'inv_white'):
    g = np.load(f'{R.R2}/out/r03_grid_{ds}.npz')
    c = g['chi2'] - g['chi2_ref']; b = g['b']; k = g['k']
    m = (b[None, :] >= 0.95) & (k[0] >= 0.19)       # D = 1, BD-sized k on the M3V (R >= 8.2 R_e)
    cc = np.where(m, c[0], np.nan)
    j, l = np.unravel_index(np.nanargmin(cc), cc.shape)
    gr[ds] = dict(dchi2_best_grazing_BD=float(cc[j, l]), rho_c=float(g['rho'][j]), b=float(b[l]), k=float(k[0, j, l]),
                  LR_exp_minus_half_dchi2=float(np.exp(-0.5*cc[j, l])))
    print(ds, gr[ds], flush=True)
json.dump(dict(method=__doc__, scenarios=out, grazing_on_M3V=gr), open(f'{R.R2}/out/r11_bd_noise.json', 'w'), indent=1, default=float)
