"""s15: final adopted planet parameters and a compact summary of all shape-analysis numbers.
Planet parameters from the s03 'single' (rho prior 9.41 +- 0.88, circular) and 'twin' (rho 15.5 +- 1.45, dilution 2)
posteriors of both photometries. Adopted value = mean of the INV and REF medians; uncertainty = larger of the two
68% half-widths added in quadrature to half the INV-REF difference (the two use the same pixels, so they are not
averaged as independent). Derived: Rp = k R* (R* ~ N(0.381, 0.012) single; N(0.289, 0.009) twin), a = aR R*,
Teq = Teff sqrt(R*/2a) (A_B = 0, full redistribution), S = L/a^2.
Output: out/final_params.json"""
import json, numpy as np
SD = 'tic214311908/analysis/shape'
rng = np.random.default_rng(3)
RSUN_REARTH = 109.076; RSUN_AU = 0.0046505
host = dict(single=dict(R=(0.3813, 0.0115), Teff=(3350, 54), L=(0.016488, 0.00034)),
            twin=dict(R=(0.2887, 0.0089), Teff=(3237, 52), L=(0.016488/2, 0.00017)))
def derived(ds, model):
    z = np.load(f'{SD}/out/fit_{ds}_{model}.npz'); n = len(z['k'])
    h = host['single' if model.startswith('single') else 'twin']
    R = rng.normal(*h['R'], n); Te = rng.normal(*h['Teff'], n); L = rng.normal(*h['L'], n)
    a_au = z['aR']*R*RSUN_AU
    d = dict(k=z['k'], b=z['b'], aR=z['aR'], inc=z['inc'], T14_h=z['T14']*24, T23_h=z['T23']*24, flat=z['T23']/z['T14'], rho=z['rho'],
             Rp_Re=z['k']*R*RSUN_REARTH, a_au=a_au, Teq_K=Te*np.sqrt(R*RSUN_AU/(2*a_au)), S_Searth=L/a_au**2, e=z['e'])
    return {kk: np.percentile(v, [15.865, 50, 84.135]) for kk, v in d.items()}
out = {}
for model in ['single', 'twin', 'single_ecc', 'twin_ecc']:
    try:
        A = derived('inv', model); B = derived('ref', model)
    except FileNotFoundError:
        continue
    res = {}
    for kk in A:
        ma, mb = A[kk][1], B[kk][1]
        ea = 0.5*(A[kk][2] - A[kk][0]); eb = 0.5*(B[kk][2] - B[kk][0])
        res[kk] = dict(adopted=float(0.5*(ma+mb)), err=float(np.hypot(max(ea, eb), 0.5*abs(ma-mb))), inv=[float(x) for x in A[kk]], ref=[float(x) for x in B[kk]])
    out[model] = res
    print(model, {kk: (round(v['adopted'], 4), round(v['err'], 4)) for kk, v in res.items()})
json.dump(out, open(f'{SD}/out/final_params.json', 'w'), indent=1)
