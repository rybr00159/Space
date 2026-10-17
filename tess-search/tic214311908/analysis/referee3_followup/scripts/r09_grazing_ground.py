"""r09: can one ground multi-band transit close the grazing (eccentric) brown-dwarf channel on the M3V by shape alone?
Truth: planet model (k 0.119, b 0.52, a/R* 25.4, P 5.6915 d, i/z quadratic LD u1 0.30 u2 0.30 assumed, r u1 0.45 u2 0.25),
2-min cadence, MuSCAT4-like white noise per 2 min (r 1.08, i 0.74, z 0.75 ppt from out/r07_snr.csv), 1.5 h baseline each side.
Alternative: grazing dark body with k >= 0.19 (BD radius prior lower bound, ../bd b03: k 0.19-0.30), free b, a/R* (any rho_c),
free depth normalisation not allowed (D = 1, undiluted, the M3V channel), linear baseline per band free.
Statistic: Delta chi2 = min chi2(grazing) - chi2(truth) for noiseless truth (Asimov), white noise, then deflated by a
correlated-noise factor beta^2 = 4 and 9 (beta = 2, 3) as a pessimistic correction. Compare with TESS: Lambda_G = 13.6 (REF) /
6.8 (INV) with beta noise (../bd).
Output: out/r09_grazing_ground.json"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, batman
from scipy.optimize import minimize
R='tic214311908/analysis/referee3_followup'
P=5.6915; t=np.arange(-(0.85+1.5)/24,(0.85+1.5)/24,2/1440)
LD={'r':[0.45,0.25],'i':[0.30,0.30],'z':[0.25,0.30]}; SIG={'r':1.08e-3,'i':0.74e-3,'z':0.75e-3}
def lc(k,b,aR,u):
    p=batman.TransitParams(); p.t0=0; p.per=P; p.rp=k; p.a=aR; p.inc=np.degrees(np.arccos(min(b/aR,1))); p.ecc=0; p.w=90; p.u=u; p.limb_dark='quadratic'
    return batman.TransitModel(p,t,supersample_factor=5,exp_time=2/1440).light_curve(p)
truth={bd:lc(0.119,0.52,25.4,LD[bd]) for bd in 'riz'}
print('truth depth', {b:round(1e3*(1-truth[b].min()),2) for b in truth}, 'T14 h', round(float(np.ptp(t[truth['i']<1])*24),2))
def chi2(x):
    k,b,aR=x
    if k<0.19 or k>0.4 or aR<3 or aR>80 or b<0 or b>1+k: return 1e12
    c=0
    for bd in 'riz':
        m=lc(k,b,aR,LD[bd]); A=np.vstack([m, m*t]).T   # baseline: multiplicative linear trend
        coef,*_=np.linalg.lstsq(A/SIG[bd],truth[bd]/SIG[bd],rcond=None)
        c+=np.sum(((truth[bd]-A@coef)/SIG[bd])**2)
    return c
best=None
for k0 in (0.19,0.22,0.26,0.30):
    for b0 in (0.95,1.0,1.05,1.1):
        for a0 in (8,12,16,20):
            r=minimize(chi2,[k0,b0,a0],method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-3,maxiter=3000))
            if best is None or r.fun<best.fun: best=r
print('best grazing', best.x, 'dchi2 white', best.fun)
out=dict(best_grazing_k_b_aR=[float(v) for v in best.x], dchi2_white_one_event_riz=float(best.fun),
         dchi2_beta2=float(best.fun/4), dchi2_beta3=float(best.fun/9), TESS_LambdaG_REF_INV=[13.6,6.8])
json.dump(out,open(f'{R}/out/r09_grazing_ground.json','w'),indent=1); print(out)
