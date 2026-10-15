"""Decompose TRICERATOPS 1.1.0 SEB vs BEB evidence into prior ratio and likelihood ratio (no contrast curve: both use 2.2 arcsec)."""
import numpy as np, pandas as pd, json
from triceratops.priors import lnprior_bound_EB, lnprior_background
from triceratops.funcs import trilegal_results
S='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
fn=S+'tri/71316387_gaia_background.csv'
Tm=trilegal_results(fn,11.358)[0]; N=len(Tm)
lb=float(lnprior_bound_EB(1.05,3.3622,np.array([2.2]))[0])
lbg=float(lnprior_background(N,np.array([2.2]))[0])
out=dict(N_comp=N,Tmag_range=[float(np.nanmin(Tm)),float(np.nanmax(Tm))],prior_bound_2p2=np.exp(lb),prior_bg_2p2=np.exp(lbg),ln_ratio=lb-lbg,ratio=float(np.exp(lb-lbg)))
for s in [0.1,0.5,1.0,1.24]:
    out[f'prior_bound_{s}']=float(np.exp(lnprior_bound_EB(1.05,3.3622,np.array([s]))[0]))
    out[f'prior_bg_{s}']=float(np.exp(lnprior_background(N,np.array([s]))[0]))
# fraction of background population that could supply 5900 ppm deep eclipse (undiluted depth <=0.5,1.0)
f=10**((11.358-Tm)/2.5); fr=f/(1+f)
for dmax in (0.5,1.0): out[f'frac_bg_can_host_d{dmax}']=float(np.mean(0.0059/fr<=dmax))
print(json.dumps(out,indent=1))
json.dump(out,open(S+'results/tri_prior_decomp.json','w'),indent=1)
