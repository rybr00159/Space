"""r03: Mira pulsation period from NEOWISE-R W1/W2 (2014-2024, ~6-month visits; own query data/neowise_mira.csv) and Gaia DR3
G epochs (ctx/lc/gaia_epoch/mira_dr3_epoch_G.csv), independent of the TESS fluxes used by the analyst.
Method: visit means (clip qual_frame>0, sigma-clip), sinusoid + per-band offset, period grid 200-400 d, chi2 with an extra
scatter term chosen so that chi2/dof = 1 at the best period (Mira cycle-to-cycle variations), Delta-chi2 of competing periods.
Then predict phase (maximum = 0) at the 2027 transit epochs from the best NEOWISE+Gaia period and compare with the analyst's
313.5 d model (last maximum BJD 2461074).
Output: out/r03_mira.json, figs/r03_mira.png"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd
B='tic214311908/analysis'
R=f'{B}/referee3_followup'
w=pd.read_csv(f'{R}/data/neowise_mira.csv'); w=w[w.qual_frame>0]
w['visit']=np.round((w.mjd-w.mjd.min())/60).astype(int)
V=[]
for v,s in w.groupby('visit'):
    for b in ('w1','w2'):
        x=s[f'{b}mpro'].dropna()
        if len(x)<3: continue
        med=np.median(x); mad=1.4826*np.median(np.abs(x-med)); x=x[np.abs(x-med)<4*max(mad,0.02)]
        V.append(dict(band=b, mjd=s.mjd.mean(), mag=x.mean(), err=max(x.std()/np.sqrt(len(x)),0.01), n=len(x)))
V=pd.DataFrame(V); print(V.pivot_table(index=np.round(V.mjd),columns='band',values='mag').round(3).to_string())
g=pd.read_csv(f'{B}/ctx/lc/gaia_epoch/mira_dr3_epoch_G.csv'); g=g[~g.rejected_by_photometry & np.isfinite(g.g_transit_mag)]
G=pd.DataFrame(dict(band='G', mjd=g.g_transit_time.values+2455197.5-2400000.5, mag=g.g_transit_mag.values, err=0.02, n=1))
D=pd.concat([V,G]); bands=sorted(D.band.unique())
def fit(P, data, extra):
    # linear least squares: per-band mean + per-band amplitude (cos, sin share phase? -> use common phase: separate amplitudes, common phase via grid on phase)
    best=None
    for ph0 in np.linspace(0,1,48,endpoint=False):
        cols=[]; 
        for b in bands:
            m=(data.band==b).values.astype(float); cols.append(m); 
            cols.append(m*np.cos(2*np.pi*((data.mjd.values/P)-ph0)))
        A=np.array(cols).T; e=np.sqrt(data.err.values**2+data.band.map(extra).values**2)
        c,*_=np.linalg.lstsq(A/e[:,None], data.mag.values/e, rcond=None)
        chi=np.sum(((data.mag.values-A@c)/e)**2)
        if best is None or chi<best[0]: best=(chi,ph0,c)
    return best
Ps=np.arange(200,400.01,0.5)
res={}
for lab,data,extra in (('NEOWISE',D[D.band!='G'],{'w1':0.12,'w2':0.10}),('NEOWISE+Gaia',D,{'w1':0.12,'w2':0.10,'G':0.15})):
    chi=np.array([fit(P,data,extra)[0] for P in Ps])
    i=int(np.argmin(chi)); dof=len(data)-(2*data.band.nunique()+2)
    # rescale to chi2/dof=1
    s=chi[i]/dof; dchi=(chi-chi[i])/s
    mins=[(float(Ps[j]),round(float(dchi[j]),1)) for j in range(1,len(Ps)-1) if chi[j]<chi[j-1] and chi[j]<chi[j+1]]
    mins=sorted(mins,key=lambda t:t[1])[:6]
    w_=np.exp(-0.5*dchi); w_/=w_.sum(); Pm=float(np.sum(w_*Ps)); Psd=float(np.sqrt(np.sum(w_*(Ps-Pm)**2)))
    at={P:round(float(np.interp(P,Ps,dchi)),1) for P in (294,313.5,326)}
    best=fit(Ps[i],data,extra)
    res[lab]=dict(P_best=float(Ps[i]), P_mean=Pm, P_sd=Psd, chi2_min=float(chi[i]), dof=int(dof), scale=float(s), local_minima=mins, dchi2_at=at,
                  phase0=float(best[1]))
    print(lab, res[lab])
# predicted phase at 2027 events from NEOWISE+Gaia best solution vs analyst
an=json.load(open(f'{B}/followup/out/mira.json'))['period_fit']
ev={'E126 2027-06-01':2461558.04,'E127 2027-06-07':2461563.73,'E137 2027-08-03':2461620.65,'E131 2027-06-29':2461586.50}
r=res['NEOWISE+Gaia']; P=r['P_best']
# epoch of maximum brightness: model mag = mean + amp*cos(2pi(t/P - ph0)); maximum light where cos term is most negative if amp>0
data=D; best=fit(P,data,{'w1':0.12,'w2':0.10,'G':0.15}); c=best[2]
ampG=c[bands.index('G')*2+1]
ph_max = best[1] + (0.5 if ampG>0 else 0.0)
out={}
for k,t in ev.items():
    mjd=t-2400000.5
    ph=((mjd/P)-ph_max)%1
    ph_an=((t-an['last_max_before_2026_5_BJD'])/an['P_best_d'])%1
    out[k]=dict(phase_neowise_gaia=round(float(ph),2), phase_analyst=round(float(ph_an),2))
print(out)
res['phase_2027']=out
json.dump(res, open(f'{R}/out/r03_mira.json','w'), indent=1)
