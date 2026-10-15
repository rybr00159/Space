"""Empirical null + injection-recovery for the SED pair test using Gaia+2MASS+WISE (9 bands) of photometric twins."""
import numpy as np, json, sys
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/scripts')
from sed_fit import *
from astropy.table import Table
T=Table.read(D+'data/twins_sel.ecsv'); M=Table.read(D+'data/twins_merged.ecsv')
tg=M[M['source_id']==6571343840265900544][0]
cols=[('phot_bp_mean_mag',None),('phot_g_mean_mag',None),('phot_rp_mean_mag',None),('j_m','j_msigcom'),('h_m','h_msigcom'),('ks_m','ks_msigcom'),('w1mpro','w1mpro_error'),('w2mpro','w2mpro_error'),('w3mpro','w3mpro_error')]
def phot(r):
    m=np.array([float(r[c]) for c,_ in cols]); e=np.array([0.003 if ec is None else float(r[ec]) for c,ec in cols])
    e=np.hypot(e,np.array([FLOOR9[b] for b in B9])); return m,e
SEL=[('K5V','M2.5V'),('K7V','M2V'),('K7V','M3V'),('M0V','M1.5V'),('M0V','M3V'),('M0V','M4.5V'),('M0.5V','M4.5V'),('M1V','M4.5V'),('M1V','M5V'),('M1.5V','M5V'),('M2V','M4V'),('M2V','M5V'),('M2.5V','M5V'),('M3V','M5V'),('K7V',None),('M0V',None),('M2V',None)]
def analyse(m,e,d):
    s=fit(B9,m,e,d); r={'single':s['chi2'],'resid':s['resid'],'Teff':s['p'][0]}
    for pr in SEL:
        f1=fit(B9,m,e,d,pair=pr); fa=fit(B9,m,e,d,pair=pr,ampfree=True)
        r['|'.join(p for p in pr if p)]=dict(dchi2_A1=f1['chi2']-s['chi2'],A=fa['p'][3],dchi2_Afree=fa['chi2']-s['chi2'])
    return r
N=int(sys.argv[1]) if len(sys.argv)>1 else 300
idx=rng.choice(len(T),size=min(N,len(T)),replace=False)
mt,et=phot(tg); dt=1000/float(tg['parallax'])
res=dict(target=analyse(mt,et,dt),null=[],inj={})
print('target:',{k:(v if not isinstance(v,dict) else {kk:round(vv,2) for kk,vv in v.items()}) for k,v in res['target'].items() if k!='resid'},flush=True)
injpairs=[('K7V','M2V'),('M0V','M3V'),('M0V','M4.5V'),('M1V','M5V'),('M1.5V','M5V'),('M2V','M5V')]
for k,i in enumerate(idx):
    r=T[i]; m,e=phot(r); d=1000/float(r['parallax'])
    if not (np.all(np.isfinite(m)) and np.all(np.isfinite(e))): continue
    res['null'].append(analyse(m,e,d))
    # injection: add MS pair at the twin's distance (A=1), with the twin's fitted extinction
    s=fit(B9,m,e,d)
    for pr in injpairs:
        mp=pair_abs(pr,B9)+5*np.log10(d/10)+np.array([TAB[b][4] for b in B9])*s['p'][2]
        mi=S.addmag(m,mp)
        s0=fit(B9,mi,e,d); f1=fit(B9,mi,e,d,pair=pr); fa=fit(B9,mi,e,d,pair=pr,ampfree=True)
        res['inj'].setdefault('|'.join(pr),[]).append(dict(dchi2_A1=f1['chi2']-s0['chi2'],A=fa['p'][3]))
    if k%25==0:
        print(k,flush=True); json.dump(res,open(D+'results/sed_twins.json','w'))
json.dump(res,open(D+'results/sed_twins.json','w'))
print('done',len(res['null']))
