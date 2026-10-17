"""r06: independent mass, K, RV epochs, TSM/ESM and ranking.
* Mass from R: (a) Chen & Kipping 2017 Neptunian segment inverted (R = 1.230 (M/2.04)^0.589, 14.6% fractional radius
  scatter), (b) Otegi+20 volatile-rich R = 0.70 M^0.63 (scatter ~ 7% in R assumed), (c) empirical: NEA planets 4-7 R_E with
  mass measured (>3 sigma), host Teff < 4000 K (own TAP query), log-log fit with intrinsic scatter.
* K = (2 pi G / P)^(1/3) Mp sin i / (M*+Mp)^(2/3) / sqrt(1-e^2), i = 88.8 deg, Rp ~ N(4.90, 0.28), M* ~ N(0.370, 0.011).
* Epochs for a 5-sigma K with random phases: N = 2 (5 sigma_eff / K)^2 (sigma_K = sigma_eff sqrt(2/N)); evaluated over the
  whole K distribution (median epochs and 84th percentile), not only at the median K.
* TSM (Kempton+18 eq. 1, scale 1.15 for 4-10 R_E), ESM (eq. 4, 7.5 um, Tday = 1.10 Teq), Teq = Teff sqrt(R*/a)(1/4)^(1/4).
* Ranking vs NEA planets 4-7 R_E, P < 10 d, Teff < 4000 K with J and K mags (own query).
Outputs: out/r06_rv_tsm.json, data/nea_mdw_neptunes.csv"""
import os, json, io
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd, requests
B='tic214311908/analysis/referee3_followup'
rng=np.random.default_rng(5); N=200000
Rp=rng.normal(4.90,0.28,N); Ms=rng.normal(0.370,0.011,N); Rs=rng.normal(0.383,0.012,N)
Teff=rng.normal(3350,100,N); aR=rng.normal(25.4,0.8,N); inc=np.radians(88.8)
P=5.6915020; G=6.674e-11; MSUN=1.989e30; ME=5.972e24
def Kms(Mp_E, Ms_):
    return (2*np.pi*G/(P*86400))**(1/3)*Mp_E*ME*np.sin(inc)/((Ms_*MSUN+Mp_E*ME)**(2/3))
# (a) CK17: Neptunian segment; R(M) median = 1.008*2.04^0.279 * (M/2.04)^0.589, fractional scatter 14.6% in R
R_T1=1.008*2.04**0.279
Rtrue=Rp/np.clip(rng.normal(1,0.146,N),0.5,None)   # radius the median relation would need
M_ck=2.04*(Rtrue/R_T1)**(1/0.589)
# (b) Otegi volatile-rich R = 0.70 M^0.63 with 7% scatter (assumed)
M_ot=(Rp/np.clip(rng.normal(1,0.07,N),0.5,None)/0.70)**(1/0.63)
# (c) empirical NEA
url='https://exoplanetarchive.ipac.caltech.edu/TAP/sync'
q=("select pl_name,pl_rade,pl_radeerr1,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_orbper,st_teff,st_rad,st_mass,sy_jmag,sy_kmag,pl_eqt,pl_bmassprov,tran_flag,pl_bmasselim "
   "from pscomppars where pl_rade between 3.5 and 8 and st_teff < 4000 and st_teff > 2500 and pl_orbper < 30 and tran_flag = 1")
r=requests.get(url, params=dict(query=q, format='csv'), timeout=120); T=pd.read_csv(io.StringIO(r.text)); T.to_csv(f'{B}/data/nea_mdw_neptunes.csv', index=False)
print('NEA rows', len(T))
E=T[(T.pl_rade>=4)&(T.pl_rade<=7)&T.pl_bmasse.notna()&(T.pl_bmassprov=='Mass')&(T.pl_bmasselim.fillna(0)==0)]
E=E[E.pl_bmasse/(0.5*(E.pl_bmasseerr1.abs()+E.pl_bmasseerr2.abs()))>3]
print(E[['pl_name','pl_rade','pl_bmasse','st_teff','pl_orbper']].to_string())
x=np.log(E.pl_rade.values); y=np.log(E.pl_bmasse.values)
A=np.vstack([np.ones_like(x),x-np.log(4.9)]).T; coef,res,_,_=np.linalg.lstsq(A,y,rcond=None)
sc=np.sqrt(np.sum((y-A@coef)**2)/(len(y)-2))
M_emp=np.exp(coef[0]+coef[1]*(np.log(Rp)-np.log(4.9))+rng.normal(0,sc,N))
out=dict(n_empirical=int(len(E)), empirical_lnM_at4p9=float(coef[0]), slope=float(coef[1]), scatter_ln=float(sc), planets=list(E.pl_name))
pc=lambda v,q=(2.5,16,50,84,97.5): [round(float(t),2) for t in np.percentile(v,q)]
SIG={'ESPRESSO':3.4,'NIRPS+HARPS':6.4,'HARPS':8.8,'HARPS (TOI-674 scaling)':12.0,'NIRPS':9.5}
res={}
for nm,M in (('CK17',M_ck),('Otegi20',M_ot),('empirical',M_emp)):
    K=Kms(M,Ms); d=dict(M=pc(M),K=pc(K))
    for ins,s in SIG.items():
        for jit in (5,10,15):
            n=2*(5*np.hypot(s,jit)/K)**2
            d[f'N5sig_{ins}_jit{jit}']=[int(np.ceil(t)) for t in np.percentile(n,(16,50,84))]
    res[nm]=d
    print(nm, d['M'], d['K'], {k:v for k,v in d.items() if k.startswith('N5sig_ESPRESSO') or k.startswith('N5sig_NIRPS+')})
# TSM / ESM
J=11.36; Ks=10.48
Teq=Teff*np.sqrt(1/aR)*0.25**0.25
def tsm(M): return 1.15*Rp**3*Teq/(M*Rs**2)*10**(-J/5)
def planck(lam_um,T): return 1/np.expm1(14387.77/(lam_um*T))
esm=4.29e6*planck(7.5,1.10*Teq)/planck(7.5,Teff)*(Rp*6371/(Rs*695700))**2*10**(-Ks/5)
res['Teq']=pc(Teq); res['ESM']=pc(esm)
for nm,M in (('CK17',M_ck),('Otegi20',M_ot),('empirical',M_emp)): res[f'TSM_{nm}']=pc(tsm(M))
print('Teq',res['Teq'],'ESM',res['ESM'],'TSM',res['TSM_CK17'],res['TSM_empirical'])
# ranking
Rk=T[(T.pl_rade>=4)&(T.pl_rade<=7)&(T.pl_orbper<10)&T.sy_jmag.notna()].copy()
def teq_row(r):
    return r.pl_eqt if np.isfinite(r.pl_eqt) else np.nan
Rk['M_use']=Rk.pl_bmasse.where(Rk.pl_bmasselim.fillna(0)==0)
miss=Rk.M_use.isna()
Rk.loc[miss,'M_use']=2.04*(Rk.loc[miss,'pl_rade']/R_T1)**(1/0.589)
Rk['TSM']=1.15*Rk.pl_rade**3*Rk.pl_eqt/(Rk.M_use*Rk.st_rad**2)*10**(-Rk.sy_jmag/5)
Rk['ESM']=4.29e6*planck(7.5,1.10*Rk.pl_eqt)/planck(7.5,Rk.st_teff)*(Rk.pl_rade*6371/(Rk.st_rad*695700))**2*10**(-Rk.sy_kmag/5)
Rk=Rk.sort_values('TSM',ascending=False)
med_tsm=np.median(tsm(M_ck)); med_esm=np.median(esm)
rank_tsm=int((Rk.TSM>med_tsm).sum())+1; rank_esm=int((Rk.ESM>med_esm).sum())+1
print(Rk[['pl_name','pl_rade','M_use','pl_eqt','st_teff','sy_jmag','TSM','ESM']].head(25).to_string())
print('rank TSM', rank_tsm, 'of', len(Rk)+1, 'ESM', rank_esm, ' n with Teq', Rk.pl_eqt.notna().sum())
cool=Rk[Rk.st_teff<3400][['pl_name','st_teff','pl_rade']]
print('hosts < 3400 K:', cool.to_string())
res['ranking']=dict(rank_TSM=rank_tsm, rank_ESM=rank_esm, n_total=int(len(Rk)+1), n_missing_teq=int(Rk.pl_eqt.isna().sum()),
                    above_by_TSM=list(Rk[Rk.TSM>med_tsm].pl_name), cooler_than_3400K=list(cool.pl_name))
res['empirical_fit']=out
res['BD_K_kms']={f'{m}MJ': float(np.median(Kms(m*317.83, Ms))/1e3) for m in (13,20,40,60,80)}
res['WD0.6_K_kms']=float(Kms(0.6*MSUN/ME, Ms).mean()/1e3)
json.dump(res, open(f'{B}/out/r06_rv_tsm.json','w'), indent=1)
