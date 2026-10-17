"""r05: (a) fraction of the PTP posterior excluded by high-resolution imaging under alternative contrast curves; (b) PTP
rejection by a two-band depth ratio under alternative colour predictions and per-event precisions; (c) FPP ladder by explicit
bookkeeping on the sign-fixed r05 scenario composition (followup/out/fpp_after_followup.csv row '0 none', which
referee3_bd reproduces: REF 1.51e-3, INV 3.59e-3 at n = 20000).
Inputs (data only): followup/out/ptp_draw_weights.npz (posterior-weighted companion draws: sep, dI, dV, dKs, dT; weights from
referee2_fpp r05, sign-fixed), followup/out/chromatic_ratio_vs_dT.csv (BT-Settl ratio vs dT), out/r08_chroma.json relations
are re-evaluated here from data/mdw_seq.csv.
Bayes factor of the two-band test for PTP: BF = <exp(-0.5((1-R_i)/sigma_R)^2)>_w (observed ratio 1, Gaussian), as r05 does.
Output: out/r05_hri_ladder.json"""
import os, json
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import numpy as np, pandas as pd
from scipy.interpolate import UnivariateSpline
B='tic214311908/analysis'
R=f'{B}/referee3_followup'
W=np.load(f'{B}/followup/out/ptp_draw_weights.npz')
sep=W['sep'].astype(float); dI=W['dI'].astype(float); dV=W['dV'].astype(float); dKs=W['dKs'].astype(float); dT=W['dT'].astype(float)
CURVES={
 'Zorro832 analyst nominal':('I',[(0.017,0.0),(0.025,2.5),(0.05,3.8),(0.1,4.3),(0.2,4.8),(0.5,5.2),(1.2,5.5)]),
 'Zorro832 cons A (dI=1.2 at 0.03")':('I',[(0.021,0.0),(0.03,1.2),(0.05,2.8),(0.1,3.8),(0.2,4.3),(0.5,4.8),(1.2,5.0)]),
 'Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")':('I',[(0.025,0.0),(0.04,1.2),(0.06,2.3),(0.1,3.0),(0.2,3.6),(0.5,4.0),(1.2,4.2)]),
 'Zorro832 cons C (dI=1.2 at 0.05")':('I',[(0.03,0.0),(0.05,1.2),(0.08,2.2),(0.1,2.6),(0.2,3.3),(0.5,3.8),(1.2,4.0)]),
 'HRCam median (analyst)':('I',[(0.04,0.0),(0.1,2.04),(0.2,4.04),(1.0,4.75),(1.5,5.28),(3.0,6.19)]),
 'HRCam 16th pct (analyst)':('I',[(0.04,0.0),(0.1,1.90),(0.2,3.43),(1.0,3.95),(1.5,4.31),(3.0,4.91)]),
 'HRCam pessimistic':('I',[(0.06,0.0),(0.1,1.4),(0.2,2.9),(1.0,3.45),(1.5,3.8),(3.0,4.4)]),
}
def lim(curve,s):
    a,v=map(np.array,zip(*curve)); return np.where(s>=a[0], np.interp(np.log10(s),np.log10(a),v), -99.)
band={'I':dI,'V':dV,'Ks':dKs}
res={'hri':{}, 'sep_cdf':{}}
for case in ('ref_beta__signfixed','inv_beta__signfixed'):
    w=W[case].astype(float)
    res['sep_cdf'][case]={str(s):round(float(np.sum(w*(sep<s))/w.sum()),3) for s in (0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.2,0.5)}
    res['hri'][case]={nm:round(float(np.sum(w*(band[b]<lim(c,sep)))/w.sum()),3) for nm,(b,c) in CURVES.items()}
    print(case, res['sep_cdf'][case]); print('  excluded', res['hri'][case])
# chromatic ratio per draw: BT-Settl table (analyst) and empirical relations (own)
T=pd.read_csv(f'{B}/followup/out/chromatic_ratio_vs_dT.csv').sort_values('dT')
bt={k: np.interp(dT, T.dT, T[f'ratio_{a}'])/np.interp(dT, T.dT, T[f'ratio_{b}']) for k,(a,b) in {'r/z':('r','z_s'),'g/z':('g','z_s'),'V/Ic':('V','Ic')}.items()}
d=pd.read_csv(f'{R}/data/mdw_seq.csv'); d=d[(d.MKs>5.3)&(d.MKs<9.3)]; dm=5*np.log10(d.parallax/100.)
xs=np.linspace(5.6,8.8,33); rel={}
for b,c in {'g':'g_sdss_mag','r':'r_sdss_mag','z':'z_sdss_mag','V':'v_jkc_mag','Ic':'i_jkc_mag'}.items():
    y=(d[c]+dm).values; x=d.MKs.values; ok=np.isfinite(y)
    med=np.array([np.median(y[ok&(np.abs(x-t)<0.25)]) for t in xs]); n=np.array([np.sum(ok&(np.abs(x-t)<0.25)) for t in xs])
    rel[b]=UnivariateSpline(xs[n>=5],med[n>=5],k=3,s=len(xs)*0.003)
MKp=6.441; MKc=np.clip(MKp+dKs,5.6,8.8)
f={b:1/(1+10**(0.4*(rel[b](MKc)-rel[b](MKp)))) for b in rel}
em={'r/z':f['r']/f['z'],'g/z':f['g']/f['z'],'V/Ic':f['V']/f['Ic']}
# r05/referee2 (Mann-based) V/I from dV, dI of the draws
mann_VI=(1/(1+10**(0.4*dV)))/(1/(1+10**(0.4*dI)))
res['twoband']={}
for case in ('ref_beta__signfixed','inv_beta__signfixed'):
    w=W[case].astype(float); out={}
    for lab,Rv in (('BT-Settl r/z',bt['r/z']),('BT-Settl g/z',bt['g/z']),('empirical r/z',em['r/z']),('empirical g/z',em['g/z']),
                   ('empirical V/Ic',em['V/Ic']),('referee2 Mann V/I (used inside r05 chroma)',mann_VI)):
        q=np.interp([0.16,0.5,0.84,0.975],np.cumsum(w[np.argsort(Rv)])/w.sum(),np.sort(Rv))
        row=dict(pct16_50_84_975=[round(float(t),3) for t in q])
        for s in (0.05,0.07,0.10,0.14):
            row[f'BF_s{s}']=float(np.sum(w*np.exp(-0.5*((1-Rv)/s)**2))/w.sum())
            row[f'rej3sig_s{s}']=round(float(np.sum(w*((1-Rv)/s>3))/w.sum()),3)
        out[lab]=row
    res['twoband'][case]=out
    for k,v in out.items(): print(case,k,v)
# ---- FPP ladder bookkeeping (sign-fixed composition, analyst's r05 row '0 none')
F=pd.read_csv(f'{B}/followup/out/fpp_after_followup.csv')
lad={}
for lc,case in (('ref_beta','ref_beta__signfixed'),('inv_beta','inv_beta__signfixed')):
    r0=F[(F.sign=='fixed')&(F.lc==lc)&(F.step=='0 none')].iloc[0]
    r1a=F[(F.sign=='fixed')&(F.lc==lc)&(F.step.str.startswith('1a'))].iloc[0]
    comp=dict(PTP=r0.PTP, BD=r0.BD, WD=r0.WD, HEB=r0.HEB, BG=r0.BG, NTP=r0.NFPP)
    # BD on companion share removed by Zorro in the analyst's own r05 run (BD 1a vs 0)
    bd_hri_factor = r1a.BD/r0.BD
    tot=lambda c: sum(c.values())
    steps={'now':tot(comp)}
    for zname in ('Zorro832 analyst nominal','Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")','Zorro832 cons C (dI=1.2 at 0.05")'):
        x=res['hri'][case][zname]
        c=dict(comp); c['PTP']*=1-x; c['BD']*=bd_hri_factor; c['HEB']*=1-x
        steps[f'speckle: {zname}']=tot(c)
    for lab in ('BT-Settl r/z','empirical r/z'):
        for s in (0.05,0.10):
            bf=res['twoband'][case][lab][f'BF_s{s}']
            c=dict(comp); c['PTP']*=bf; c['HEB']*=bf; c['NTP']=0; c['WD']=0
            steps[f'two-band {lab} sigma {s}']=tot(c)
            x=res['hri'][case]['Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")']
            c2=dict(c); c2['PTP']*=1-x; c2['HEB']*=1-x; c2['BD']*=bd_hri_factor
            steps[f'Zorro consB + two-band {lab} sigma {s}']=tot(c2)
            c3=dict(c2); c3['BD']=0
            steps[f'Zorro consB + two-band {lab} sigma {s} + low-res RV (BD closed)']=tot(c3)
    # f_BD = 1.3e-4 variant (bd/referee3_bd): BD x 0.13
    c=dict(comp); c['BD']*=0.13; steps['now, f_BD 1.3e-4']=tot(c)
    lad[lc]=dict(composition_now=comp, bd_factor_from_zorro=bd_hri_factor, steps=steps)
    for k,v in steps.items(): print(lc, f'{k:85s} {v:.2e}')
res['ladder']=lad
json.dump(res, open(f'{R}/out/r05_hri_ladder.json','w'), indent=1, default=float)
