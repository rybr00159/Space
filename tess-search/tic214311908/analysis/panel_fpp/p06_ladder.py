"""p06: FPP ladder re-done with the ADOPTED sign-fixed composition (referee3_bd s12: kipping planet e, Ferreira short-period
BD e law, f_BD 1.3e-4; BD with the Gaia RV factor, and with RV + calibrated shape), using referee3_followup's own
bookkeeping rules and factors (r05_hri_ladder.py lines 63-90). Output out/p06_ladder.json"""
import json, numpy as np
B='tic214311908/analysis/'
s12=json.load(open(B+'referee3_bd/out/s12_summary.json'))
L=json.load(open(B+'referee3_followup/out/r05_hri_ladder.json'))
out={}
for ds,lc in (('ref_beta','ref_beta'),('inv_beta','inv_beta')):
    row=[r for r in s12['table'] if r['ds']==ds and r['planet_e']=='kipping' and r['bd_e']=='bdsp'][0]
    old=L['ladder'][lc]['composition_now']
    x=L['hri'][f'{lc}__signfixed']['Zorro832 cons B (faint, poor seeing; dI=1.2 at 0.04")']
    bdf=L['ladder'][lc]['bd_factor_from_zorro']
    # WD/HEB from s02 adopted row (same as ladder's to within 2%)
    import pandas as pd
    d=pd.read_csv(B+'referee3_bd/out/s02_fpp_fixed.csv'); a=d[(d.ds==ds)&(d.ecc=='kipping+bdsp')&(d.f_BD==1.3e-4)].iloc[0]
    res={}
    for bdlab,bdval in (('BD noRV',row['BD_total_noRV']),('BD RVmed (README table)',row['BD_total_RVmed']),('BD RVmed+shapecal',row['BD_total_RVmed_shapecal']),('ladder as published (f_BD 1e-3, BD e = planet e)',old['BD'])):
        comp=dict(PTP=row['PTP'],BD=bdval,WD=float(a.WD),HEB=float(a.HEB),BG=float(a.BG),NTP=row['NTP'])
        tot=lambda c: sum(c.values()); st={'now':tot(comp)}
        c=dict(comp); c['PTP']*=1-x; c['BD']*=bdf; c['HEB']*=1-x; st['Zorro consB']=tot(c)
        for lab in ('empirical g/z','empirical r/z','BT-Settl r/z'):
            for s in ('0.05','0.1'):
                bf=L['twoband'][f'{lc}__signfixed'][lab][f'BF_s{s}']
                c=dict(comp); c['PTP']*=bf; c['HEB']*=bf; c['NTP']=0; c['WD']=0
                c2=dict(c); c2['PTP']*=1-x; c2['HEB']*=1-x; c2['BD']*=bdf
                st[f'Zorro+two-band {lab} s{s}']=tot(c2)
                c3=dict(c2); c3['BD']=0; st[f'Zorro+two-band {lab} s{s} + BD closed']=tot(c3)
        res[bdlab]=dict(comp=comp,steps=st)
    out[ds]=res
json.dump(out,open('out/p06_ladder.json','w'),indent=1)
for ds in out:
    for bdlab in out[ds]:
        st=out[ds][bdlab]['steps']
        zt=[v for k,v in st.items() if k.startswith('Zorro+two-band') and 'closed' not in k]
        zc=[v for k,v in st.items() if 'closed' in k]
        print(f"{ds} {bdlab:48s} BD={out[ds][bdlab]['comp']['BD']:.2e} now {st['now']:.2e} Zorro {st['Zorro consB']:.2e} Z+2band {min(zt):.2e}-{max(zt):.2e}  +BDclosed {min(zc):.2e}-{max(zc):.2e}")
