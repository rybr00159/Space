import pandas as pd, numpy as np, json
B='tic214311908/analysis/'
d=pd.read_csv(B+'referee3_bd/out/s02_fpp_fixed.csv')
f=d[(d.sign=='fixed')&d.f_BD.notna()].copy()
f['plaw']=f.ecc.str.split('+').str[0]; f['bdlaw']=f.ecc.apply(lambda s: s.split('+')[1] if '+' in s else 'same')
print('n fixed rows with f_BD', len(f)); print(f.groupby('bdlaw').size())
f['NFPP_tri']=f.NTP_pairA            # TRICERATOPS-like: resolved neighbours only
f['NFPP_wbg']=f.NTP_pairA+f.BG
f['FPP_strict']=f.FPP+f.TP_twin       # twin-hosted planet counted as FP (TRICERATOPS STP convention, upper bound)
f['P_any_incl_pairA']=f.P_planet_any_star+f.NTP_pairA
f['sumcheck']=1-(f.TP_single+f.TP_diluted+f.TP_twin)-f.FPP
f['sys_check']=f.FPP_system-(f.BD_total+f.NTP_pairA+f.WD+f.HEB+f.BG)
print('max |sumcheck|',f.sumcheck.abs().max(),' max |sys_check|',f.sys_check.abs().max())
out={}
for name,sub in [('all192',f),('bdsp64',f[f.bdlaw=='bdsp']),('same64',f[f.bdlaw=='same']),('bdlp64',f[f.bdlaw=='bdlp'])]:
    r={c:[float(sub[c].min()),float(sub[c].max())] for c in ['FPP','P_planet_any_star','NTP_pairA','NFPP_wbg','PTP','BD_total','FPP_strict','P_any_incl_pairA','WD']}
    out[name]=r
    print(name, {k:['%.3g'%v[0],'%.3g'%v[1]] for k,v in r.items()})
# rows where NFPP >1e-3
print(f[f.NFPP_wbg>1e-3][['ds','ecc','f_BD','NTP_pairA','NFPP_wbg','FPP']].to_string())
# adopted
ad=f[(f.plaw.isin(['kipping','rayleigh19','circ']))&(f.bdlaw=='bdsp')&(f.f_BD==1.3e-4)&(f.ds.isin(['ref_beta','inv_beta']))]
print(ad[['ds','plaw','FPP','FPP_strict','PTP','NTP_pairA','NFPP_wbg','BD_total','WD','BG','HEB','TP_twin','P_planet_any_star','P_any_incl_pairA']].to_string(float_format=lambda x:'%.3g'%x))
# PTP circular/eccentric ratio
for ds in ['ref_beta','inv_beta','ref_white','inv_white']:
    s=f[(f.ds==ds)&(f.bdlaw=='bdsp')&(f.f_BD==1.3e-4)].set_index('plaw')
    print(ds,'PTP kip/circ %.1f  ray/circ %.1f'%(s.PTP['kipping']/s.PTP['circ'], s.PTP['rayleigh19']/s.PTP['circ']))
json.dump(out,open('out/p02_ranges.json','w'),indent=1)
