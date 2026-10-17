import json
L=json.load(open('tic214311908/analysis/referee3_followup/out/r05_hri_ladder.json'))
P=json.load(open('out/p06_ladder.json'))
for ds in ('ref_beta','inv_beta'):
    for bdlab in ('BD RVmed (README table)','ladder as published (f_BD 1e-3, BD e = planet e)'):
        comp=P[ds][bdlab]['comp']; vals=[]
        for lab in ('empirical g/z','empirical r/z','BT-Settl r/z'):
            for s in ('0.05','0.1'):
                bf=L['twoband'][f'{ds}__signfixed'][lab][f'BF_s{s}']
                c=dict(comp); c['PTP']*=bf; c['HEB']*=bf; c['NTP']=0; c['WD']=0; vals.append(sum(c.values()))
        print(ds, bdlab, 'two-band only: %.2e - %.2e'%(min(vals),max(vals)))
