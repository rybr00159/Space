import json, glob, numpy as np, sys
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/tri/'
EBs=['EB','EBx2P','PEB','PEBx2P','SEB','SEBx2P','DEB','DEBx2P','BEB','BEBx2P']
out={}
for grp in sys.argv[1:]:
    fs=sorted(glob.glob(D+grp+'*.json')); rows=[]
    for f in fs:
        d=json.load(open(f)); sc={s['scenario']:s for s in d['scenarios'] if str(s['ID'])=='71316387'}
        lnZ={k:(sc[k]['lnZ'] if k in sc else -np.inf) for k in EBs}
        L=np.array([lnZ[k] for k in EBs]); L=np.where(np.isfinite(L),L,-1e9); p=np.exp(L-L.max()); p/=p.sum()
        r={k:float(v) for k,v in zip(EBs,p)}
        r['lnZ_SEBx2P']=lnZ['SEBx2P']; r['lnZ_BEBx2P']=lnZ['BEBx2P']; r['lnZ_SEB']=lnZ['SEB']; r['lnZ_BEB']=lnZ['BEB']
        r['bound_frac']=r['SEB']+r['SEBx2P']+r['PEB']+r['PEBx2P']; r['bg_frac']=r['BEB']+r['BEBx2P']+r['DEB']+r['DEBx2P']
        r['tdepth']=d['tdepth']; r['FPP_raw']=d['FPP']
        best={k:dict(M_s=sc[k]['M_s'],R_s=sc[k]['R_s'],inc=sc[k]['inc'],b=sc[k]['b'],M_EB=sc[k]['M_EB'],R_EB=sc[k]['R_EB'],fr_EB=sc[k].get('fluxratio_EB'),fr_comp=sc[k].get('fluxratio_comp')) for k in ['SEBx2P','BEBx2P','SEB','BEB'] if k in sc}
        r['best']=best; rows.append(r)
    if not rows: continue
    dl=np.array([r['lnZ_SEBx2P']-r['lnZ_BEBx2P'] for r in rows])
    summ=dict(n=len(rows),files=[f.split('/')[-1] for f in fs],
        median_probs={k:float(np.median([r[k] for r in rows])) for k in EBs},
        min_max_probs={k:(float(np.min([r[k] for r in rows])),float(np.max([r[k] for r in rows]))) for k in EBs},
        dlnZ_SEBx2P_minus_BEBx2P=dict(values=[round(float(x),2) for x in dl],median=float(np.median(dl)),std=float(np.std(dl))),
        bound_frac=[round(r['bound_frac'],4) for r in rows], bg_frac=[round(r['bg_frac'],4) for r in rows],
        # evidence pooled over seeds: mean of Z across runs (each run is an unbiased MC estimate of Z)
        pooled={},tdepth=rows[0]['tdepth'],best_example=rows[0]['best'])
    for k in EBs:
        pass
    Zs={}
    for f in fs:
        d=json.load(open(f)); sc={s['scenario']:s for s in d['scenarios'] if str(s['ID'])=='71316387'}
        for k in EBs:
            Zs.setdefault(k,[]).append(sc[k]['lnZ'] if k in sc and sc[k]['lnZ'] is not None else -np.inf)
    lz={k:(np.log(np.mean(np.exp(np.array(v,float)-np.max(v))))+np.max(v) if np.isfinite(np.max(v)) else -np.inf) for k,v in Zs.items()}
    L=np.array([lz[k] for k in EBs]); L=np.where(np.isfinite(L),L,-1e9); p=np.exp(L-L.max()); p/=p.sum()
    summ['pooled']={k:float(v) for k,v in zip(EBs,p)}
    summ['pooled_bound']=summ['pooled']['SEB']+summ['pooled']['SEBx2P']+summ['pooled']['PEB']+summ['pooled']['PEBx2P']
    summ['pooled_bg']=summ['pooled']['BEB']+summ['pooled']['BEBx2P']+summ['pooled']['DEB']+summ['pooled']['DEBx2P']
    out[grp]=summ
    print('=====',grp,'n=',len(rows),'tdepth',round(summ['tdepth'],5))
    print(' pooled:',{k:round(v,4) for k,v in summ['pooled'].items() if v>1e-4},' bound',round(summ['pooled_bound'],4),'bg',round(summ['pooled_bg'],4))
    print(' per-run SEBx2P:',[round(r['SEBx2P'],3) for r in rows],' BEBx2P:',[round(r['BEBx2P'],3) for r in rows])
    print(' dlnZ(SEBx2P-BEBx2P):',summ['dlnZ_SEBx2P_minus_BEBx2P'])
json.dump(out,open(D+'tri_summary_all.json','w'),indent=1,default=str)
