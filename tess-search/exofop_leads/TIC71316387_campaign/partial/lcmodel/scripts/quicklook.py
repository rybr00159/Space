import numpy as np, sys, json
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import lclib as L
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def event_box(lc, ws):
    res=[]
    for w in ws:
        i=w['idx']; dt=lc['t'][i]-w['tc']; y=lc['f'][i]
        oot=np.abs(dt)>0.75*L.T14
        c=np.polyfit(dt[oot],y[oot],2); base=np.polyval(c,dt)
        core=np.abs(dt)<0.2*L.T14
        if core.sum()<2: continue
        res.append((w['ph'],1-np.mean(y[core]/base[core]),w['n']))
    return np.array(res)
out={}
lcs={}
for sec in (1,28,68,95,105,106):
    prods=[]
    if sec in (1,28,68): prods.append(L.load_spoc(sec))
    prods.append(L.load_tc(sec))
    if sec in (1,28,68,95,105): prods.append(L.load_qlp(sec))
    for lc in prods:
        lc=L.clip_outliers(lc)
        ws=L.windows(lc)
        r=event_box(lc,ws)
        line=f"S{sec} {lc['prod']:9s} pt-sigma {lc['sig_pt']*1e6:5.0f} ppm clip {lc['nclip']:3d} "
        for ph in (0.0,0.5):
            d=r[r[:,0]==ph,1]
            line+=f" ph{ph}: N={len(d):2d} depth {np.mean(d)*1e6:6.0f}+-{np.std(d,ddof=1)/np.sqrt(len(d))*1e6:4.0f}"
        out[f"{sec}_{lc['prod']}"]=dict(prim=float(np.mean(r[r[:,0]==0,1])),prim_e=float(np.std(r[r[:,0]==0,1],ddof=1)/np.sqrt((r[:,0]==0).sum())),
            sec=float(np.mean(r[r[:,0]==0.5,1])),sec_e=float(np.std(r[r[:,0]==0.5,1],ddof=1)/np.sqrt((r[:,0]==0.5).sum())))
        print(line, flush=True)
        lcs[(sec,lc['prod'])]=lc
json.dump(out,open(L.HERE+'/results/quicklook_depths.json','w'),indent=1)
# folded plot of default set
fig,ax=plt.subplots(2,1,figsize=(9,8))
for sec,prod,col in [(1,'SPOC-PDC','C0'),(28,'SPOC-PDC','C1'),(68,'SPOC-PDC','C2'),(95,'TC-AP','C3'),(105,'TC-AP','C4'),(106,'TC-AP','C5')]:
    lc=lcs[(sec,prod)]
    for k,ph in enumerate((0.0,0.5)):
        ws=L.windows(lc,phases=(ph,))
        X=[];Y=[]
        for w in ws:
            i=w['idx']; dt=lc['t'][i]-w['tc']; y=lc['f'][i]; oot=np.abs(dt)>0.75*L.T14
            c=np.polyfit(dt[oot],y[oot],2); X.append(dt); Y.append(y/np.polyval(c,dt))
        X=np.concatenate(X)*24; Y=np.concatenate(Y)
        b=np.arange(-4.8,4.81,0.2); idx=np.digitize(X,b)
        bx=[X[idx==j].mean() for j in range(1,len(b)) if (idx==j).sum()>3]; by=[Y[idx==j].mean() for j in range(1,len(b)) if (idx==j).sum()>3]
        ax[k].plot(bx,by,'.-',color=col,label=f'S{sec} {prod}')
for k in range(2): ax[k].set_xlabel('hours from predicted eclipse'); ax[k].legend(fontsize=7); ax[k].set_title(['primary (phase 0)','secondary (phase 0.5)'][k])
plt.tight_layout(); plt.savefig(L.HERE+'/plots/quicklook_fold.png',dpi=110)
