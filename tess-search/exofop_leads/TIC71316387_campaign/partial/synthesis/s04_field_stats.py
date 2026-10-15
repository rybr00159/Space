import numpy as np, json
from astropy.table import Table
B='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
t=Table.read(B+'synthesis/field_Mdwarfs_1p5deg.ecsv')
ra0,de0=330.887948,-42.835511
ra=np.radians(np.array(t['ra'])); de=np.radians(np.array(t['dec']))
sep=np.degrees(np.arccos(np.clip(np.sin(np.radians(de0))*np.sin(de)+np.cos(np.radians(de0))*np.cos(de)*np.cos(ra-np.radians(ra0)),-1,1)))
k=sep>1/60.
br=np.array(t['bp_rp'],float); G=np.array(t['phot_g_mean_mag'],float)
Tm=G-0.00522555*br**3+0.0891337*br**2-0.633923*br+0.0324473
f=10**(-0.4*(Tm-11.358))
plx=np.array(t['parallax'].filled(np.nan),float); pmra=np.array(t['pmra'].filled(np.nan),float); pmde=np.array(t['pmdec'].filled(np.nan),float)
area=np.pi*(1.5**2-(1/60)**2)*3600**2
circ=np.pi*1.24**2
r9=json.load(open(B+'referee_all/r09_msprior.json'))
out={}
dmu=np.hypot(pmra-38.503,pmde+3.134)
for lab,(c1,c2) in {'2.1-2.5':(2.1,2.5),'1.8-2.7':(1.8,2.7)}.items():
    m=k&(br>=c1)&(br<=c2)
    w=m&(f>=0.022)&(f<=0.065)
    o={'n_window':int(w.sum()),'N_per_1.24circle_window':float(w.sum()/area*circ)}
    for s in ('A_Mann','B_PM_table'):
        tab=r9[s]; ff=np.array([x['f'] for x in tab]); tot=np.array([x['total'] for x in tab]); tmin=tot.min()
        Lf=np.exp(-0.5*(np.interp(np.log(np.clip(f[m],1e-4,1)),np.log(ff),tot,left=1e3,right=1e3)-tmin))
        o[f'Neff_Lweighted_{s}']=float(np.nansum(Lf)); o[f'Neff_per_circle_{s}']=float(np.nansum(Lf)/area*circ)
    ww=w&np.isfinite(dmu)
    o['frac_dmu_lt3']=float((dmu[ww]<3).mean()) if ww.sum() else None
    o['frac_dmu_lt5']=float((dmu[ww]<5).mean()) if ww.sum() else None
    o['n_dmu_lt3']=int((dmu[ww]<3).sum())
    o['frac_dmu_lt3_and_plx_within_1mas']=float(((dmu[ww]<3)&(np.abs(plx[ww]-3.3622)<1.0)).mean()) if ww.sum() else None
    o['dist_pc_quantiles_window']=[float(x) for x in np.nanpercentile(1000/np.clip(plx[w],0.05,None),[5,16,50,84,95])]
    o['frac_window_265_345pc']=float(((1000/np.clip(plx[w],0.05,None)>265)&(1000/np.clip(plx[w],0.05,None)<345)).mean())
    o['dmu_quantiles_window']=[float(x) for x in np.nanpercentile(dmu[ww],[5,16,50,84,95])]
    out[lab]=o
# bound q-window for the wider colour window: pair mass ~0.50-0.80 vs 0.57-0.72 -> (0.80-0.50)/0.9 vs 0.15/0.9
out['note']='like-for-like odds scale as q_window/N; widen colour window 1.8-2.7 <-> q_tot 0.48-0.78 (0.30/0.9)'
r6=json.load(open(B+'referee_all/r06_priors.json'))
P0=0.58*0.2716531478730537
for lab,qw in {'2.1-2.5':0.15/0.9,'1.8-2.7':0.30/0.9}.items():
    o=out[lab]
    out[lab]['odds_hardwindow']=P0*qw/o['N_per_1.24circle_window'] if o['n_window'] else None
    for s in ('A_Mann','B_PM_table'):
        out[lab][f'odds_Lweighted_{s}_noLb']=P0*qw/o[f'Neff_per_circle_{s}']
print(json.dumps(out,indent=1))
json.dump(out,open(B+'synthesis/s04_field_stats.json','w'),indent=1)
