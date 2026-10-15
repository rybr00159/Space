"""SED analysis of TIC 71316387: single-star fit, G star + K/M pair grid, empirical twin null and injection-recovery.
Model: BT-Settl CIFIST synthetic photometry (SVO filters), CCM89 extinction, distance from Gaia DR3 parallax.
Pair components on the Mamajek (2022.04) dwarf sequence (Teff, R)."""
import numpy as np, json, sys
sys.path.insert(0,'/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/scripts')
import sedlib as S
from scipy.optimize import least_squares
from astropy.table import Table
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
rng=np.random.default_rng(7)
# ---------------- photometry ----------------
# band key -> (mag, err, system, floor, use_in_fit)
PH={'GALEX__GALEX.NUV':(16.8866,0.014,'AB',0.10,False),
 'SkyMapper__SkyMapper.u':(13.4222,0.0146,'AB',0.03,False),'SkyMapper__SkyMapper.v':(13.0779,0.0136,'AB',0.03,False),
 'SkyMapper__SkyMapper.g':(12.0213,0.0087,'AB',0.02,True),'SkyMapper__SkyMapper.r':(11.7585,0.0083,'AB',0.02,True),
 'SkyMapper__SkyMapper.i':(11.6204,0.0048,'AB',0.02,True),'SkyMapper__SkyMapper.z':(11.5991,0.0073,'AB',0.02,True),
 'Generic__Johnson.B':(12.493,0.019,'Vega',0.03,True),'Generic__Johnson.V':(11.857,0.026,'Vega',0.03,True),
 'SLOAN__SDSS.gprime_filter':(12.141,0.020,'AB',0.03,True),'SLOAN__SDSS.rprime_filter':(11.745,0.038,'AB',0.03,True),'SLOAN__SDSS.iprime_filter':(11.627,0.028,'AB',0.03,True),
 'TYCHO__TYCHO.B':(13.172,0.364,'Vega',0.0,True),'TYCHO__TYCHO.V':(11.781,0.169,'Vega',0.0,True),
 'GAIA__GAIA3.Gbp':(12.066102,0.002858,'Vega',0.02,True),'GAIA__GAIA3.G':(11.775332,0.002771,'Vega',0.02,True),'GAIA__GAIA3.Grp':(11.280263,0.003799,'Vega',0.02,True),
 'DENIS__DENIS.I':(11.171,0.04,'Vega',0.03,True),
 '2MASS__2MASS.J':(10.728,0.023,'Vega',0.02,True),'2MASS__2MASS.H':(10.389,0.024,'Vega',0.02,True),'2MASS__2MASS.Ks':(10.341,0.023,'Vega',0.02,True),
 'WISE__WISE.W1':(10.304,0.023,'Vega',0.03,True),'WISE__WISE.W2':(10.333,0.021,'Vega',0.03,True),'WISE__WISE.W3':(10.242,0.064,'Vega',0.03,True)}
B9=['GAIA__GAIA3.Gbp','GAIA__GAIA3.G','GAIA__GAIA3.Grp','2MASS__2MASS.J','2MASS__2MASS.H','2MASS__2MASS.Ks','WISE__WISE.W1','WISE__WISE.W2','WISE__WISE.W3']
FLOOR9={'GAIA__GAIA3.Gbp':0.02,'GAIA__GAIA3.G':0.02,'GAIA__GAIA3.Grp':0.02,'2MASS__2MASS.J':0.02,'2MASS__2MASS.H':0.02,'2MASS__2MASS.Ks':0.02,'WISE__WISE.W1':0.03,'WISE__WISE.W2':0.03,'WISE__WISE.W3':0.03}
SYS={k:v[2] for k,v in PH.items()}
ALLB=list(PH.keys())+['TESS__TESS.Red','Generic__Cousins.I','Paranal__VISTA.Ks']
for b in ['TESS__TESS.Red','Generic__Cousins.I','Paranal__VISTA.Ks']: SYS[b]='Vega'
# precompute surface-mag tables
TAB={}
for b in ALLB:
    T45,m45,_,_=S.G45; T50,m50,_,_=S.G50
    TAB[b]=(T45,np.array([m[b][SYS[b]] for m in m45]),T50,np.array([m[b][SYS[b]] for m in m50]),S.ccm_AlamAv(S.FILT[b]['leff'])*3.1)
def smag(b,teff):
    T45,v45,T50,v50,_=TAB[b]
    return np.interp(teff,T45,v45) if teff>=4500 else np.interp(teff,T50,v50)
def gmodel(bands,teff,R,d,ebv): return np.array([smag(b,teff)+S.dist_mod_R(R,d)+TAB[b][4]*ebv for b in bands])
# pair grid on Mamajek sequence
EEMr={sp:(T,R) for sp,T,R,_ in S.eem()}
PRIM=['K5V','K7V','M0V','M0.5V','M1V','M1.5V','M2V','M2.5V','M3V','M3.5V','M4V']
SEC=PRIM+['M4.5V','M5V','M5.5V']
PAIRS=[(a,None) for a in PRIM]+[(a,b) for i,a in enumerate(PRIM) for b in SEC[i:]]
def pair_abs(pair,bands):
    """absolute (10 pc) magnitudes of the pair in each band (no extinction)"""
    ms=[]
    for sp in pair:
        if sp is None: continue
        T,R=EEMr[sp]; ms.append(np.array([smag(b,T)+S.dist_mod_R(R,10.0) for b in bands]))
    return S.addmag(*ms)
def fit(bands,m,e,d,pair=None,ampfree=False,ebv0=0.009,sebv=0.01,x0=(5850,1.1,0.01)):
    A_ext=np.array([TAB[b][4] for b in bands])
    Mp=pair_abs(pair,bands) if pair is not None else None
    def model(p):
        mg=gmodel(bands,p[0],p[1],d,p[2])
        if Mp is None: return mg
        amp=p[3] if ampfree else 1.0
        if amp<=0: return mg
        mp=Mp+5*np.log10(d/10)+A_ext*p[2]-2.5*np.log10(amp)
        return S.addmag(mg,mp)
    def res(p):
        r=(m-model(p))/e
        return np.concatenate([r,[(p[2]-ebv0)/sebv]])
    lo=[4800,0.5,0.0]; hi=[6800,2.5,0.3]; p0=list(x0)
    if ampfree: lo+= [0.0]; hi+=[30.0]; p0+=[0.3]
    best=None
    for start in ([p0] if not ampfree else [p0[:3]+[0.01],p0[:3]+[0.5],p0[:3]+[2.0]]):
        r=least_squares(res,start,bounds=(lo,hi),x_scale=[100,0.05,0.01]+([0.3] if ampfree else []))
        if best is None or r.cost<best.cost: best=r
    return dict(p=best.x.tolist(),chi2=float(2*best.cost),resid=((m-model(best.x))).tolist())
def amp_profile(bands,m,e,d,pair,grid=np.array([0,0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0])):
    """chi2 profile vs pair amplitude (A=1: main sequence at the star's distance)"""
    out=[]
    for A in grid:
        if A==0: r=fit(bands,m,e,d)
        else:
            Mp=pair_abs(pair,bands); A_ext=np.array([TAB[b][4] for b in bands])
            def res(p):
                mg=gmodel(bands,p[0],p[1],d,p[2]); mp=Mp+5*np.log10(d/10)+A_ext*p[2]-2.5*np.log10(A)
                return np.concatenate([(m-S.addmag(mg,mp))/e,[(p[2]-0.009)/0.01]])
            rr=least_squares(res,[5850,1.1,0.01],bounds=([4800,0.5,0],[6800,2.5,0.3]),x_scale=[100,0.05,0.01]); r=dict(chi2=float(2*rr.cost))
        out.append(r['chi2'])
    return grid,np.array(out)
if __name__=='__main__':
    plx=3.3622; d=1000/plx
    out={}
    # ---- full-band single fit ----
    bands=[b for b,v in PH.items() if v[4]]
    m=np.array([PH[b][0] for b in bands]); e=np.array([np.hypot(PH[b][1],PH[b][3]) for b in bands])
    s1=fit(bands,m,e,d); out['single_allbands']=dict(Teff=s1['p'][0],R=s1['p'][1],EBV=s1['p'][2],chi2=s1['chi2'],nband=len(bands),
        resid={b:round(r,3) for b,r in zip(bands,s1['resid'])})
    allb=list(PH.keys()); mall=np.array([PH[b][0] for b in allb])
    pred=gmodel(allb,*s1['p'][:2],d,s1['p'][2]); out['single_allbands']['resid_all']={b:round(float(x),3) for b,x in zip(allb,mall-pred)}
    # bootstrap-ish uncertainty: refit with parallax +-1sigma and floors perturbed
    Ts=[];Rs=[]
    for k in range(200):
        dd=1000/(plx+rng.normal(0,0.017)); mm=m+rng.normal(0,e)
        r=fit(bands,mm,e,dd); Ts.append(r['p'][0]); Rs.append(r['p'][1])
    out['single_allbands']['Teff_std']=float(np.std(Ts)); out['single_allbands']['R_std']=float(np.std(Rs))
    print('single all-band:',json.dumps(out['single_allbands'],indent=0)[:1500],flush=True)
    # ---- pair grid on target (all bands) ----
    grid=[]
    for pr in PAIRS:
        r=fit(bands,m,e,d,pair=pr); ra=fit(bands,m,e,d,pair=pr,ampfree=True)
        # pair light fractions
        bb=['TESS__TESS.Red','GAIA__GAIA3.G','Generic__Cousins.I','2MASS__2MASS.Ks','WISE__WISE.W1']
        mg=gmodel(bb,r['p'][0],r['p'][1],d,r['p'][2]); mp=pair_abs(pr,bb)+5*np.log10(d/10)+np.array([TAB[b][4] for b in bb])*r['p'][2]
        frac={b.split('.')[-1]:float(1/(1+10**(-0.4*(mg[i]-mp[i])))) for i,b in enumerate(bb)}
        dm={b.split('.')[-1]:float(mp[i]-mg[i]) for i,b in enumerate(bb)}
        grid.append(dict(pair=[p for p in pr if p],chi2=r['chi2'],dchi2=r['chi2']-s1['chi2'],Teff_G=r['p'][0],R_G=r['p'][1],amp_best=ra['p'][3],chi2_ampfree=ra['chi2'],frac=frac,dmag=dm))
    out['pair_grid_allbands']=grid
    json.dump(out,open(D+'results/sed_target.json','w'),indent=1)
    gs=sorted(grid,key=lambda g:g['dchi2'])
    print('best pairs:',[(g['pair'],round(g['dchi2'],1),round(g['amp_best'],2),round(g['frac']['Red'],3)) for g in gs[:8]])
    print('worst pairs:',[(g['pair'],round(g['dchi2'],1)) for g in gs[-5:]])
