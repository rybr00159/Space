"""Proper-motion anomaly tests for the target vs an empirical null from Tycho-2/TGAS/Gaia DR3 field stars."""
import numpy as np, json
from astropy.table import Table
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/'
T=Table.read(D+'data/tyc_field.ecsv')
sid=6571343840265900544
def arr(c): return np.array(T[c].filled(np.nan) if hasattr(T[c],'filled') else T[c],float)
dec=np.radians(arr('dec')); cd=np.cos(dec)
ra3,de3=arr('ra'),arr('dec'); pmra3,pmde3=arr('pmra'),arr('pmdec')
# Tycho-2 mean position (J2000) -> back to mean epoch with Tycho-2 pm
raJ,deJ=arr('ra_mdeg'),arr('de_mdeg'); pmT_ra,pmT_de=arr('pm_ra'),arr('pm_de')
epra,epde=arr('ep_ra_m'),arr('ep_de_m')
eraM,edeM=arr('e_ra_mdeg'),arr('e_de_mdeg')   # mas at mean epoch
# offsets (mas) of DR3(2016.0) relative to Tycho J2000 mean position
dRA=(ra3-raJ)*3.6e6*cd; dDE=(de3-deJ)*3.6e6
# position at mean epoch relative to J2000 point: -pmT*(2000-ep)
dRA_m=dRA+pmT_ra*(2000-epra); dDE_m=dDE+pmT_de*(2000-epde)
mu_long_ra=dRA_m/(2016.0-epra); mu_long_de=dDE_m/(2016.0-epde)
e_long_ra=eraM/(2016.0-epra); e_long_de=edeM/(2016.0-epde)
# Tycho observed epoch (~1991.5) positions
ra91,de91=arr('ra_deg'),arr('de_deg'); ep91r=1990+arr('ep_ra1990'); ep91d=1990+arr('ep_de1990')
mu91_ra=(ra3-ra91)*3.6e6*cd/(2016-ep91r); mu91_de=(de3-de91)*3.6e6/(2016-ep91d)
e91_ra=arr('e_ra_deg')/(2016-ep91r); e91_de=arr('e_de_deg')/(2016-ep91d)
res={}
tests={'Tycho2_meanpos_to_DR3pos':(mu_long_ra,mu_long_de,e_long_ra,e_long_de),
       'Tycho2_obs1991_to_DR3pos':(mu91_ra,mu91_de,e91_ra,e91_de),
       'Tycho2_pm':(pmT_ra,pmT_de,arr('e_pm_ra'),arr('e_pm_de')),
       'TGAS_pm':(arr('tgas_pmra'),arr('tgas_pmdec'),arr('tgas_pmra_error'),arr('tgas_pmdec_error'))}
vt=arr('vt_mag'); ruwe=arr('ruwe'); G=arr('phot_g_mean_mag')
it=np.where(np.array(T['source_id'])==sid)[0][0]
print('target VT',vt[it],'G',G[it])
for nm,(a,b,ea,eb) in tests.items():
    da=pmra3-a; db=pmde3-b
    za=da/np.hypot(ea,arr('pmra_error')); zb=db/np.hypot(eb,arr('pmdec_error'))
    chi2=za**2+zb**2
    sel=np.isfinite(chi2)&(ruwe<1.4)&(np.abs(vt-vt[it])<0.75)
    sel2=np.isfinite(chi2)&(np.abs(vt-vt[it])<0.75)
    # robust scale of z in field
    sa=1.4826*np.nanmedian(np.abs(za[sel]-np.nanmedian(za[sel]))); sb=1.4826*np.nanmedian(np.abs(zb[sel]-np.nanmedian(zb[sel])))
    pct=np.mean(chi2[sel2]>=chi2[it])
    r=dict(dmu_ra=round(float(da[it]),3),dmu_de=round(float(db[it]),3),err_ra=round(float(np.hypot(ea[it],arr('pmra_error')[it])),3),err_de=round(float(np.hypot(eb[it],arr('pmdec_error')[it])),3),
           z_ra=round(float(za[it]),2),z_de=round(float(zb[it]),2),chi2=round(float(chi2[it]),2),
           field_n=int(sel.sum()),field_med_z=(round(float(np.nanmedian(za[sel])),2),round(float(np.nanmedian(zb[sel])),2)),field_robust_sigma_z=(round(float(sa),2),round(float(sb),2)),
           empirical_p_chi2_ge_target=round(float(pct),3),
           z_ra_cal=round(float((za[it]-np.nanmedian(za[sel]))/sa),2),z_de_cal=round(float((zb[it]-np.nanmedian(zb[sel]))/sb),2),
           field_dmu_robust_scatter_masyr=(round(float(1.4826*np.nanmedian(np.abs(da[sel]-np.nanmedian(da[sel])))),2),round(float(1.4826*np.nanmedian(np.abs(db[sel]-np.nanmedian(db[sel])))),2)))
    res[nm]=r; print(nm,r)
# DR2 vs DR3 (target only; DR2 pm from catalogue)
dr2=(38.425,-3.090,0.059,0.059); d3=(pmra3[it],pmde3[it],arr('pmra_error')[it],arr('pmdec_error')[it])
res['DR3_minus_DR2']=dict(dmu_ra=round(d3[0]-dr2[0],3),dmu_de=round(d3[1]-dr2[1],3),err=round(float(np.hypot(0.059,0.013)),3),note='DR2/DR3 not independent; DR2 bright-star frame spin ~0.1 mas/yr not included')
json.dump(res,open(D+'results/pma_tests.json','w'),indent=1)
