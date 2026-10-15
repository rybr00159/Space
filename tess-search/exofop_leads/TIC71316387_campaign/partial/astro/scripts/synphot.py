"""Convert BT-Settl CIFIST surface-flux spectra (SVO) into (a) synthetic magnitudes at the stellar surface for all filters,
(b) a compact low-resolution spectrum (300-1100 nm at 1 nm; 1.1-30 um coarse) ; then delete the big XML."""
import numpy as np, glob, os, re, json, sys
from astropy.io.votable import parse
S='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/svo/'
O='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/models/'
c_A=2.99792458e18
filt={}
for fn in sorted(glob.glob(S+'filters/*_Vega.xml')):
    name=os.path.basename(fn).replace('_Vega.xml','')
    out={}
    for sysn in ['Vega','AB']:
        v=parse(fn.replace('_Vega',f'_{sysn}')); t=v.get_first_table(); p={q.name:q.value for q in t.params}
        tb=t.to_table(); out['wl']=np.array(tb['Wavelength'],float); out['tr']=np.array(tb['Transmission'],float)
        out['zp_'+sysn]=float(p['ZeroPoint']); out['det']=int(p.get('DetectorType',1) or 1); out['lpiv']=float(p['WavelengthPivot']); out['leff']=float(p['WavelengthEff'])
    filt[name]=out
np.save(O+'filters.npy',filt,allow_pickle=True)
print('filters',list(filt))
grid_wl=np.concatenate([np.arange(3000,11000,10.),np.geomspace(11000,300000,800)])
for fn in sorted(glob.glob(S+'btsettl/bts_*.xml')):
    m=re.search(r'bts_(\d+)_([\d.]+)\.xml',fn); te,lg=int(m.group(1)),float(m.group(2))
    of=O+f'bts_{te:05d}_{lg:.1f}.npz'
    if os.path.exists(of): os.remove(fn); continue
    tb=parse(fn).get_first_table().to_table()
    wl=np.array(tb['WAVELENGTH'],float); fl=np.array(tb['FLUX'],float)
    ok=(wl>0); wl,fl=wl[ok],fl[ok]; o=np.argsort(wl); wl,fl=wl[o],fl[o]
    u,iu=np.unique(wl,return_index=True); wl,fl=u,fl[iu]
    fbol=np.trapezoid(fl,wl); teff_chk=(fbol/5.670374e-5)**0.25
    mags={}
    for name,f in filt.items():
        T=np.interp(wl,f['wl'],f['tr'],left=0,right=0); s=T>0
        if s.sum()<3: continue
        if f['det']==1: num=np.trapezoid(fl[s]*wl[s]*T[s],wl[s]); den=np.trapezoid(wl[s]*T[s],wl[s])
        else: num=np.trapezoid(fl[s]*T[s],wl[s]); den=np.trapezoid(T[s],wl[s])
        flam=num/den; fnu_jy=flam*f['lpiv']**2/c_A*1e23
        mags[name]=dict(Vega=-2.5*np.log10(fnu_jy/f['zp_Vega']),AB=-2.5*np.log10(fnu_jy/f['zp_AB']))
    # binned spectrum (mean over bins)
    edges=np.concatenate([[grid_wl[0]-(grid_wl[1]-grid_wl[0])/2],(grid_wl[1:]+grid_wl[:-1])/2,[grid_wl[-1]*1.01]])
    cs=np.concatenate([[0],np.cumsum(0.5*(fl[1:]+fl[:-1])*np.diff(wl))]); cw=np.concatenate([[0],np.cumsum(np.diff(wl))])
    Fi=np.interp(edges,wl,cs); Wi=np.interp(edges,wl,cw)
    spec=np.diff(Fi)/np.maximum(np.diff(edges),1e-9)
    np.savez(of,teff=te,logg=lg,teff_chk=teff_chk,fbol=fbol,mags=json.dumps(mags),wl=grid_wl,flux=spec)
    os.remove(fn)
    print(te,lg,'Teff_chk',round(teff_chk),flush=True)
