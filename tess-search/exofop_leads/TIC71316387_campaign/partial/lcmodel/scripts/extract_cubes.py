"""Copy an 11x11 sub-cube (centred on the target) of each TESScut 21x21 cutout into a compact npz,
so that our photometry does not depend on the sibling folder's large files staying on disk."""
import glob, numpy as np, sys
from astropy.io import fits
from astropy.wcs import WCS
SRC='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut'
OUT='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/data'
RA,DEC=330.887948,-42.835511
import os; os.makedirs(OUT,exist_ok=True)
for fn in sorted(glob.glob(SRC+'/s*/*.fits')):
    h=fits.open(fn, memmap=True)
    d=h[1].data; w=WCS(h[2].header)
    x0,y0=w.all_world2pix(RA,DEC,0); ix,iy=int(round(float(x0))),int(round(float(y0)))
    hw=5
    sl=(slice(None),slice(iy-hw,iy+hw+1),slice(ix-hw,ix+hw+1))
    F=np.array(d['FLUX'][sl],np.float32); Fe=np.array(d['FLUX_ERR'][sl],np.float32)
    # full-frame median image for background estimation using the whole 21x21 stamp
    full=np.array(d['FLUX'],np.float32)
    t=np.array(d['TIME'],float); q=np.array(d['QUALITY'])
    sec=int(h[0].header['SECTOR'])
    # store background proxies: median of the faintest 25% pixels of full stamp per cadence
    med=np.nanmedian(full,axis=0)
    faint=med<np.nanpercentile(med,25)
    bk=np.nanmedian(full[:,faint],axis=1)
    bk2=np.nanpercentile(full.reshape(len(t),-1),10,axis=1)
    hdr2=h[2].header
    # WCS of subcube: shift CRPIX
    wsub=w.deepcopy(); wsub.wcs.crpix=wsub.wcs.crpix-np.array([ix-hw,iy-hw])
    np.savez_compressed(f'{OUT}/cube_s{sec}.npz',t=t,q=q,F=F,Fe=Fe,bk=bk,bk10=bk2,ix=ix,iy=iy,x0=float(x0),y0=float(y0),
        med_full=med,cam=h[0].header['CAMERA'],ccd=h[0].header['CCD'],col0=h[1].header.get('1CRV4P',0),row0=h[1].header.get('2CRV4P',0),
        wcs_sub=wsub.to_header_string(),wcs_full=w.to_header_string(),timedel=h[1].header.get('TIMEDEL'))
    print(sec, fn.split('/')[-1], len(t), 'target pix', x0, y0, 'col0,row0', h[1].header.get('1CRV4P'), h[1].header.get('2CRV4P'))
