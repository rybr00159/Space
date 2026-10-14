import glob, numpy as np
from astropy.io import fits
R='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/ctoi_vet/recent/data'
for f in sorted(glob.glob(R+'/mastDownload/HLSP/*71316387*/*.fits'))+sorted(glob.glob(R+'/qlp_71316387*.fits')):
    h=fits.open(f)
    h0=h[0].header; h1=h[1].header
    print(f.split('/')[-1], h0.get('SECTOR'), h0.get('CAMERA'), h0.get('CCD'), 'CROWD',h1.get('CROWDSAP'),'FLFR',h1.get('FLFRCSAP'), 'CDPP', h1.get('CDPP1_0'), 'TIMEDEL',h1.get('TIMEDEL'))
    print('   cols', h[1].columns.names)
    d=h[1].data; t=d['TIME']; print('   t', np.nanmin(t), np.nanmax(t), len(t))
