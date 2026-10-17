import sys,time; sys.path.insert(0,'tic214311908/analysis/referee_loc/scripts')
t=time.time()
import builtins
src=open('tic214311908/analysis/referee_loc/scripts/rq02_wcs_mira.py').read().split('for sec in SECTORS:')[0]
exec(src)
print('import',time.time()-t)
S=Sec(12); print('load',time.time()-t)
psf=S.prf(); O=np.median(S.flux[S.good],0); print('med',time.time()-t)
r=global_fit(S,O,psf); print('glob',time.time()-t, r['off'], r['nstar'], r['nfree'])
