import time; t=time.time()
from fpplib import *
print('import', time.time()-t); t=time.time()
for ds in ('ref_beta','inv_beta'):
    r=run(dict(ADOPT, ds=ds)); print(ds, {k:'%.3g'%v for k,v in r.items()}, time.time()-t); t=time.time()
