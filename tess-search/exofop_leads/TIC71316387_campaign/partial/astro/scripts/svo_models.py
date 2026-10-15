import subprocess, os
from astropy.io.votable import parse
S='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/svo/'
t=parse(S+'list_bts.xml').get_first_table().to_table()
n=0
for r in t:
    te,lg=float(r['teff']),float(r['logg'])
    if not (2600<=te<=7000 and lg in (4.5,5.0)): continue
    fn=S+f'btsettl/bts_{int(te):05d}_{lg:.1f}.xml'
    if os.path.exists(fn) and os.path.getsize(fn)>100000: continue
    subprocess.run(['curl','-sS','--max-time','120','-o',fn,str(r['Access.Reference'])]); n+=1
print('downloaded',n)
