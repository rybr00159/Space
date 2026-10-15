import numpy as np
from astropy.table import Table, join
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/'
G=Table.read(D+'twins_gaia.ecsv'); X=Table.read(D+'twins_xm.ecsv')
T=join(G,X,keys='source_id',join_type='left')
for c in T.colnames:
    if hasattr(T[c],'filled'):
        k=T[c].dtype.kind
        T[c]=T[c].filled(np.nan if k=='f' else (-1 if k in 'iu' else ('' if k in 'US' else False)))
T['MG']=T['phot_g_mean_mag']+5*np.log10(T['parallax']/100)
sid=6571343840265900544
it=np.where(T['source_id']==sid)[0][0]; tg=T[it]
print('target',tg['bp_rp'],tg['MG'],tg['tm_qual'],tg['aw_qual'],tg['cc_flags'],tg['ipd_frac_multi_peak'])
qual=(T['tm_qual']=='AAA')&np.array([str(q)[:3]=='AAA' or str(q)[:2]=='AA' for q in T['aw_qual']])&(T['cc_flags']=='0000')&(T['ext_flag']==0)
clean=(T['ruwe']<1.3)&(T['ipd_frac_multi_peak']<=5)
sel=qual&clean&(np.abs(T['bp_rp']-tg['bp_rp'])<0.05)&(np.abs(T['MG']-tg['MG'])<0.45)&(T['source_id']!=sid)
print('quality',qual.sum(),'clean',(qual&clean).sum(),'twins',sel.sum())
T.write(D+'twins_merged.ecsv',overwrite=True)
S=T[sel]; S.write(D+'twins_sel.ecsv',overwrite=True)
# XP subset
xs=S[S['has_xp_sampled']==True]
print('with XP',len(xs))
open(D+'xp/ids_twins.txt','w').write('\n'.join(str(x) for x in xs['source_id'][:800]))
