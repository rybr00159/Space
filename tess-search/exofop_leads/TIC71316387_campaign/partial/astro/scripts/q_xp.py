import sys, warnings, os; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
ids=[int(x) for x in open(sys.argv[1]).read().split()]
out=sys.argv[2]
ids=ids[:450]
for i in range(0,len(ids),50):
    dl=Gaia.load_data(ids=ids[i:i+50],data_release='Gaia DR3',retrieval_type='XP_SAMPLED',data_structure='DATAMODEL_GAIA',format='votable',output_file=None)
    for k,v in dl.items():
        for j,t in enumerate(v):
            tb=t.to_table(); tb.write(f"{out}_{i}_{j}.ecsv",overwrite=True); print(k,len(tb),tb.colnames[:8],flush=True)
