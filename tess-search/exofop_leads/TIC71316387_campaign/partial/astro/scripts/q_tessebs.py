import warnings; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
Vizier.TIMEOUT=200
v=Vizier(columns=['TIC','Tmag','Per','Morph','Dp-pf','Ds-pf','Wp-pf','Ws-pf','Dp-2g','Ds-2g'],row_limit=-1)
r=v.get_catalogs('J/ApJS/258/16/tess-ebs'); t=r[0]; print(len(t), t.colnames); t.write('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/tess_ebs_prsa2022.ecsv',overwrite=True)
