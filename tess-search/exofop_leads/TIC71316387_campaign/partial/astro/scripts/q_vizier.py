import sys, warnings
warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
out=sys.argv[1]
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=60
v=Vizier(columns=['**','+_r'],row_limit=-1)
cats={'galex_guv':'II/335/galex_ais','galex_ais':'II/312/ais','apass9':'II/336/apass9','sm_dr2':'II/358/smss','sm_dr4':'II/379/smssdr4',
'gaia3':'I/355/gaiadr3','gspc':'I/360/syntphot','tyc2':'I/259/tyc2','tmass':'II/246/out','vhs':'II/367/vhs_dr5','allwise':'II/328/allwise',
'catwise':'II/365/catwise','unwise':'II/363/unwise','denis':'B/denis/denis','ucac4':'I/322A/out','ucac5':'I/340/ucac5','rave6':'III/283',
'galah3':'J/MNRAS/506/150','apogee17':'III/286','tic82':'IV/39/tic82','ppmxl':'I/317/sample','urat':'I/329','hsoy':'I/339/hsoy','ucac3':'I/315/out','nomad':'I/297/out','ascc':'I/280B/ascc','gaia2':'I/345/gaia2','dr1':'I/337/gaia','tgas':'I/337/tgas','spm4':'I/327/cat','galah4':'J/PASA/42/E051','rave5':'III/279','denisw':'B/denis','sdss':'V/154','des':'II/371/des_dr2','ps1':'II/349','gaiaesp':'I/356','skym4':'II/379'}
r={'galex_guv':10,'galex_ais':10}
for k,cat in cats.items():
    try:
        res=v.query_region(c,radius=r.get(k,8)*u.arcsec,catalog=cat)
        if len(res)==0: print(k,cat,'none'); continue
        for i,t in enumerate(res):
            fn=f"{out}/vz_{k}_{i}.ecsv"; t.write(fn,overwrite=True)
            print(k,cat,t.meta.get('name',''),len(t))
    except Exception as e: print(k,cat,'ERR',str(e)[:200])
