from astroquery.vizier import Vizier
import warnings; warnings.filterwarnings('ignore')
for kw in ['SOAR speckle', 'speckle interferometry SOAR', 'SOAR TESS', 'Speckle TESS', 'speckle Gemini TESS', 'Tokovinin speckle']:
    try:
        cats = Vizier.find_catalogs(kw)
    except Exception as e:
        print(kw, 'ERR', e); continue
    print('==', kw, len(cats))
    for k, v in cats.items():
        print('  ', k, '|', v.description[:110])
