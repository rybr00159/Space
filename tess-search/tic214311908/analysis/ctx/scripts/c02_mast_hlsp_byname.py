"""Query MAST by target name (TIC and Gaia DR3 IDs) and by provenance for HLSPs that the coordinate query may miss (TGLC, TESS-SPOC, CDIPS, PATHOS, DIAmante, etc.)."""
import warnings; warnings.filterwarnings("ignore")
from astroquery.mast import Observations
import astropy.units as u
from astropy.coordinates import SkyCoord
names = {"214311906": "5925373268528684288", "214311908": "5925373272850784768", "214311910": "5925373272850782208"}
for tic, gaia in names.items():
    for nm in (tic, "TIC " + tic, "TIC-" + tic, gaia, "Gaia DR3 " + gaia, "gaiaid-" + gaia):
        try:
            r = Observations.query_criteria(target_name=nm)
            if len(r): print(nm, len(r), sorted(set(r["provenance_name"])))
        except Exception as e:
            print(nm, "ERR", e)
c = SkyCoord(260.977650, -51.323690, unit="deg")
for prov in ("TGLC", "TESS-SPOC", "SPOC", "CDIPS", "PATHOS", "DIAMANTE", "TASOC", "T16", "GSFC-ELEANOR-LITE", "TARS", "TEQUILA", "QLP"):
    try:
        r = Observations.query_criteria(coordinates=c, radius=0.02*u.deg, provenance_name=prov)
        print(prov, len(r), sorted(set(r["target_name"]))[:20])
    except Exception as e:
        print(prov, "ERR", e)
