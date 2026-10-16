import os, sys, time, glob, warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
OUT = "/home/user/work/tesscut5"
T = {90086266: (276.27770, -42.64403, [13, 104]),
     61332759: (287.61974, -44.73357, [13, 27, 67, 94, 104]),
     214311906: (260.97765, -51.32369, [12, 39, 66, 93, 103, 104]),
     55829034: (272.03163, -35.69140, [13, 66, 93, 104]),
     61141412: (286.96474, -43.00969, [13, 27, 67, 94, 104])}
for tic, (ra, dec, secs) in T.items():
    d = f"{OUT}/{tic}"; os.makedirs(d, exist_ok=True)
    for s in secs:
        if glob.glob(f"{d}/*s{s:04d}*astrocut.fits"): continue
        for a in range(5):
            try:
                Tesscut.download_cutouts(coordinates=SkyCoord(ra, dec, unit="deg"), size=21, sector=s, path=d)
                print(tic, s, "ok", flush=True); break
            except Exception as e:
                print(tic, s, "retry", a, str(e)[:100], flush=True); time.sleep(10 * 2 ** a)
