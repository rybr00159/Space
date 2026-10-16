"""QLP time stamps vs TESScut (barycentred) for TIC 4206066 in S98: nearest-cadence residual test.
If QLP TIME = TESScut TIME - TIMECORR, the residual to the nearest TESScut cadence is about -(TIMECORR mod cadence)."""
import sys, glob, io, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
import numpy as np, requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
import fetch_qlp
for tic, ra, dec, sectors in ((4206066, 80.8598409616627, -1.32727827489794, (98,)),):
    for s in sectors:
        d = f"/home/user/work/tesscut_small/{tic}"
        if not glob.glob(f"{d}/*s{s:04d}*astrocut.fits"):
            Tesscut.download_cutouts(coordinates=SkyCoord(ra, dec, unit="deg"), size=3, sector=s, path=d)
        f = glob.glob(f"{d}/*s{s:04d}*astrocut.fits")[0]
        with fits.open(f) as h:
            tc = np.array(h[1].data["TIME"], float); tcorr = np.array(h[1].data["TIMECORR"], float)
        r = requests.get(fetch_qlp.url(tic, s), timeout=120)
        with fits.open(io.BytesIO(r.content)) as h:
            tq = np.array(h[1].data["TIME"], float)
        tq = tq[np.isfinite(tq)]; ok = np.isfinite(tc); tc, tcorr = tc[ok], tcorr[ok]
        j = np.clip(np.searchsorted(tc, tq), 1, len(tc) - 1); j = np.where(np.abs(tc[j - 1] - tq) < np.abs(tc[j] - tq), j - 1, j)
        res = (tq - tc[j]) * 86400; cad = np.median(np.diff(tc)) * 86400; tcm = np.median(tcorr) * 86400
        exp_unc = -(((tcm + cad / 2) % cad) - cad / 2)
        print(f"TIC {tic} S{s}: cadence {cad:.0f} s; TIMECORR median {tcm:.1f} s; QLP - nearest TESScut = {np.median(res):+.1f} s "
              f"(5-95% {np.percentile(res,5):+.1f}..{np.percentile(res,95):+.1f}); expected if QLP uncorrected: {exp_unc:+.1f} s, if corrected: ~0 s")
