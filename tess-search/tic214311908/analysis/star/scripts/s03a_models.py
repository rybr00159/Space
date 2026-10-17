"""Download BT-Settl CIFIST2011 model spectra (SVO theory service) for Teff 2800-4000 K (100 K steps), logg 4.5/5.0/5.5,
[M/H]=0, rebin to a 2 A grid from 0.25 to 30 micron, save as data/models/btsettl_<teff>_<logg>.npy (lam_A, flam).
Also downloads the filter curves used in the SED fit (SVO FPS) into data/models/filters.npz."""
import os, io, time, re, warnings
warnings.filterwarnings('ignore')
import numpy as np, requests
from astropy.io.votable import parse
from astroquery.svo_fps import SvoFps
M = 'tic214311908/analysis/star/data/models'
grid = np.concatenate([np.arange(2500, 3000, 100), np.arange(2900, 4101, 100)])
grid = np.arange(2600, 4101, 100)
lam_out = np.arange(2500., 300000., 2.0)
def fetch(teff, logg):
    fn = f'{M}/btsettl_{teff}_{logg:.1f}.npy'
    if os.path.exists(fn): return
    u = f"http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=bt-settl-cifist&teff_min={teff}&teff_max={teff}&logg_min={logg}&logg_max={logg}&REQUEST=queryData"
    txt = requests.get(u, timeout=120).text
    m = re.search(r'CDATA\[(http[^\]]+fid=\d+)\]', txt)
    if not m: print('no model', teff, logg); return
    r = requests.get(m.group(1).replace('&amp;', '&'), timeout=300)
    vt = parse(io.BytesIO(r.content)).get_first_table().to_table()
    lam = np.array(vt.columns[0], float); fl = np.array(vt.columns[1], float)
    o = np.argsort(lam); lam, fl = lam[o], fl[o]
    # integrate onto 2 A bins (flux-conserving via cumulative trapezoid)
    edges = np.concatenate([[lam_out[0] - 1], lam_out + 1])
    cum = np.concatenate([[0], np.cumsum(0.5 * (fl[1:] + fl[:-1]) * np.diff(lam))])
    ce = np.interp(edges, lam, cum)
    fb = np.diff(ce) / 2.0
    # total bolometric flux of the full model (for the Fbol tail beyond 30 micron)
    np.save(fn, np.vstack([lam_out, fb]).astype(np.float32))
    np.save(fn.replace('.npy', '_fbol.npy'), np.array([np.trapezoid(fl, lam)]))
    print('saved', teff, logg, len(lam), flush=True)
for logg in (5.0, 4.5, 5.5):
    for teff in grid:
        if logg != 5.0 and not (3000 <= teff <= 3700): continue
        for k in range(3):
            try:
                fetch(int(teff), logg); break
            except Exception as e:
                print('retry', teff, logg, e); time.sleep(5)
filters = ['GAIA/GAIA3.G', 'GAIA/GAIA3.Gbp', 'GAIA/GAIA3.Grp', '2MASS/2MASS.J', '2MASS/2MASS.H', '2MASS/2MASS.Ks',
           'SkyMapper/SkyMapper.g', 'SkyMapper/SkyMapper.r', 'SkyMapper/SkyMapper.i', 'SkyMapper/SkyMapper.z', 'SkyMapper/SkyMapper.v',
           'WISE/WISE.W1', 'WISE/WISE.W2', 'WISE/WISE.W3', 'DENIS/DENIS.I', 'TESS/TESS.Red', 'Paranal/VISTA.J', 'Paranal/VISTA.Ks']
out = {}
for f in filters:
    t = SvoFps.get_transmission_data(f)
    out[f.replace('/', '__')] = np.vstack([np.array(t['Wavelength'], float), np.array(t['Transmission'], float)])
np.savez(f'{M}/filters.npz', **out)
info = {}
for fac in ['GAIA', '2MASS', 'SkyMapper', 'WISE', 'DENIS', 'TESS', 'Paranal']:
    try:
        fl = SvoFps.get_filter_list(facility=fac)
        for r in fl:
            if r['filterID'] in filters:
                info[r['filterID']] = dict(zp_jy=float(r['ZeroPoint']), lpivot=float(r['WavelengthPivot']), magsys=str(r['MagSys']), zptype=str(r['ZeroPointType']))
    except Exception as e:
        print(fac, e)
import json; json.dump(info, open(f'{M}/filters_info.json', 'w'), indent=1); print(info)
