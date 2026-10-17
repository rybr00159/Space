"""f06: one TRICERATOPS run (venv_tri) for TIC 214311908 at P = 5.6915087 d.
Usage: f06_tri_run.py <lc_tag> <seed> <N> <cc: none|cons|opt> <nearby: keep|loc|drop> <aperture radius px> <err: emp|prop> [out_tag]
  nearby: keep = all TIC stars that could host the dip (TIC duplicate of pair A removed); loc = only the stars inside the
  pixel-localisation 99.7% region (the 3.1\" and 5.5\" Gaia stars); drop = on-target scenarios only.
* Light curve: data/tri/lc_<lc_tag>.npz from f05 (deblended, detrended, folded, binned). exptime = 200 s (tags *200)
  or 600 s (tags *600); nsamples = 7 (200 s) or 11 (600 s).
* flux_err_0: 'emp' = rms of the out-of-transit bins (|dt| > 0.1 d); 'prop' = median propagated bin error.
* Target overrides: M = 0.370, R = 0.381, Teff = 3350 K, plx = 15.549 mas (star/star_params.json).
* Pixel data: the local TESScut cubes of phase 1 (/home/user/work/tesscut5/214311906/, 21x21 px, sectors 12-104) are
  served to TRICERATOPS through a stand-in for lightkurve.search_tesscut (no download).
* Apertures: per sector, all pixels whose centres lie within <radius> px of the M3V's pixel position (phot/ fitted the
  M3V PRF over pixels within 3 px of the M3V, p05_phot.py RFIT = 3). TRICERATOPS models every star as a Gaussian
  PSF (sigma 0.75 px) and uses TIC Tmags; the Mira (TIC Tmag 11.26) is therefore treated as a normal star.
* tdepth = 0.0157 is the on-target (deblended) depth; dilution_corrected = True.
* Background population: TRICERATOPS 'gaia' option (Gaia DR3 within 0.1 deg^2; file written to data/tri/).
Writes out/tri/<out_tag>.json."""
import sys, os, json, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
from astropy.io import fits
FPP = 'tic214311908/analysis/fpp'
CUBES = '/home/user/work/tesscut5/214311906'
lc_tag, seed, N, ccopt, nearby, aprad, errmode = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], float(sys.argv[6]), sys.argv[7]
out_tag = sys.argv[8] if len(sys.argv) > 8 else f'{lc_tag}_cc{ccopt}_{nearby}_ap{aprad:g}_{errmode}_N{N}_s{seed}'
os.makedirs(f'{FPP}/out/tri', exist_ok=True)
if os.path.exists(f'{FPP}/out/tri/{out_tag}.json') or os.path.exists(f'{FPP}/out/tri/{out_tag}.running'):
    print('exists or running, skipped:', out_tag); sys.exit(0)
open(f'{FPP}/out/tri/{out_tag}.running', 'w').close()
os.chdir(f'{FPP}/data/tri')   # TRICERATOPS writes its Gaia background file into the cwd
np.random.seed(seed)

import triceratops.triceratops as tr
class _Cut:
    def __init__(self, path): self.hdu = fits.open(path, memmap=True)
class _Search:
    def __init__(self, sector): self.sector = sector
    def download_all(self, cutout_size=None):
        f = [x for x in os.listdir(CUBES) if x.startswith(f'tess-s{self.sector:04d}-')]
        return [_Cut(f'{CUBES}/{f[0]}')]
def _search_tesscut(target=None, sector=None, **kw):
    return _Search(int(sector))
tr.lightkurve.search_tesscut = _search_tesscut
# Optional (env TRI_EBCIRC=1): circular orbits for all eclipsing-binary draws (planets keep TRICERATOPS' eccentricity
# prior). Physical motivation: TRICERATOPS' x2P likelihood fits only one eclipse shape and ignores the second eclipse's
# phase (shifted by ~(2/pi) e cos w) and depth (changed by e sin w); our odd/even depths and times agree, and M-dwarf
# binaries at P = 5.7-11.4 d are tidally circularised.
EBCIRC = os.environ.get('TRI_EBCIRC') == '1'
if EBCIRC:
    import triceratops.marginal_likelihoods as _ml
    _orig_sample_ecc = _ml.sample_ecc
    def _sample_ecc(x, planet, P_orb):
        return _orig_sample_ecc(x, planet, P_orb) if planet else np.zeros_like(x)
    _ml.sample_ecc = _sample_ecc

SECTORS = [12, 39, 66, 93, 103, 104]
z = np.load(f'{FPP}/data/tri/lc_{lc_tag}.npz')
tb, fb, eb = z['t'], z['f'], z['e']
oot = np.abs(tb) > 0.10
sig = float(np.std(fb[oot])) if errmode == 'emp' else float(np.median(eb))
exptime = (600. if lc_tag.endswith('600') else 200.)/86400.
nsamp = 11 if lc_tag.endswith('600') else 7
P = 5.6915087
t_start = time.time()
target = tr.target(ID=214311908, sectors=np.array(SECTORS), search_radius=10)
row0 = target.stars.iloc[0]
assert str(row0['ID']) == '214311908'
target.stars.loc[0, 'mass'] = 0.370; target.stars.loc[0, 'rad'] = 0.381
target.stars.loc[0, 'Teff'] = 3350.; target.stars.loc[0, 'plx'] = 15.549
aps = []
for k in range(len(target.pix_coords)):
    x0, y0 = target.pix_coords[k][0]
    xs = np.arange(np.floor(x0 - aprad) - 1, np.ceil(x0 + aprad) + 2)
    ys = np.arange(np.floor(y0 - aprad) - 1, np.ceil(y0 + aprad) + 2)
    XX, YY = np.meshgrid(xs, ys)
    m = np.hypot(XX - x0, YY - y0) <= aprad
    aps.append(np.array([XX[m], YY[m]]).T)
target.calc_depths(tdepth=0.0157, all_ap_pixels=aps)
st = target.stars
nfp = st[st.tdepth > 0][['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth', 'mass', 'rad', 'Teff']]
print(nfp.to_string(), flush=True)
# TIC 1515709367 is a duplicate of TIC 214311910 (pair A; same Gaia DR3 source 5925373272850782208): always removed.
target.remove_star([1515709367])
if nearby == 'drop':      # on-target scenarios only
    target.remove_nearby_stars()
elif nearby == 'loc':     # keep only catalogued stars inside the pixel-localisation 99.7% region (referee_loc, 8.1"):
    # TIC 1515709327 = Gaia DR3 5925373272827260544 (3.1" star), TIC 1515709330 = Gaia DR3 5925373272827262208 (5.5" star)
    keep = ['214311908', '1515709327', '1515709330']
    target.remove_star([i for i in target.stars.ID.astype(str).values if i not in keep])
kw = {}
if ccopt != 'none':
    kw = dict(contrast_curve_file=f'{FPP}/data/tri/cc_{ccopt}.txt', filt='TESS')
target.calc_probs(time=tb, flux_0=fb, flux_err_0=sig, P_orb=P, parallel=True, N=N, verbose=0,
                  exptime=exptime, nsamples=nsamp, **kw)
df = target.probs.copy(); df['lnZ'] = target.lnZ
res = dict(FPP=float(target.FPP), NFPP=float(target.NFPP), FPP_degenerate=bool(getattr(target, 'FPP_degenerate', False)),
           lc_tag=lc_tag, seed=seed, N=N, cc=ccopt, nearby=nearby, aperture_radius_px=aprad, flux_err_0=sig, errmode=errmode,
           exptime_s=exptime*86400, nsamples=nsamp, eb_circular=EBCIRC, tdepth=0.0157, n_aperture_pixels=[int(len(a)) for a in aps],
           nfp_candidates=nfp.astype({'ID': str}).to_dict('records'),
           scenarios=df.astype({'ID': str}).to_dict('records'), runtime_s=time.time() - t_start)
json.dump(res, open(f'{FPP}/out/tri/{out_tag}.json', 'w'), indent=1, default=str)
os.remove(f'{FPP}/out/tri/{out_tag}.running')
print(df.sort_values('prob', ascending=False).head(10).to_string())
print('FPP', target.FPP, 'NFPP', target.NFPP, 'runtime', time.time() - t_start)
