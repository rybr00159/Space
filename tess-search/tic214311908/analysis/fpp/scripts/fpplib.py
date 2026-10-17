"""fpplib: shared constants, paths and data loaders for the false-positive-probability (FPP) analysis of the
5.69-d signal on TIC 214311908. Read-only use of the phase-1/phase-2 folders (phot, referee_phot, loc, referee_loc,
star, referee_star, ctx, referee_ctx, shape, search). Nothing is written outside ../fpp/."""
import numpy as np, json, os
B = 'tic214311908/analysis'
F = f'{B}/fpp'
G = 6.674e-8; RSUN = 6.957e10; MSUN = 1.989e33; AU = 1.496e13; DAY = 86400.0; REARTH = 6.371e8
# ephemeris: mean of the INV/REF free fits (shape/out/final_params.json); P from phase-1 transits.json
T0, P = 3840.9155, 5.6915087
DIST_PC = 64.31
# M3V single-star parameters (star/star_params.json; referee_star r06)
M3V = dict(M=0.370, R=0.381, rho=9.42, Teff=3350., MKs=6.4406, Tmag=13.0135, F_prf=979.0)
TWIN = dict(M=0.264, R=0.289, rho=15.45, Teff=3237.)
FLUX_LOSS = 15.3      # e-/s (15.2-15.4, phase 1)
DEPTH = FLUX_LOSS/979.  # 1.56 % on the M3V (PRF-unit normalisation); fitted depth with LD 15.8-15.9 ppt
SC = 935.0/979.0
LD_M3V = (0.125, 0.533)   # Claret 2017 TESS, 3350 K, logg 4.84 (shape/out/ld_claret17.json)

def jload(path):
    return json.load(open(path if path.startswith('/') else f'{B}/{path}'))

def jdump(obj, path):
    p = path if path.startswith('/') else f'{F}/{path}'
    json.dump(obj, open(p, 'w'), indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))

def inv_clean_windows(halfwin=0.30, sectors=None):
    """INV (phot model B) light curve with the RRab leak removed (search/data/m3v_clean.npz, det_clean, already on the
    979 e-/s scale), flare-flagged cadences removed, per-point sigma = err * white scale * beta (shape/out/noise.json, INV).
    Returns a dict of arrays for |dt| < halfwin."""
    import pandas as pd
    z = np.load(f'{B}/search/data/m3v_clean.npz')
    t, f, e, sec, cad = z['t'], z['det_clean'], z['err'], z['sector'], z['cad']
    d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv', usecols=['time_btjd_tdb', 'quality_ok', 'flare_flag'])
    d = d[d.quality_ok == 1]
    fl = dict(zip(np.round(d.time_btjd_tdb.values, 7), d.flare_flag.values))
    flare = np.array([fl.get(round(x, 7), 0) for x in t]) == 1
    noise = jload('shape/out/noise.json')['INV']
    sig = np.array(e, float)
    for s in np.unique(sec):
        m = sec == s
        sig[m] = e[m]*noise[str(int(s))]['scale']*noise[str(int(s))]['beta']
    ep = np.round((t - T0)/P).astype(int); dt = t - (T0 + ep*P)
    m = (np.abs(dt) < halfwin) & np.isfinite(f) & ~flare
    if sectors is not None: m &= np.isin(sec, sectors)
    return dict(t=t[m], f=f[m], sig=sig[m], sector=sec[m], cad=cad[m].astype(int), epoch=ep[m], dt=dt[m])

def ref_windows(sectors=None):
    """REF (referee local-window photometry) as prepared by shape/s01 (sig includes the fake-window beta)."""
    z = np.load(f'{B}/shape/data/ref_win.npz')
    ep = np.round((z['t'] - T0)/P).astype(int); dt = z['t'] - (T0 + ep*P)
    m = np.isfinite(z['f'])
    if sectors is not None: m &= np.isin(z['sector'], sectors)
    return dict(t=z['t'][m], f=z['f'][m], sig=z['sig'][m], sector=z['sector'][m], cad=z['cad'][m].astype(int),
                epoch=ep[m], dt=dt[m], win=z['win'][m])
