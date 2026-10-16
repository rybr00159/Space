"""TIC 214311906 / 55829034: difference-image PRF localisation per sector,
combined offset, and implied intrinsic depth for each candidate neighbour.
Usage: python -I t34_loc.py <TIC>"""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from reflib import *
from diffimg import *

TIC = int(sys.argv[1])
CFG = {214311906: dict(P=5.69149, T0=3840.913, T14=1.55 / 24, secs=(12, 39, 66, 93, 103, 104), rmax=25),
       55829034: dict(P=0.616146, T0=3837.8667, T14=2.0 / 24, secs=(13, 66, 93, 104), rmax=20)}[TIC]
P, T0, T14 = CFG['P'], CFG['T0'], CFG['T14']
ZP = 20.44  # TESS zero point (e-/s)
rng = np.random.default_rng(7)
OUT = dict(sectors={})
offs = []
for sec in CFG['secs']:
    cut = load_cut(TIC, sec)
    t = cut['time']
    # robust per-cadence sky: 20th percentile of the cutout
    sky = np.percentile(cut['flux'].reshape(len(t), -1), 20, axis=1)
    cut['flux'] = cut['flux'] - sky[:, None, None]
    n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
    tcs = T0 + np.arange(n0, n1 + 1) * P
    tcs = tcs[(tcs > t.min() + 0.1) & (tcs < t.max() - 0.1)]
    evs, direct = event_diff_images(cut, tcs, T14, fin=0.35, gap=0.65, wout=max(0.12, 1.5 * T14))
    if len(evs) < 2:
        continue
    D, E = combine(evs)
    g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
    tgt = g.iloc[0]
    xt, yt = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    prf = PRFModel(cut)
    derr = np.sqrt((0.003 * np.abs(direct)) ** 2 + 10 ** 2)
    ft = fit_prf(direct, derr, prf, xt, yt, box=1)
    fd = fit_prf(D, E, prf, xt, yt, box=2)
    o = pix_to_EN(cut, fd['x'], fd['y'], ft['x'], ft['y'])
    ow = pix_to_EN(cut, fd['x'], fd['y'], xt, yt)
    # bootstrap over events
    bs = []
    for b in range(150):
        idx = rng.integers(0, len(evs), len(evs))
        Db, Eb = combine(evs, idx)
        fb = fit_prf(Db, Eb, prf, fd['x'], fd['y'], box=2)
        bs.append(pix_to_EN(cut, fb['x'], fb['y'], ft['x'], ft['y']))
    bs = np.array(bs)
    sd = np.std(bs, axis=0)
    # robust sd (bootstrap can have outliers when S/N is low)
    sdr = 1.4826 * np.median(np.abs(bs - np.median(bs, axis=0)), axis=0)
    dF = fd['A']
    dFerr = np.sqrt(fd['cov'][2, 2])
    rec = dict(n_events=len(evs), offset_vs_direct=list(o), offset_vs_wcs=list(ow), boot_sd=sd.tolist(), boot_rsd=sdr.tolist(),
               deltaF_eps=dF, deltaF_err=dFerr, target_flux_prf=ft['A'], depth_vs_target_prf_ppm=dF / ft['A'] * 1e6,
               wcs_minus_direct_arcsec=list(pix_to_EN(cut, xt, yt, ft['x'], ft['y'])))
    OUT['sectors'][sec] = rec
    offs.append((sec, o[0], o[1], sdr[0], sdr[1], dF, dFerr, ft['A']))
    print(f"S{sec:3d}: n={len(evs)} diff-image source at ({o[0]:+5.1f},{o[1]:+5.1f})\" from target direct-image fit (boot robust sd {sdr[0]:.1f},{sdr[1]:.1f}); "
          f"vs WCS ({ow[0]:+5.1f},{ow[1]:+5.1f}); dF={dF:.1f}+-{dFerr:.1f} e/s; target PRF flux {ft['A']:.0f} e/s -> depth vs target {dF/ft['A']*1e6:.0f} ppm")
offs = np.array(offs)
w = 1 / np.maximum(offs[:, 3:5], 0.5) ** 2
mE = np.sum(offs[:, 1] * w[:, 0]) / np.sum(w[:, 0]); sE = 1 / np.sqrt(np.sum(w[:, 0]))
mN = np.sum(offs[:, 2] * w[:, 1]) / np.sum(w[:, 1]); sN = 1 / np.sqrt(np.sum(w[:, 1]))
chiE = np.sum(((offs[:, 1] - mE) / np.maximum(offs[:, 3], 0.5)) ** 2); chiN = np.sum(((offs[:, 2] - mN) / np.maximum(offs[:, 4], 0.5)) ** 2)
print(f'combined offset from target: dE={mE:+.1f}+-{sE:.1f}" dN={mN:+.1f}+-{sN:.1f}" (chi2 E {chiE:.1f}, N {chiN:.1f} for {len(offs)-1} dof)')
OUT['combined'] = dict(dE=mE, dE_err=sE, dN=mN, dN_err=sN, chi2E=chiE, chi2N=chiN, n=len(offs))
# deltaF weighted mean (e-/s)
wf = 1 / offs[:, 6] ** 2
dFm = np.sum(offs[:, 5] * wf) / np.sum(wf); dFs = 1 / np.sqrt(np.sum(wf))
OUT['deltaF_mean'] = [dFm, dFs]
# candidate table
g = gaia_table(TIC, btjd_to_jyear(4000.))
cand = g[(g.sep_t < CFG['rmax'])].copy()
cand['flux_eps'] = 10 ** (-0.4 * (cand['tmag'] - ZP))
cand['intrinsic_depth'] = dFm / cand['flux_eps']
cand['chi2_pos'] = ((cand['dE'] - mE) / np.hypot(sE, 0.5)) ** 2 + ((cand['dN'] - mN) / np.hypot(sN, 0.5)) ** 2
cand = cand.sort_values('chi2_pos')
print(f'mean deltaF = {dFm:.1f} +- {dFs:.1f} e-/s; candidates (sep<{CFG["rmax"]}"), sorted by positional chi2:')
rows = []
for _, c in cand.head(15).iterrows():
    print(f"  {int(c['source_id'])} sep={c['sep_t']:5.1f}\" (dE,dN)=({c['dE']:+5.1f},{c['dN']:+5.1f}) G={c['phot_g_mean_mag']:.2f} BP-RP={c['bp_rp']:.2f} T~{c['tmag']:.2f} "
          f"plx={c['parallax']:.2f} ruwe={c['ruwe']:.1f} -> intrinsic depth {c['intrinsic_depth']*100:6.1f}%  pos chi2={c['chi2_pos']:.1f}")
    rows.append(dict(source_id=int(c['source_id']), sep=c['sep_t'], dE=c['dE'], dN=c['dN'], G=c['phot_g_mean_mag'], bp_rp=c['bp_rp'], tmag=c['tmag'],
                     parallax=c['parallax'], ruwe=c['ruwe'], intrinsic_depth=c['intrinsic_depth'], chi2_pos=c['chi2_pos']))
OUT['candidates'] = rows
save_json(OUT, f'{BASE}/results/t34_loc_{TIC}.json')
