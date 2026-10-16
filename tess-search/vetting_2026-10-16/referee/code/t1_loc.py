"""TIC 90086266: independent difference-image localisation.
Offset of the eclipse source relative to the target, measured as
(PRF position in difference image) - (PRF position in direct image), with
bootstrap errors, injection tests on the target and on each Gaia neighbour.
"""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 90086266
T14 = 3.0 / 24
P = 2.339474
EV = json.load(open(f'{BASE}/results/t1_events.json'))
rng = np.random.default_rng(42)
RES = {}

for sec in (13, 104):
    cut = load_cut(TIC, sec)
    prf = PRFModel(cut)
    tcs = np.array([e['tc'] for e in EV[f'S{sec}_CUT:r2.0']['events']])
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    tgt = g.iloc[0]
    xg, yg = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    evs, direct = event_diff_images(cut, tcs, T14)
    D, E = combine(evs)
    # direct-image error: use robust scatter of OOT frames / sqrt(n) plus Poisson-ish floor
    derr = np.sqrt(np.abs(direct) * 0.0 + (0.002 * np.abs(direct)) ** 2 + 5.0 ** 2)
    fd = fit_prf(direct, derr, prf, xg, yg, box=3)
    # direct fit including the Gaia neighbours as fixed extra sources (relative)
    near = g[(g.sep_t > 0.01) & (g.sep_t < 60) & (g.tmag < 19)]
    extra = []
    for _, s in near.iterrows():
        xs, ys = cut['wcs'].world_to_pixel_values(s['ra_ep'], s['dec_ep'])
        extra.append((xs - xg, ys - yg, 10 ** (-0.4 * (s['tmag'] - tgt['tmag']))))
    fdn = fit_prf(direct, derr, prf, xg, yg, box=3, extra_sources=extra)
    fdiff = fit_prf(D, E, prf, xg, yg, box=3)
    dE, dN = pix_to_EN(cut, fdiff['x'], fdiff['y'], fd['x'], fd['y'])
    dEn, dNn = pix_to_EN(cut, fdiff['x'], fdiff['y'], fdn['x'], fdn['y'])
    dEg, dNg = pix_to_EN(cut, fdiff['x'], fdiff['y'], xg, yg)
    # bootstrap over events
    bs = []
    for b in range(200):
        idx = rng.integers(0, len(evs), len(evs))
        Db, Eb = combine(evs, idx)
        fb = fit_prf(Db, Eb, prf, fdiff['x'], fdiff['y'], box=3)
        bs.append(pix_to_EN(cut, fb['x'], fb['y'], fd['x'], fd['y']))
    bs = np.array(bs)
    # injection on the target at phase 0.5 (no eclipse there), same depth
    depth = fdiff['A'] / fd['A']
    tmpl_t = prf(fd['x'], fd['y']) * fd['A']
    tcs_inj = tcs + P / 2.
    tcs_inj = tcs_inj[(tcs_inj > cut['time'].min() + 0.3) & (tcs_inj < cut['time'].max() - 0.3)]
    evs_i, direct_i = event_diff_images(cut, tcs_inj, T14, inject=(tmpl_t, depth * 2.0))  # V peak = 2x mean in-window depth approx
    Di, Ei = combine(evs_i)
    fi = fit_prf(Di, Ei, prf, xg, yg, box=3)
    fdi = fit_prf(direct_i, derr, prf, xg, yg, box=3)
    inj_t = pix_to_EN(cut, fi['x'], fi['y'], fdi['x'], fdi['y'])
    # injections on neighbours (each with the same absolute flux deficit as observed)
    inj_n = []
    for _, s in g[(g.sep_t > 0.01) & (g.sep_t < 12)].iterrows():
        xs, ys = cut['wcs'].world_to_pixel_values(s['ra_ep'], s['dec_ep'])
        tmpl = prf(xs, ys) * fd['A']
        evs_n, direct_n = event_diff_images(cut, tcs_inj, T14, inject=(tmpl, depth * 2.0))
        Dn, En = combine(evs_n)
        fn = fit_prf(Dn, En, prf, xg, yg, box=3)
        fdn2 = fit_prf(direct_n, derr, prf, xg, yg, box=3)
        rec = pix_to_EN(cut, fn['x'], fn['y'], fdn2['x'], fdn2['y'])
        fluxfrac = 10 ** (-0.4 * (s['tmag'] - tgt['tmag']))
        inj_n.append(dict(source_id=int(s['source_id']), sep=s['sep_t'], dE_gaia=s['dE'], dN_gaia=s['dN'], G=s['phot_g_mean_mag'],
                          Tmag=s['tmag'], flux_ratio_T=fluxfrac, max_depth_if_100pct=fluxfrac / (1 + fluxfrac),
                          recovered_dE=rec[0], recovered_dN=rec[1]))
    # per-pixel depth map near the target (deficit / direct)
    yc, xc = int(round(fd['y'])), int(round(fd['x']))
    dmap = (D / direct)[yc - 2:yc + 3, xc - 2:xc + 3]
    emap = (E / np.abs(direct))[yc - 2:yc + 3, xc - 2:xc + 3]
    RES[f'S{sec}'] = dict(n_events=len(evs), target_pix_gaia=[xg, yg], direct_fit=[fd['x'], fd['y']], direct_fit_with_nbrs=[fdn['x'], fdn['y']],
                          diff_fit=[fdiff['x'], fdiff['y']], diff_chi2_dof=[fdiff['chi2'], fdiff['dof']],
                          depth_from_prf=depth,
                          offset_vs_direct_arcsec=[dE, dN], offset_vs_direct_with_nbrs=[dEn, dNn], offset_vs_gaia_wcs=[dEg, dNg],
                          bootstrap_std_arcsec=bs.std(axis=0).tolist(), bootstrap_mean=bs.mean(axis=0).tolist(),
                          inj_target_recovered_offset=list(inj_t), injections_neighbours=inj_n,
                          depth_map_5x5=dmap.tolist(), depth_map_err_5x5=emap.tolist())
    print(f'S{sec}: nev={len(evs)} diff-direct offset dE={dE:+.2f}" dN={dN:+.2f}" (bootstrap sd {bs.std(axis=0)[0]:.2f}, {bs.std(axis=0)[1]:.2f}); '
          f'with nbrs in direct fit: {dEn:+.2f},{dNn:+.2f}; vs Gaia WCS {dEg:+.2f},{dNg:+.2f}; inj-on-target recovers {inj_t[0]:+.2f},{inj_t[1]:+.2f}; depth(PRF)={depth:.4f}')
    for r in inj_n:
        print(f"   nbr {r['source_id']} sep={r['sep']:.2f}\" G={r['G']:.2f} T={r['Tmag']:.2f} fluxratio={r['flux_ratio_T']:.4f} -> max depth {r['max_depth_if_100pct']*100:.2f}%  "
              f"Gaia (dE,dN)=({r['dE_gaia']:+.2f},{r['dN_gaia']:+.2f}) injected-recovered=({r['recovered_dE']:+.2f},{r['recovered_dN']:+.2f})")
    print('   depth map (5x5 around target, %):')
    for row in dmap[::-1]:
        print('    ' + ' '.join('%6.2f' % (v * 100) for v in row))
save_json(RES, f'{BASE}/results/t1_loc.json')
