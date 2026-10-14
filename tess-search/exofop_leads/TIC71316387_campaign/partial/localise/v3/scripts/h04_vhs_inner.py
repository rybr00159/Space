"""VHS Ks: inner-separation sensitivity and photometric bias of the companion fit (h01 machinery).
For the three Ks frames (2 pawprints + tile): scanner with rmin = 0.3"; real maximum recomputed; inject companions at
0.4-1.2" (8 random PAs), dKs = 2.0-4.0; 'recovered' = best-fit position within 0.35" of the injected one AND
dchi2 > real maximum.  Also: recovered dKs for injections at the real feature's position region (1.1", PA ~350) to
measure the flux bias of the fitted companion next to the masked core.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from h01_vhs import load, build_epsf, psf, fit_single, Scanner, stamp, local_sky, H, GAIN, RN, SATLIM, g
from v3lib import V3, to_json
import warnings; warnings.filterwarnings('ignore')

out = []
for f in [V3 + '/imaging/vhs_ADP.2014-12-22T14:52:14.147.fits.fz', V3 + '/imaging/vhs_ADP.2014-12-22T14:59:09.927.fits.fz',
          V3 + '/imaging/vhs_ADP.2014-12-22T14:57:07.660.fits.fz']:
    I = load(f); d = I['d']; x, y = I['x'], I['y']
    it = int(np.where(g.source_id.values == 6571343840265900544)[0][0])
    pix = float(np.sqrt(abs(np.linalg.det(I['J'])))); Jn = I['J'] / pix
    cand = []
    for i in range(len(g)):
        if i == it or g.phot_g_mean_mag[i] > 17.3:
            continue
        if not (H + 45 < x[i] < d.shape[1] - H - 45 and H + 45 < y[i] < d.shape[0] - H - 45):
            continue
        dd = np.hypot((g.ra - g.ra[i]) * np.cos(np.radians(g.dec[i])), g.dec - g.dec[i]) * 3600
        if np.any((dd > 0.01) & (dd < 6) & (g.phot_g_mean_mag - g.phot_g_mean_mag[i] < 5)):
            continue
        st_, _, _ = stamp(d, x[i], y[i])
        if st_.max() > 20000:
            continue
        cand.append(i)
    E = build_epsf(I, cand)
    raw, _, _ = stamp(d, x[it], y[it]); st = raw - local_sky(d, x[it], y[it])
    sig = np.sqrt(np.clip(raw, 1, None) / GAIN + RN ** 2)
    yy, xx = np.mgrid[-H:H + 1, -H:H + 1]
    mask = (raw < SATLIM) & (np.hypot(xx, yy) <= H)
    p0, _, _ = fit_single(st, sig, mask, E)
    SC = Scanner(E, p0, pix, rmin=0.3, rmax=3.0, step=0.1)
    b, chi0, c0 = SC.run(st, sig, mask)
    F = c0[0]
    rng = np.random.default_rng(21)
    rec = []
    for sep in (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2):
        for dm in (2.0, 2.5, 3.0, 3.5, 4.0):
            ok = []; dmr = []
            for pa in np.arange(0, 360, 45.0) + rng.uniform(0, 45):
                dpx = np.linalg.solve(Jn, sep * np.array([np.sin(np.radians(pa)), np.cos(np.radians(pa))]))
                comp = F * 10 ** (-0.4 * dm) * psf(E, p0[0] + dpx[0] / pix, p0[1] + dpx[1] / pix)
                raw2 = raw + comp; mk2 = (raw2 < SATLIM) & (np.hypot(xx, yy) <= H)
                sg2 = np.sqrt(np.clip(raw2, 1, None) / GAIN + RN ** 2)
                bb, _, _ = SC.run(st + comp, sg2, mk2)
                hit = np.hypot(bb[1][0] - dpx[0], bb[1][1] - dpx[1]) < 0.35
                ok.append(bool(hit and bb[0] > b[0]))
                if hit and bb[2] > 0:
                    dmr.append(-2.5 * np.log10(bb[2]))
            rec.append(dict(sep=sep, dm=dm, frac=float(np.mean(ok)), dm_rec_med=float(np.median(dmr)) if dmr else None))
    so = Jn @ np.array(b[1])
    r = dict(file=os.path.basename(f), real_best=dict(dchi2=b[0], sep=float(np.hypot(*b[1])), pa=float(np.degrees(np.arctan2(*so)) % 360),
                                                     dm=float(-2.5 * np.log10(b[2]))), inj=rec, mask_r=float(np.sqrt((~mask & (np.hypot(xx, yy) <= H)).sum() / np.pi) * pix))
    out.append(r)
    print(f"{r['file']}: mask r {r['mask_r']:.2f}\"; real best (rmin 0.3) dchi2 {b[0]:.0f} sep {r['real_best']['sep']:.2f} PA {r['real_best']['pa']:.0f} dKs {r['real_best']['dm']:.2f}", flush=True)
    for sep in (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2):
        print(f"   {sep}\": " + '  '.join(f"dm{q['dm']}:{q['frac']:.2f}(rec {q['dm_rec_med'] if q['dm_rec_med'] is None else round(q['dm_rec_med'], 2)})" for q in rec if q['sep'] == sep), flush=True)
to_json(out, V3 + '/results/h04_vhs_inner.json')
