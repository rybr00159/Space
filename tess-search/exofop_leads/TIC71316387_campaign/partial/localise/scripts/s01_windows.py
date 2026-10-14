"""Per-eclipse difference images for real deep (phase 0) and shallow (phase 0.5) eclipses and for fake
windows at phases 0.16-0.34 and 0.66-0.84 (never overlapping real eclipse cadences).
For every window: D (baseline - in-eclipse), sigma(D), B (baseline = direct image at in-eclipse times),
sbar = mean of the exposure-integrated V-shaped eclipse profile over the in-eclipse cadences.
Optional per-pixel pointing regressors: LOC_POINTING=1 (from s01b_pointing.py).
Output: results/win_s{sec}{tag}.npz"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loclib import *

TAG = os.environ.get('LOC_TAG', '')
FAKE_SHIFTS = np.r_[np.arange(0.16, 0.341, 0.02), np.arange(0.66, 0.841, 0.02)]


def run(sec):
    S = Sector(sec)
    reg = None
    if os.environ.get('LOC_POINTING', '0') == '1':
        pz = np.load(f'{L}/results/pointing_s{sec}.npz')
        reg = np.column_stack([pz['dx'], pz['dy']])
    out = dict(D=[], sig=[], B=[], sbar=[], kind=[], tc=[], epoch=[], shift=[], n_in=[], n_l=[], n_r=[])
    all_real = []
    for kind, ph in [('deep', 0.0), ('shallow', 0.5)]:
        tcs, ns = S.eclipse_epochs(phase=ph)
        all_real += list(tcs)
        for tc, n in zip(tcs, ns):
            w = S.window(tc)
            if not w['ok']:
                continue
            D, sig, B, _ = difference_image(S.t, S.flux, w, regressors=reg)
            sh = eclipse_shape(S.t[w['inn']], tc, S.exptime).mean()
            for k, v in [('D', D), ('sig', sig), ('B', B), ('sbar', sh), ('kind', kind), ('tc', tc), ('epoch', n),
                         ('shift', 0.0), ('n_in', w['n_in']), ('n_l', w['n_l']), ('n_r', w['n_r'])]:
                out[k].append(v)
    all_real = np.array(all_real)
    # fake windows: around every deep-eclipse epoch (phase 0) shifted by FAKE_SHIFTS
    tcs, ns = S.eclipse_epochs(phase=0.0)
    for tc, n in zip(tcs, ns):
        for f in FAKE_SHIFTS:
            tf = tc + f * P
            if np.min(np.abs(tf - all_real)) < OOT_OUT + 0.5 * T14 + 0.002:
                continue
            w = S.window(tf)
            if not w['ok']:
                continue
            D, sig, B, _ = difference_image(S.t, S.flux, w, regressors=reg)
            sh = eclipse_shape(S.t[w['inn']], tf, S.exptime).mean()
            for k, v in [('D', D), ('sig', sig), ('B', B), ('sbar', sh), ('kind', 'fake'), ('tc', tf), ('epoch', n),
                         ('shift', f), ('n_in', w['n_in']), ('n_l', w['n_l']), ('n_r', w['n_r'])]:
                out[k].append(v)
    arr = {k: np.array(v) for k, v in out.items()}
    kinds = arr['kind']
    print(f'S{sec}: cam{S.camera} ccd{S.ccd} exptime {S.exptime*1440:.1f} min, good {S.good.sum()}, '
          f'centroid-cut {S.nbad_centroid}; windows deep {np.sum(kinds=="deep")}, shallow {np.sum(kinds=="shallow")}, '
          f'fake {np.sum(kinds=="fake")}', flush=True)
    np.savez_compressed(f'{L}/results/win_s{sec}{TAG}.npz', **arr, jac=S.jac, xt_wcs=S.xt_wcs, yt_wcs=S.yt_wcs,
                        camera=S.camera, ccd=S.ccd, exptime=S.exptime, epoch_year=S.epoch, nbad_centroid=S.nbad_centroid)


if __name__ == '__main__':
    secs = [int(a) for a in sys.argv[1:]] or SECTORS
    for s in secs:
        run(s)
