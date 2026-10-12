"""Step 3: injection-recovery calibration on the real cutouts.

Synthetic transits (trapezoid, T14 = 0.1875 d, T12 = 0.045 d, exposure-integrated; amplitude = the
real per-sector difference-image amplitude) are injected as scaled PSF dips on the TARGET or on the
NEIGHBOUR at fake epochs tc + phi*P (phi = 0.25..0.75, windows never overlapping real transit windows).
The identical pipeline (window selection, per-pixel linear baselines, stacking, fits) is run.
Because the difference image is linear in the cube, D(cube + inj) = D(cube) + D(inj) exactly; the
script verifies this identity on a literal injected cube for one window per sector.

Two kinds of simulated data sets per sector and host:
  'literal'  : the real ephemeris shifted by phi (one data set per phi);
  'resampled': for each real transit a randomly chosen fake window (no repeats) -> NSIM data sets.
Also 'null' data sets (no injection) give the pure-noise behaviour.

usage: python s03_inject.py <inj_psf: prf|gauss> <ana_psf: prf|gauss> <real_tag of the analysis model> [rfit] [nsim]
(the injected amplitude is the real amplitude measured with the injection model, diff_<inj_psf>_main.json)
Outputs: inj_<inj>_<ana>.json, inj_<inj>_<ana>.npz
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pipeline import *

inj_psf = sys.argv[1] if len(sys.argv) > 1 else 'prf'
ana_psf = sys.argv[2] if len(sys.argv) > 2 else 'prf'
real_tag = sys.argv[3] if len(sys.argv) > 3 else 'prf_main'
rfit = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
NSIM = int(sys.argv[5]) if len(sys.argv) > 5 else 600
REGV = {'prf': 'local_r4', 'gauss': 'local_r4_ratiofix'}
rng = np.random.default_rng(7)
KEYS = ['offT_E', 'offT_N', 's_par', 'perp_arcsec', 'f_N', 'a_T', 'a_N', 'dchi2_TminusN', 'free_amp']


def run_sector(sec, real, real_inj):
    S = Sector(sec)
    ictx = SectorContext(sec, psfkind=inj_psf, reg_variant=REGV[inj_psf], rfit=rfit, S=S)
    actx = ictx if (ana_psf == inj_psf) else SectorContext(sec, psfkind=ana_psf, reg_variant=REGV[ana_psf], rfit=rfit, S=S)
    rw = actx.real_windows()
    nT = len(rw)
    # amplitude measured on the real data WITH THE INJECTION PSF model (so injected S/N matches reality)
    A_real = real_inj[str(sec)]['stacked']['free_amp']
    shape_mean = np.mean([ictx.in_window_shape_mean(w) for w in rw])
    amp = A_real / shape_mean
    fw = ictx.fake_windows()
    # per-window products
    Dn, Sg, DT, DN, phis, eps = [], [], [], [], [], []
    for w in fw:
        D, sg, _ = actx.diff(w)
        Dn.append(D); Sg.append(sg)
        DT.append(ictx.inject_diff(w, ictx.pT, amp))
        DN.append(ictx.inject_diff(w, ictx.pN, amp))
        phis.append(w['phi']); eps.append(w['epoch'])
    Dn, Sg, DT, DN = map(np.array, (Dn, Sg, DT, DN))
    phis = np.array(phis); eps = np.array(eps)
    # literal-injection identity check on one window
    w = fw[len(fw) // 2]
    idx = np.concatenate([w['lft'], w['inn'], w['rgt']])
    cube = S.flux.copy()
    sub = np.linspace(-0.5, 0.5, 11) * S.exptime
    shp = np.mean([ictx.trapezoid(S.t[idx] + s, w['tc']) for s in sub], axis=0)
    cube[idx] -= amp * shp[:, None, None] * ictx.psf(*ictx.pT)[None]
    Dlit, _, _ = difference_image(S.t, cube, w)
    ident = float(np.max(np.abs(Dlit - (Dn[len(fw) // 2] + DT[len(fw) // 2]))))
    out = dict(n_real=nT, n_fake_windows=len(fw), A_real=A_real, amp_injected=amp, identity_maxabs=ident)
    sims = {}
    for mode in ['literal', 'resampled']:
        if mode == 'literal':
            sets = []
            for phi in np.unique(phis):
                k = np.where(np.isclose(phis, phi))[0]
                if len(k) >= max(2, nT - 2):
                    sets.append(k)
        else:
            sets = [rng.choice(len(fw), size=min(nT, len(fw)), replace=False) for _ in range(NSIM)]
        for host in ['T', 'N', 'null']:
            vals = []
            for k in sets:
                Dk = Dn[k] + (DT[k] if host == 'T' else DN[k] if host == 'N' else 0.0)
                Ds, ss = stack(Dk, Sg[k], mask=actx.mask)
                a = actx.analyse(Ds, ss, free=True, extra=False)
                vals.append([a[q] for q in KEYS])
            sims[f'{mode}_{host}'] = np.array(vals)
    out['summary'] = {}
    for key, v in sims.items():
        out['summary'][key] = {q: dict(mean=float(np.nanmean(v[:, j])), std=float(np.nanstd(v[:, j], ddof=1)),
                                       median=float(np.nanmedian(v[:, j])), n=int(np.isfinite(v[:, j]).sum()))
                               for j, q in enumerate(KEYS)}
    return out, sims


def main():
    real = json.load(open(f'{HERE}/diff_{real_tag}.json'))       # real result for the ANALYSIS model
    real_inj = json.load(open(f'{HERE}/diff_{inj_psf}_main.json'))  # real amplitude for the INJECTION model
    res = {}
    arrays = {}
    for sec in SECTORS:
        t0 = time.time()
        out, sims = run_sector(sec, real, real_inj)
        res[str(sec)] = out
        for k, v in sims.items():
            arrays[f's{sec}_{k}'] = v
        sm = out['summary']
        r = real[str(sec)]['stacked']
        print(f'S{sec} inj={inj_psf} ana={ana_psf} A={out["A_real"]:.1f} e/s  identity={out["identity_maxabs"]:.2e}  '
              f'nfake={out["n_fake_windows"]}', flush=True)
        for key in ['resampled_T', 'resampled_N', 'resampled_null', 'literal_T', 'literal_N']:
            s = sm[key]
            print(f'   {key:15s} n={s["s_par"]["n"]:4d}  s_par={s["s_par"]["mean"]:+.3f}+-{s["s_par"]["std"]:.3f}  '
                  f'f_N={s["f_N"]["median"]:+.3f}+-{s["f_N"]["std"]:.3f}  dchi2={s["dchi2_TminusN"]["mean"]:+.1f}+-{s["dchi2_TminusN"]["std"]:.1f}'
                  f'  offT=({s["offT_E"]["mean"]:+.2f},{s["offT_N"]["mean"]:+.2f})', flush=True)
        print(f'   REAL            s_par={r["s_par"]:+.3f}  f_N={r["f_N"]:+.3f}  dchi2={r["dchi2_TminusN"]:+.1f}  ({time.time()-t0:.0f}s)', flush=True)
    to_json(res, f'{HERE}/inj_{inj_psf}_{ana_psf}.json')
    np.savez_compressed(f'{HERE}/inj_{inj_psf}_{ana_psf}.npz', keys=np.array(KEYS), **arrays)


if __name__ == '__main__':
    main()
