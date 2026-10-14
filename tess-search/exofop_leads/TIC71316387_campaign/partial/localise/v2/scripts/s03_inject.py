"""Injection-recovery calibration of the difference-image offsets.

For each sector and eclipse type (deep / shallow):
  * draw n_real fake windows (one per eclipse epoch, random phase shift 0.16-0.34 or 0.66-0.84; never overlapping
    real eclipses) and stack them exactly as the real windows are stacked;
  * inject an eclipse of the real fitted amplitude A_real.  Because the difference image is linear in the cube,
    D(cube + injection) = D(cube) + A * sbar * image, so the injection is added at the stacked level (stacked units
    are already normalised by sbar);
  * injected image = EMPIRICAL unit-flux target image E (from the mean direct image of all fake windows, background
    and Gaia neighbours removed), shifted by delta using the PRF gradient: E_delta = E + PRF(x_T + delta) - PRF(x_T);
    delta = 0, and 1", 3", 5" at 8 position angles;
  * variant 'prf': injected image = PRF model itself at x_T + delta (isolates noise-only behaviour);
  * fit with the identical pipeline (direct-image PRF fit for x_T, free PRF + constant for D).
Output: results/inj_s{sec}.npz and a summary json.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastfit import *
from loclib import R, SECTORS, to_json

NT = int(os.environ.get('INJ_NT', 300))
SEPS = [0.0, 1.0, 3.0, 5.0]
PAS = np.arange(0, 360, 45.0)
rng = np.random.default_rng(int(os.environ.get('INJ_SEED', 12345)))


def run(sec):
    t0 = time.time()
    F = Fast(sec); ctx = F.ctx; z = ctx.z; m = F.mask
    kind = ctx.kind
    real = json.load(open(f'{R}/results/real_fast.json'))[str(sec)]
    fk = np.where(kind == 'fake')[0]
    ep = z['epoch'][fk]
    sb = z['sbar']; sg = z['sig']; D = z['D']; B = z['B']
    med = np.array([np.median(s[m] ** 2) for s in sg])
    # empirical target image from all fake windows' direct images
    wb_all = 1.0 / med[fk]
    B_all = (wb_all[:, None, None] * B[fk]).sum(0) / wb_all.sum()
    s_all = np.sqrt(np.median(sg[fk] ** 2, axis=0) / len(fk))
    E, rB0 = F.empirical_image(B_all, s_all)
    P0 = F.prf_img(rB0['dx'], rB0['dy'])
    Ji = np.linalg.inv(F.jac)
    configs = [(0.0, 0.0)] + [(s, pa) for s in SEPS[1:] for pa in PAS]
    out = {}
    for name in ['deep', 'shallow']:
        n_real = int(real[name]['n']); A = real[name]['amp']
        uniq = np.unique(ep)
        res = {mode: np.full((len(configs), NT, 2), np.nan) for mode in ['emp', 'prf']}
        amps = np.full((len(configs), NT), np.nan)
        null_off = np.full((NT, 2), np.nan)
        for t in range(NT):
            eps = rng.choice(uniq, size=min(n_real, len(uniq)), replace=False)
            idx = np.array([rng.choice(fk[ep == e]) for e in eps])
            a = sb[idx]
            wk = a ** 2 / med[idx]
            W = wk.sum()
            Ds = (wk[:, None, None] * D[idx] / a[:, None, None]).sum(0) / W
            ss = np.sqrt((wk[:, None, None] ** 2 * (sg[idx] / a[:, None, None]) ** 2).sum(0)) / W
            wb = 1.0 / med[idx]
            Bs = (wb[:, None, None] * B[idx]).sum(0) / wb.sum()
            rB = F.fit_direct(Bs, ss)
            for ci, (sep, pa) in enumerate(configs):
                dsky = np.array([sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))])
                dpx = Ji @ dsky
                Pd = F.prf_img(rB0['dx'] + dpx[0], rB0['dy'] + dpx[1])
                for mode, img in [('emp', E + Pd - P0), ('prf', Pd)]:
                    rD = F.fit_diff(Ds + A * img, ss)
                    sky, _ = F.offset_sky(rD, rB)
                    res[mode][ci, t] = sky - dsky      # recovered minus injected
                    if mode == 'emp':
                        amps[ci, t] = rD['amp']
        summ = []
        for ci, (sep, pa) in enumerate(configs):
            row = dict(sep=sep, pa=pa)
            for mode in ['emp', 'prf']:
                r = res[mode][ci]
                row[mode] = dict(bias=np.nanmean(r, 0).tolist(), cov=np.cov(r.T).tolist(),
                                 sd=np.nanstd(r, 0).tolist(),
                                 mad_sd=(1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r, 0)), 0)).tolist(),
                                 p95_radius=float(np.nanpercentile(np.hypot(*(r - np.nanmean(r, 0)).T), 95)))
            row['amp_mean'] = float(np.nanmean(amps[ci]))
            summ.append(row)
        out[name] = dict(n_real=n_real, A=A, summary=summ)
        np.savez_compressed(f'{R}/results/inj_s{sec}_{name}.npz', emp=res['emp'], prf=res['prf'], amps=amps,
                            configs=np.array(configs))
        s0 = summ[0]
        print(f"S{sec} {name:7s} n={n_real} A={A:.1f}: at target emp bias ({s0['emp']['bias'][0]:+.2f},{s0['emp']['bias'][1]:+.2f})\" "
              f"sd ({s0['emp']['sd'][0]:.2f},{s0['emp']['sd'][1]:.2f}) | prf bias ({s0['prf']['bias'][0]:+.2f},{s0['prf']['bias'][1]:+.2f}) "
              f"sd ({s0['prf']['sd'][0]:.2f},{s0['prf']['sd'][1]:.2f})  [{time.time()-t0:.0f}s]", flush=True)
        for sep in SEPS[1:]:
            rows = [r for r in summ if r['sep'] == sep]
            b = np.array([r['emp']['bias'] for r in rows]); sd = np.array([r['emp']['sd'] for r in rows])
            print(f"      inj {sep:.0f}\": |bias| max {np.max(np.hypot(*b.T)):.2f}\" mean bias ({b[:,0].mean():+.2f},{b[:,1].mean():+.2f}) sd mean ({sd[:,0].mean():.2f},{sd[:,1].mean():.2f})", flush=True)
    out['direct_fit_all_fakes'] = dict(dx=rB0['dx'], dy=rB0['dy'], F=rB0['F'], c=rB0['c'])
    to_json(out, f'{R}/results/inj_s{sec}.json')


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
