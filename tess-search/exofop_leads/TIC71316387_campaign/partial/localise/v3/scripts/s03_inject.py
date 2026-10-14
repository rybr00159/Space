"""v3 injection-recovery calibration of the difference-image offsets (per sector, per eclipse type).

Each trial: choose n_real fake windows (distinct epochs, off-eclipse phases 0.15-0.35 / 0.65-0.85, real out-of-eclipse
frames), stack them exactly as the real windows (same normalisation), fit the direct image, then inject an eclipse of
the real stacked amplitude A at sky offset delta and fit the difference image with the identical pipeline.
Injection is done on the stacked difference image, which is exactly equivalent to injecting
A * vshape(t) * image into every in-eclipse frame (the difference operator is linear and the baseline cadences carry
no injected signal); s03 --check verifies this equivalence for one trial at the cube level.
Injected image: 'emp' = empirical target image E (from all fake windows' direct images, background and Gaia neighbours
removed) shifted by delta via the PRF:  E + PRF(x0 + delta) - PRF(x0);   'prf' = PRF(x0 + delta).
delta = 0, and 1", 3", 5" at PA 0, 45, ..., 315 deg.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3lib import *

NT = int(os.environ.get('V3_NT', 400))
RATIO = float(os.environ.get('V3_RATIO', 0.48))
SEPS = [1.0, 3.0, 5.0]
PAS = np.arange(0, 360, 45.0)
CONFIGS = [(0.0, 0.0)] + [(s, p) for s in SEPS for p in PAS]


def run(sec, seed=0):
    t0 = time.time()
    rng = np.random.default_rng(1000 + sec + seed)
    C = Cube(sec, centroid_cut=False)
    M = Model(C)
    z = np.load(f'{V3}/results/win_s{sec}.npz', allow_pickle=True); W = {k: z[k] for k in z.files}
    real = json.load(open(f'{V3}/results/real_s{sec}.json'))
    fk = np.where(W['kind'] == 'fake')[0]
    ep = W['n'][fk]
    Ball, Sall, _ = None, None, None
    Ds_, Ss_, Bs_ = stack(W, fk, np.ones(len(fk)))
    E, rB0 = M.empirical(Bs_, Ss_)
    P0 = M.img(rB0['x'], rB0['y'])
    Ji = np.linalg.inv(C.jac)
    out = dict(sector=sec, configs=CONFIGS, nt=NT, n_fake=int(len(fk)), n_fake_epochs=int(len(np.unique(ep))))
    arrs = {}
    for name, rel in (('deep', 1.0), ('shallow', RATIO)):
        n_real = int(real[name]['n']); A = float(real[name]['amp'])
        uniq = np.unique(ep)
        rec = {m: np.full((len(CONFIGS), NT, 2), np.nan) for m in ('emp', 'prf')}
        fcov = np.full((NT, 2, 2), np.nan)
        amps = np.full((len(CONFIGS), NT), np.nan)
        for t in range(NT):
            es = rng.choice(uniq, size=min(n_real, len(uniq)), replace=False)
            idx = np.array([rng.choice(fk[ep == e]) for e in es])
            Ds, Ss, Bs = stack(W, idx, np.full(len(idx), rel))
            rB = M.fit(Bs, Ss, direct=True)
            for ci, (sep, pa) in enumerate(CONFIGS):
                dsky = np.array([sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))])
                dpx = Ji @ dsky
                Pd = M.img(rB0['x'] + dpx[0], rB0['y'] + dpx[1])
                for mode, im in (('emp', E + Pd - P0), ('prf', Pd)):
                    if mode == 'prf' and sep > 0 and sep != 3.0:
                        continue   # PRF-only injections at 0" and 3" are enough to isolate noise-only behaviour
                    rD = M.fit(Ds + A * im, Ss, x0=rB['x'], y0=rB['y'])
                    sky = C.jac @ np.array([rD['x'] - rB['x'], rD['y'] - rB['y']])
                    rec[mode][ci, t] = sky - dsky
                    if mode == 'emp':
                        amps[ci, t] = rD['amp']
                        if ci == 0:
                            fcov[t] = C.jac @ rD['cov'] @ C.jac.T
        summ = []
        for ci, (sep, pa) in enumerate(CONFIGS):
            row = dict(sep=sep, pa=pa)
            for mode in ('emp', 'prf'):
                r = rec[mode][ci]
                if np.all(np.isnan(r)):
                    continue
                b = np.nanmean(r, 0)
                row[mode] = dict(bias=b, median=np.nanmedian(r, 0), sd=np.nanstd(r, 0), cov=np.cov(r[np.all(np.isfinite(r), 1)].T),
                                 r95=float(np.nanpercentile(np.hypot(*(r - b).T), 95)))
            row['amp_ratio'] = float(np.nanmean(amps[ci]) / A)
            summ.append(row)
        s0 = summ[0]
        fsd = np.sqrt(np.nanmedian(fcov[:, [0, 1], [0, 1]], 0))
        out[name] = dict(n_real=n_real, A=A, summary=summ, formal_sd_median=fsd)
        arrs[f'{name}_emp'] = rec['emp']; arrs[f'{name}_prf'] = rec['prf']; arrs[f'{name}_amps'] = amps
        print(f"S{sec} {name:7s} n={n_real} A={A:.1f}: target emp bias ({s0['emp']['bias'][0]:+.2f},{s0['emp']['bias'][1]:+.2f}) "
              f"sd ({s0['emp']['sd'][0]:.2f},{s0['emp']['sd'][1]:.2f}) | prf bias ({s0['prf']['bias'][0]:+.2f},{s0['prf']['bias'][1]:+.2f}) "
              f"sd ({s0['prf']['sd'][0]:.2f},{s0['prf']['sd'][1]:.2f}) | formal sd ({fsd[0]:.2f},{fsd[1]:.2f}) [{time.time() - t0:.0f}s]", flush=True)
        for sep in SEPS:
            rows = [r for r in summ if r['sep'] == sep]
            b = np.array([r['emp']['bias'] for r in rows]); sd = np.array([r['emp']['sd'] for r in rows])
            print(f"     {sep:.0f}\": bias |max| {np.max(np.hypot(*b.T)):.2f}\" mean ({b[:, 0].mean():+.2f},{b[:, 1].mean():+.2f}) "
                  f"sd mean ({sd[:, 0].mean():.2f},{sd[:, 1].mean():.2f}) amp ratio {np.mean([r['amp_ratio'] for r in rows]):.2f}", flush=True)
    out['empirical_direct'] = dict(x=rB0['x'], y=rB0['y'], F=rB0['amp'])
    np.savez_compressed(f'{V3}/results/inj_s{sec}.npz', configs=np.array(CONFIGS), **arrs)
    to_json(out, f'{V3}/results/inj_s{sec}.json')


def check(sec):
    """Cube-level injection for one trial at 3" PA 45 deg: inject into the frames, rebuild the windows, compare D."""
    C = Cube(sec)
    M = Model(C)
    z = np.load(f'{V3}/results/win_s{sec}.npz', allow_pickle=True)
    fk = np.where(z['kind'] == 'fake')[0][:5]
    dsky = np.array([3 * np.sin(np.radians(45)), 3 * np.cos(np.radians(45))])
    dpx = np.linalg.inv(C.jac) @ dsky
    im = M.img(C.xt + dpx[0], C.yt + dpx[1])
    A = 30.0
    maxdiff = 0
    for i in fk:
        tc = float(z['tc'][i])
        w = C.window(tc)
        D0, _, _, sb = diff_image(C, w)
        f_save = C.f.copy()
        C.f = C.f - A * vshape(C.t, tc, C.dt)[:, None, None] * im[None]
        D1, _, _, _ = diff_image(C, w)
        C.f = f_save
        maxdiff = max(maxdiff, np.max(np.abs((D1 - D0) - A * sb * im)))
    print(f'S{sec}: cube-level vs stacked-level injection, max |difference| = {maxdiff:.2e} e/s (signal peak {A * im.max():.1f})')


if __name__ == '__main__':
    if sys.argv[1] == '--check':
        check(int(sys.argv[2]))
    else:
        for s in [int(a) for a in sys.argv[1:]]:
            run(s)
