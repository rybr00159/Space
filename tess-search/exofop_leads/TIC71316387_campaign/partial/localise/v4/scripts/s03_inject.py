"""v4 injection-recovery calibration of the per-sector offsets.

Trial: n_real fake windows (distinct orbits, off-eclipse phases, real frames) are stacked exactly like the real ones;
the direct image is fitted; an eclipse with the real stacked amplitude a_real is injected at sky offset delta from the
G star and the amplitude image is fitted with the identical code.  Injection is done on the per-window amplitude images,
which is exactly equivalent to injecting a_real*S(t)*I into every cadence (the matched filter is linear and its
coefficient on S is exactly 1); --check verifies this at the cube level.
Injected image I: 'emp' = empirical target image E (direct image of all fake windows, background and tied Gaia
neighbours removed), shifted by delta through the PRF difference: E + PRF(x0+d) - PRF(x0);  'prf' = PRF(x0+d).
Configurations: 0", and 1", 3", 5" at PA 0..315 step 45 deg.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4lib import *

NT = int(os.environ.get('V4_NT', 500))
TAG = os.environ.get('V4_TAG', '')
SEPS = [1.0, 3.0, 5.0]
PAS = np.arange(0, 360, 45.0)
CONFIGS = [(0.0, 0.0)] + [(s, p) for s in SEPS for p in PAS]


def run(sec, seed=0):
    t0 = time.time()
    rng = np.random.default_rng(7000 + sec + 100 * seed)
    S = Sector(sec)
    G = PRFGrid(S)
    z = np.load(f'{V4}/results/win_s{sec}{TAG}.npz', allow_pickle=True); W = {k: z[k] for k in z.files}
    real = json.load(open(f'{V4}/results/real_s{sec}{TAG}.json'))
    fk = np.where(W['kind'] == 'fake')[0]
    ep = W['ep'][fk]
    uniq = np.unique(ep)
    Af, sf, Bf = stack(W, fk)
    E, rE, cE = G.empirical(Bf, sf)
    x0, y0 = rE['x'], rE['y']
    P0 = G.img(x0, y0)
    Ji = np.linalg.inv(S.J)
    out = dict(sector=sec, configs=CONFIGS, nt=NT, n_fake=int(len(fk)), n_fake_epochs=int(len(uniq)), E_pos=[x0, y0])
    arrs = {}
    for name in ('deep', 'shallow'):
        n_real = int(real[name]['n']); a_real = float(real[name]['amp'])
        nsel = min(n_real, len(uniq))
        rec = {m: np.full((len(CONFIGS), NT, 2), np.nan) for m in ('emp', 'prf')}
        amp = np.full((len(CONFIGS), NT), np.nan)
        for t in range(NT):
            es = rng.choice(uniq, size=nsel, replace=False)
            idx = np.array([rng.choice(fk[ep == e]) for e in es])
            A, sA, B = stack(W, idx)
            rB = G.fit(B, sA, direct=True)
            for ci, (sep, pa) in enumerate(CONFIGS):
                dsky = np.array([sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))])
                dpx = Ji @ dsky
                Pd = G.img(x0 + dpx[0], y0 + dpx[1])
                for mode in ('emp', 'prf'):
                    if mode == 'prf' and sep not in (0.0, 3.0):
                        continue
                    I = E + Pd - P0 if mode == 'emp' else Pd
                    r = G.fit(A + a_real * I, sA)
                    rec[mode][ci, t] = S.J @ np.array([r['x'] - rB['x'], r['y'] - rB['y']]) - dsky
                    if mode == 'emp':
                        amp[ci, t] = r['a'] / a_real
        summ = []
        for ci, (sep, pa) in enumerate(CONFIGS):
            row = dict(sep=sep, pa=pa, amp_ratio=float(np.nanmean(amp[ci])))
            for mode in ('emp', 'prf'):
                r = rec[mode][ci]; ok = np.all(np.isfinite(r), 1)
                if ok.sum() < 10:
                    continue
                r = r[ok]
                b = r.mean(0)
                row[mode] = dict(bias=b, median=np.median(r, 0), sd=r.std(0), cov=np.cov(r.T),
                                 r68=float(np.percentile(np.hypot(*(r - b).T), 68.27)), r95=float(np.percentile(np.hypot(*(r - b).T), 95.45)))
            summ.append(row)
        out[name] = dict(n_real=n_real, a_real=a_real, summary=summ)
        arrs[f'{name}_emp'] = rec['emp']; arrs[f'{name}_prf'] = rec['prf']; arrs[f'{name}_amp'] = amp
        s0 = summ[0]
        print(f"S{sec} {name:7s} n={n_real} a={a_real:.1f} | target emp bias ({s0['emp']['bias'][0]:+.2f},{s0['emp']['bias'][1]:+.2f}) "
              f"sd ({s0['emp']['sd'][0]:.2f},{s0['emp']['sd'][1]:.2f}) r95 {s0['emp']['r95']:.2f} | prf bias ({s0['prf']['bias'][0]:+.2f},"
              f"{s0['prf']['bias'][1]:+.2f}) sd ({s0['prf']['sd'][0]:.2f},{s0['prf']['sd'][1]:.2f}) [{time.time() - t0:.0f}s]", flush=True)
        for sep in SEPS:
            rows = [r for r in summ if r['sep'] == sep]
            b = np.array([r['emp']['bias'] for r in rows]); sd = np.array([r['emp']['sd'] for r in rows])
            print(f"    {sep:.0f}\": |bias| max {np.max(np.hypot(*b.T)):.2f}\" mean ({b[:, 0].mean():+.2f},{b[:, 1].mean():+.2f}) "
                  f"sd mean ({sd[:, 0].mean():.2f},{sd[:, 1].mean():.2f}) amp ratio {np.mean([r['amp_ratio'] for r in rows]):.3f}", flush=True)
    np.savez_compressed(f'{V4}/results/inj_s{sec}{TAG}.npz', configs=np.array(CONFIGS), **arrs)
    to_json(out, f'{V4}/results/inj_s{sec}{TAG}.json')


def check(sec):
    """Cube-level injection into 6 fake windows (3" at PA 45): rebuild amplitude images from the modified cube."""
    S = Sector(sec)
    G = PRFGrid(S, half=0.3)
    z = np.load(f'{V4}/results/win_s{sec}{TAG}.npz', allow_pickle=True)
    fk = np.where(z['kind'] == 'fake')[0][::7][:6]
    dsky = np.array([3 * np.sin(np.radians(45)), 3 * np.cos(np.radians(45))])
    dpx = np.linalg.inv(S.J) @ dsky
    I = G.img(S.xt + dpx[0], S.yt + dpx[1])
    a = 30.0
    md = 0
    for i in fk:
        tc = float(z['tc'][i])
        w0 = S.window_fit(tc)
        f0 = S.f.copy()
        S.f = S.f - a * template(S.t - tc, S.texp)[:, None, None] * I[None]
        w1 = S.window_fit(tc)
        S.f = f0
        md = max(md, float(np.max(np.abs((w1['A'] - w0['A']) - a * I))))
    print(f'S{sec}: cube-level vs image-level injection: max |dA - a*I| = {md:.2e} e/s (peak {a * I.max():.1f} e/s)')


if __name__ == '__main__':
    if sys.argv[1] == '--check':
        check(int(sys.argv[2]))
    else:
        for s in [int(a) for a in sys.argv[1:]]:
            run(s)
