"""v4 lateral-colour analysis of the f01 field-star tables.

Per FFI, residuals r = fitted - WCS(Gaia) (px -> arcsec with the local Jacobian at each star's position, approximated by
the Jacobian at the target) are modelled as
    r = poly2(x, y) [each sky component]  +  m * (G - 13) [each component]
        + (C - C0) * (k_rad * u_rad + k_tan * u_tan)
with C = BP-RP, C0 = 0.786 (the G star), u_rad = unit sky vector from the star TOWARDS the camera boresight.
k_rad < 0  <=>  redder stars are displaced AWAY from the boresight.
Stars used: within RFIT deg of the target (default 1.6 deg: local, no field-angle scaling; variant 'theta' uses 2.6 deg
with the colour term scaled by theta/theta_target).  Robust iterative 4-sigma clipping; errors from a bootstrap over
stars (and, separately, from the agreement of the two FFIs of each sector).
Also: binned radial-residual colour curve (after removing poly + magnitude terms) per sector and combined, relative to
the G star's colour, for dwarfs (parallax-based M_G > 4.5) and all stars, to map lever arm -> colour at the red end.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v4lib import V4, to_json

C0 = 0.786
RFIT = float(os.environ.get('F02_RFIT', 1.6))
SCALE = os.environ.get('F02_SCALE', 'none')       # none | theta
GMAXFIT = float(os.environ.get('F02_GMAX', 15.5))
TAG = os.environ.get('F02_TAG', '')
FLOOR = 0.10


def load(f):
    z = np.load(f)
    a = z['stars']
    ok = (a['chi2r'] < 25) & (a['ruwe'] < 1.4) & (a['flux'] > 0) & (a['ex'] < 0.2) & (a['ey'] < 0.2) & (a['peak'] < 1.5e5)
    a = a[ok]
    J = z['J']
    dsky = J @ np.vstack([a['xf'] - a['xw'], a['yf'] - a['yw']])
    esky = np.sqrt((J ** 2) @ np.vstack([a['ex'] ** 2, a['ey'] ** 2]))
    pa = np.radians(a['pa_bs'])
    ur = np.vstack([np.sin(pa), np.cos(pa)])
    ut = np.vstack([ur[1], -ur[0]])
    return dict(a=a, d=dsky, e=esky, ur=ur, ut=ut, th=a['theta'], tht=float(z['theta_t']), sec=int(z['sec']), xt=float(z['xt']),
                yt=float(z['yt']), file=os.path.basename(f), pa_t=float(z['pa_t']))


def design(F, scale):
    a = F['a']; N = len(a)
    x = (a['xw'] - F['xt']) / 500.0; y = (a['yw'] - F['yt']) / 500.0
    poly = [np.ones(N), x, y, x * x, x * y, y * y]
    Z = np.zeros(N)
    sc = F['th'] / F['tht'] if scale == 'theta' else np.ones(N)
    dc = (a['C'] - C0) * sc
    dg = a['G'] - 13.0
    cols = [np.r_[p, Z] for p in poly] + [np.r_[Z, p] for p in poly] + [np.r_[dg, Z], np.r_[Z, dg]]
    cols += [np.r_[dc * F['ur'][0], dc * F['ur'][1]], np.r_[dc * F['ut'][0], dc * F['ut'][1]]]
    return np.column_stack(cols), np.r_[F['d'][0], F['d'][1]], np.sqrt(np.r_[F['e'][0], F['e'][1]] ** 2 + FLOOR ** 2), sc


def rfit(A, b, e, use):
    n = len(b) // 2
    ok = use.copy()
    for _ in range(15):
        w = 1 / e[ok]
        p, *_ = np.linalg.lstsq(A[ok] * w[:, None], b[ok] * w, rcond=None)
        r = (b - A @ p) / e
        s = 1.4826 * np.median(np.abs(r[ok]))
        g = np.abs(r) < 4 * s
        g = g[:n] & g[n:]
        new = use & np.r_[g, g]
        if np.array_equal(new, ok):
            break
        ok = new
    return p, ok, s


def main():
    frames = [load(f) for f in sorted(glob.glob(f'{V4}/results/f01_tess*.npz'))]
    rng = np.random.default_rng(5)
    res, pts = [], []
    for F in frames:
        rmax = RFIT if SCALE == 'none' else 2.6
        A, b, e, sc = design(F, SCALE)
        N = len(F['a'])
        sel = (F['a']['rt'] < rmax) & (F['a']['G'] < GMAXFIT)
        use = np.r_[sel, sel]
        p, ok, s = rfit(A, b, e, use)
        boots = []
        for _ in range(300):
            i = rng.integers(0, N, N); ii = np.r_[i, i + N]
            pb, _, _ = rfit(A[ii], b[ii], e[ii], use[ii])
            boots.append(pb[-2:])
        boots = np.array(boots)
        rr = b - A[:, :-2] @ p[:-2]
        rad = (rr[:N] * F['ur'][0] + rr[N:] * F['ur'][1]) / sc
        tan = (rr[:N] * F['ut'][0] + rr[N:] * F['ut'][1]) / sc
        okk = ok[:N]
        res.append(dict(file=F['file'], sec=F['sec'], n_used=int(okk.sum()), n_red2=int((okk & (F['a']['C'] > 2)).sum()),
                        theta_t=F['tht'], pa_t=F['pa_t'], k_rad=p[-2], k_rad_err=boots[:, 0].std(), k_tan=p[-1], k_tan_err=boots[:, 1].std(),
                        m=[p[12], p[13]], clip_scale=s))
        a = F['a']
        absG = a['G'] + 5 * np.log10(np.clip(a['plx'], 1e-3, None) / 100.0)
        pts.append(dict(sid=a['sid'], C=a['C'], G=a['G'], rad=rad, tan=tan, ok=okk | ((a['G'] >= GMAXFIT) & (a['rt'] < rmax)),
                        fitused=okk, sec=np.full(N, F['sec']), dwarf=(a['plx'] / np.clip(a['eplx'], 1e-3, None) > 5) & (absG > 4.5),
                        rt=a['rt'], th=a['theta'], e=np.hypot(F['e'][0], F['e'][1]) / np.sqrt(2)))
        r = res[-1]
        print(f"{F['file'][4:24]} S{F['sec']:3d} th_t={F['tht']:.2f} n={r['n_used']:4d} (C>2: {r['n_red2']:2d}) k_rad={r['k_rad']:+.3f}+-{r['k_rad_err']:.3f} "
              f"k_tan={r['k_tan']:+.3f}+-{r['k_tan_err']:.3f} \"/mag  m=({p[12]:+.3f},{p[13]:+.3f}) clip {s:.2f}", flush=True)
    # per sector (two FFIs): mean and FFI-to-FFI difference
    secs = sorted(set(r['sec'] for r in res))
    persec = {}
    for s in secs:
        rs = [r for r in res if r['sec'] == s]
        k = np.array([r['k_rad'] for r in rs]); ek = np.array([r['k_rad_err'] for r in rs])
        kt = np.array([r['k_tan'] for r in rs]); ekt = np.array([r['k_tan_err'] for r in rs])
        w = 1 / ek ** 2
        persec[s] = dict(k_rad=float((k * w).sum() / w.sum()), k_rad_err=float(1 / np.sqrt(w.sum())), k_rad_ffi=k.tolist(),
                         ffi_diff_sig=float(abs(k[0] - k[-1]) / np.hypot(ek[0], ek[-1])) if len(k) > 1 else None,
                         k_tan=float((kt / ekt ** 2).sum() / (1 / ekt ** 2).sum()), k_tan_err=float(1 / np.sqrt((1 / ekt ** 2).sum())),
                         pa_t=rs[0]['pa_t'], theta_t=rs[0]['theta_t'])
        print(f"  S{s:3d}: k_rad = {persec[s]['k_rad']:+.3f} +- {persec[s]['k_rad_err']:.3f} (FFIs {np.round(k, 3)}, diff {persec[s]['ffi_diff_sig']:.1f} sig)  "
              f"k_tan = {persec[s]['k_tan']:+.3f} +- {persec[s]['k_tan_err']:.3f}")
    K = np.array([persec[s]['k_rad'] for s in secs]); EK = np.array([persec[s]['k_rad_err'] for s in secs])
    kw = (K / EK ** 2).sum() / (1 / EK ** 2).sum()
    chi = (((K - kw) / EK) ** 2).sum()
    print(f'  weighted mean k_rad = {kw:+.3f} +- {1 / np.sqrt((1 / EK ** 2).sum()):.3f}; chi2 between sectors = {chi:.1f}/{len(K) - 1}')
    # colour curve (stars averaged over the two FFIs of a sector first)
    allp = {k: np.concatenate([p[k] for p in pts]) for k in pts[0]}
    bins = [(-0.2, 0.55), (0.55, 0.7), (0.7, 0.85), (0.85, 1.0), (1.0, 1.2), (1.2, 1.45), (1.45, 1.75), (1.75, 2.05), (2.05, 2.35),
            (2.35, 2.65), (2.65, 3.0), (3.0, 3.5), (3.5, 5.0)]
    curve = []
    key = allp['sid'].astype(np.int64) * 1000 + allp['sec']
    print('\n  radial colour curve (arcsec; + towards boresight; rel. to linear-model zero at C0); stars x sectors:')
    for lo, hi in bins:
        m = allp['ok'] & (allp['C'] >= lo) & (allp['C'] < hi)
        if m.sum() < 2:
            continue
        # average the two FFIs of the same star/sector
        u, inv = np.unique(key[m], return_inverse=True)
        rad = np.bincount(inv, allp['rad'][m]) / np.bincount(inv)
        tan = np.bincount(inv, allp['tan'][m]) / np.bincount(inv)
        cc = np.bincount(inv, allp['C'][m]) / np.bincount(inv)
        dw = np.bincount(inv, allp['dwarf'][m].astype(float)) / np.bincount(inv) > 0.5
        med = float(np.median(rad)); mad = 1.4826 * np.median(np.abs(rad - med))
        err = float(1.2533 * mad / np.sqrt(len(u))) if len(u) > 2 else float(np.std(rad) / np.sqrt(max(len(u) - 1, 1)))
        row = dict(lo=lo, hi=hi, n=int(len(u)), C=float(np.median(cc)), rad=med, err=err, tan=float(np.median(tan)), lin=float(kw * (np.median(cc) - C0)),
                   n_dwarf=int(dw.sum()), rad_dwarf=float(np.median(rad[dw])) if dw.sum() >= 2 else None,
                   rad_giant=float(np.median(rad[~dw])) if (~dw).sum() >= 2 else None)
        curve.append(row)
        print(f"   C {lo:4.2f}-{hi:4.2f} n={row['n']:4d} medC={row['C']:.2f} rad={med:+.3f}+-{err:.3f} (lin {row['lin']:+.3f}) tan={row['tan']:+.3f} "
              f"dwarfs {row['n_dwarf']}: {row['rad_dwarf'] if row['rad_dwarf'] is None else round(row['rad_dwarf'], 3)}  "
              f"others: {row['rad_giant'] if row['rad_giant'] is None else round(row['rad_giant'], 3)}")
    # reddest individual stars
    red = allp['ok'] & (allp['C'] > 2.6)
    u, inv = np.unique(key[red], return_inverse=True)
    reds = []
    for j, kk in enumerate(u):
        mm = np.where(red)[0][inv == j]
        reds.append(dict(sid=int(kk // 1000), sec=int(kk % 1000), C=float(allp['C'][mm[0]]), G=float(allp['G'][mm[0]]),
                         rad=float(allp['rad'][mm].mean()), tan=float(allp['tan'][mm].mean()), e=float(allp['e'][mm].mean() / np.sqrt(len(mm))),
                         dwarf=bool(allp['dwarf'][mm[0]])))
    reds.sort(key=lambda r: r['C'])
    print('\n  stars with C > 2.6:')
    for r in reds:
        print(f"   {r['sid']} S{r['sec']} C={r['C']:.2f} G={r['G']:.2f} rad={r['rad']:+.2f}+-{r['e']:.2f} tan={r['tan']:+.2f} dwarf={r['dwarf']}")
    to_json(dict(scale=SCALE, rfit=RFIT, gmaxfit=GMAXFIT, frames=res, sectors=persec, k_rad_wmean=kw, chi2_sectors=chi, curve=curve, reds=reds),
            f'{V4}/results/f02_colour{TAG}.json')
    np.savez_compressed(f'{V4}/results/f02_points{TAG}.npz', **allp)


if __name__ == '__main__':
    main()
