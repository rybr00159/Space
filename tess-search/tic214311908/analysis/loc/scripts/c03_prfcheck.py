"""How well does the TESS PRF model match isolated field stars in these cutouts?  For each bright calibrator
(T < 12.5, > 4 px from the Mira, > 2 px from the edges) fit PRF (free position, flux, constant) to the median image
within r <= 2.6 px and report the residual rms relative to the peak pixel and relative to total flux.  This is the
scale of PRF-model mismatch expected for the transit-signal image (whose noise per pixel is ~5-10% of its peak)."""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
from scipy.optimize import least_squares
g = stars(); out = {}
for sec in SECTORS:
    S = Sector(sec); psf = S.prf(); O = np.median(S.flux[S.good], 0)
    from l1img import background
    xy = S.star_xy(g)
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    im = int(np.where(g.source_id == IDS['Mira'])[0][0])
    rows = []
    for j in np.where(g['T'].values < 12.5)[0]:
        x0, y0 = xy[j]
        if not (2 <= x0 <= S.nx - 3 and 2 <= y0 <= S.ny - 3) or j == im: continue
        if np.hypot(x0 - xy[im, 0], y0 - xy[im, 1]) < 4: continue
        nbs = [k for k in np.where(g['T'].values < 16)[0] if k != j and np.hypot(*(xy[k] - xy[j])) < 4.5]
        mm = np.hypot(xx - x0, yy - y0) <= 2.6
        fixed = sum(g.F.values[k] * psf(*xy[k]) for k in nbs) if nbs else 0 * O
        def r1(p):
            A = np.column_stack([psf(x0 + p[0], y0 + p[1])[mm], np.ones(mm.sum())])
            c, *_ = np.linalg.lstsq(A, (O - fixed)[mm], rcond=None)
            return A @ c - (O - fixed)[mm]
        rr = least_squares(r1, [0, 0], diff_step=1e-3)
        A = np.column_stack([psf(x0 + rr.x[0], y0 + rr.x[1])[mm], np.ones(mm.sum())])
        c, *_ = np.linalg.lstsq(A, (O - fixed)[mm], rcond=None)
        res = rr.fun
        peak = (c[0] * A[:, 0]).max()
        rows.append(dict(source_id=int(g.source_id.values[j]), T=float(g['T'].values[j]), C=float(g.C.values[j]), flux=float(c[0]),
                         rms_over_peak=float(np.sqrt(np.mean(res ** 2)) / peak), maxabs_over_peak=float(np.abs(res).max() / peak),
                         rms_over_flux=float(np.sqrt(np.mean(res ** 2)) / c[0]), photon_rms_over_peak=float(np.sqrt(np.mean(O[mm])) / peak / np.sqrt(S.good.sum()) * 1.25),
                         contam=float(np.sum(fixed[mm]) / c[0]) if nbs else 0.0))
    df = pd.DataFrame(rows)
    out[str(sec)] = df.to_dict('records')
    print(f'S{sec}: ' + '; '.join(f"T{r['T']:.1f} C{r['C']:.2f}: rms/peak {r['rms_over_peak']*100:.1f}% max {r['maxabs_over_peak']*100:.1f}%" for r in rows))
to_json(out, f'{BASE}/out/prf_match_fieldstars.json')
