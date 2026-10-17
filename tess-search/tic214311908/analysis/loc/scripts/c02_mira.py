"""The Mira's apparent TESS position and PSF match, with the WCS calibrated on normal-colour stars (out/wcs_cal.json).

Direct (median) image, pixels within R px of the Mira's Gaia position: Mira as a PRF at a free position with free
flux; every other star at its calibrated apparent position (colour curve applied) with flux fixed to the catalogue
value (pair A+B, M3V scaled by a common free factor in variant 'nbfree'); the reference star free; background
constant.  sigma^2 = O + (f*O)^2.  Variants: R = 2.5/3.5/4.5 px, f = 0.01/0.02/0.05.  Also per-epoch fits
(median image of each 1-d chunk) to see whether the Mira's apparent position moves with its pulsation.
Output: out/mira_pos.json, figs/mira_direct_s{sec}.png
"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

g = stars()
out = {}


def fit_mira(O, S, psf, xy, mira_gaia, im, R=3.5, f=0.02, nbfree=False):
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mask = np.hypot(xx - mira_gaia[0], yy - mira_gaia[1]) <= R
    near = (np.hypot(xy[:, 0] - mira_gaia[0], xy[:, 1] - mira_gaia[1]) < R + 4) & (g['T'].values < 18)
    fixed = np.zeros_like(O); nbimg = np.zeros_like(O); ref = None
    for j in np.where(near)[0]:
        if j == im: continue
        img = psf(*xy[j])
        sid = g.source_id.values[j]
        if sid == IDS['ref']: ref = img
        elif nbfree and sid in (IDS['M3V'], IDS['pairA'], IDS['pairB']): nbimg += g.F.values[j] * img
        else: fixed += g.F.values[j] * img
    extra = [ref] if ref is not None else []
    if nbfree: extra.append(nbimg)
    def solve(p):
        A = np.array([psf(mira_gaia[0] + p[0], mira_gaia[1] + p[1])] + extra + [np.ones_like(O)])
        w = 1 / np.sqrt(np.abs(O[mask]) + (f * O[mask]) ** 2)
        M = A[:, mask].T * w[:, None]; b = (O - fixed)[mask] * w
        c, *_ = np.linalg.lstsq(M, b, rcond=None)
        return c, M @ c - b, A
    best = None
    for gx in np.linspace(-1.0, 1.0, 11):
        for gy in np.linspace(-1.0, 1.0, 11):
            c2 = (solve([gx, gy])[1] ** 2).sum()
            if best is None or c2 < best[0]: best = (c2, [gx, gy])
    r = least_squares(lambda p: solve(p)[1], best[1], diff_step=1e-3, x_scale=0.05)
    c, rr, A = solve(r.x)
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * max(1, (rr ** 2).sum() / (mask.sum() - len(c) - 2))
    mod = np.tensordot(c, A, axes=1) + fixed
    return dict(p=r.x, cov=cov, c=c, mod=mod, mask=mask, chi2=float((rr ** 2).sum()), n=int(mask.sum()))


for sec in SECTORS:
    S = Sector(sec); psf = S.prf()
    xy = S.star_xy(g, colour=True)
    xy_nc = S.star_xy(g, colour=False)
    im = int(np.where(g.source_id == IDS['Mira'])[0][0])
    mira_gaia = xy_nc[im]
    O = np.median(S.flux[S.good], 0)
    res = {}
    for R in (2.5, 3.5, 4.5):
        for f in (0.01, 0.02, 0.05):
            for nbf in (False, True):
                F = fit_mira(O, S, psf, xy, mira_gaia, im, R, f, nbf)
                dEN = S.J @ F['p']; e = np.sqrt(np.diag(S.J @ F['cov'] @ S.J.T))
                res[f'R{R}_f{f}_nb{int(nbf)}'] = dict(dEN=dEN.tolist(), err=e.tolist(), rad=float(dEN @ S.u_rad),
                                                      tan=float(dEN @ np.array([S.u_rad[1], -S.u_rad[0]])),
                                                      flux=float(F['c'][0]), nbscale=float(F['c'][2]) if nbf else 1.0,
                                                      chi2r=F['chi2'] / (F['n'] - 6),
                                                      resid_rms_frac=float(np.sqrt(np.mean(((O - F['mod'])[F['mask']] / O[F['mask']]) ** 2))))
    base = res['R3.5_f0.02_nb0']
    allE = np.array([v['dEN'] for v in res.values()])
    # per-chunk (1-d) fits: does the apparent position move?
    gd = S.good; t = S.t[gd]
    chunks = []
    for t0 in np.arange(t.min(), t.max(), 2.0):
        m = (t >= t0) & (t < t0 + 2.0)
        if m.sum() < 10: continue
        Oc = np.median(S.flux[gd][m], 0)
        F = fit_mira(Oc, S, psf, xy, mira_gaia, im, 3.5, 0.02, False)
        d = S.J @ F['p']
        chunks.append([float(t[m].mean()), float(d[0]), float(d[1]), float(F['c'][0])])
    chunks = np.array(chunks)
    out[str(sec)] = dict(variants=res, base=base, variant_range_E=[float(allE[:, 0].min()), float(allE[:, 0].max())],
                         variant_range_N=[float(allE[:, 1].min()), float(allE[:, 1].max())],
                         chunks=chunks.tolist(), chunk_sd=[float(chunks[:, 1].std()), float(chunks[:, 2].std())],
                         chunk_flux_range=[float(chunks[:, 3].min()), float(chunks[:, 3].max())],
                         predicted_colour_shift_rad=float(colour_shift(4.68, S.colscale)), u_rad=S.u_rad.tolist())
    print(f"S{sec}: Mira apparent offset from Gaia (base) ({base['dEN'][0]:+.1f},{base['dEN'][1]:+.1f})\" +- ({base['err'][0]:.1f},{base['err'][1]:.1f}); "
          f"radial {base['rad']:+.1f}\" tangential {base['tan']:+.1f}\"; flux {base['flux']:.0f}; resid rms {base['resid_rms_frac']*100:.1f}%; "
          f"variants E [{allE[:,0].min():+.1f},{allE[:,0].max():+.1f}] N [{allE[:,1].min():+.1f},{allE[:,1].max():+.1f}]; "
          f"2-d chunks sd ({chunks[:,1].std():.1f},{chunks[:,2].std():.1f})\", flux {chunks[:,3].min():.0f}-{chunks[:,3].max():.0f}; "
          f"colour-curve extrapolation predicts radial {colour_shift(4.68, S.colscale):+.1f}\"")
    # figure
    F = fit_mira(O, S, psf, xy, mira_gaia, im, 3.5, 0.02, False)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    sl = (slice(max(0, int(mira_gaia[1]) - 5), int(mira_gaia[1]) + 6), slice(max(0, int(mira_gaia[0]) - 5), int(mira_gaia[0]) + 6))
    ext = [sl[1].start - 0.5, sl[1].stop - 0.5, sl[0].start - 0.5, sl[0].stop - 0.5]
    ax[0].imshow(np.log10(O[sl]), origin='lower', extent=ext); ax[0].set_title(f'S{sec} median image (log)')
    rimg = ((O - F['mod']) / O)[sl] * 100
    v = np.nanmax(np.abs(rimg[F['mask'][sl]]))
    h = ax[1].imshow(np.where(F['mask'][sl], rimg, np.nan), origin='lower', cmap='RdBu_r', vmin=-v, vmax=v, extent=ext)
    plt.colorbar(h, ax=ax[1]); ax[1].set_title('residual / data (%)')
    for a in ax[:2]:
        for k, c in (('Mira', 'r'), ('M3V', 'c'), ('pairA', 'm'), ('ref', 'y')):
            j = int(np.where(g.source_id == IDS[k])[0][0]); a.plot(*xy_nc[j], '+', color=c, ms=10, mew=2)
        a.plot(mira_gaia[0] + F['p'][0], mira_gaia[1] + F['p'][1], 'x', color='k', ms=10, mew=2)
    ax[2].plot(chunks[:, 0], chunks[:, 1], 'o-', label='dE'); ax[2].plot(chunks[:, 0], chunks[:, 2], 's-', label='dN')
    ax2 = ax[2].twinx(); ax2.plot(chunks[:, 0], chunks[:, 3], 'k:', label='flux'); ax[2].legend(); ax[2].set_xlabel('BTJD')
    ax[2].set_ylabel('Mira offset from Gaia (arcsec)')
    fig.tight_layout(); fig.savefig(f'{BASE}/figs/mira_direct_s{sec}.png', dpi=80); plt.close(fig)
to_json(out, f'{BASE}/out/mira_pos.json')
