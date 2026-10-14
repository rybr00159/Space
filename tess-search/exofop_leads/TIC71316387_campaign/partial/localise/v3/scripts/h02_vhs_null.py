"""VHS null test: run the identical saturated-core companion scan (h01) on the target and on other stars of similar
brightness on the same VIRCAM chip (N228: Ks 10.43 / J 10.96, the closest brightness match; N169: Ks 11.77; N314: Ks 11.98,
unsaturated), using larger SODA cutouts of the four pawprints (2 J, 2 Ks).  ePSF for each analysed star is built from
unsaturated isolated Gaia stars within 200" of it.  If the same 'companion' (sep, PA) appears around the null stars,
the target's best-fit companion is a PSF/saturation artefact.
Also: injection-recovery around N228 (fake companions in its real stamp) as a second sensitivity estimate.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from h01_vhs import build_epsf, psf, fit_single, Scanner, stamp, local_sky, H, GAIN, RN, SATLIM, EP
from v3lib import V3, to_json
import warnings; warnings.filterwarnings('ignore')

gw = Table.read(V3 + '/gaia/wide_2p7deg_G18p5.fits').to_pandas()
c0 = np.cos(np.radians(-42.84))
gw = gw[(np.abs(gw.dec + 42.86) < 0.3) & (np.abs((gw.ra - 330.85) * c0) < 0.3)].reset_index(drop=True)
STARS = {'target': 6571343840265900544}
for name, (ra, de) in {'N228': (330.9009, -42.7729), 'N169': (330.8423, -42.8027), 'N314': (330.8518, -42.7524),
                       'N132': (330.8750, -42.8708), 'N217': (330.8302, -42.8786)}.items():
    d = np.hypot((gw.ra - ra) * c0, gw.dec - de) * 3600
    STARS[name] = int(gw.source_id[np.argmin(d)])


def analyse(f, name, sid, inject=False):
    h = fits.open(f); ext = [e for e in h if e.data is not None and e.data.shape[1] > 100][0]
    d = ext.data.astype(float); w = WCS(ext.header)
    ra = gw.ra + np.nan_to_num(gw.pmra) * (EP - 2016) / 3.6e6 / np.cos(np.radians(gw.dec))
    de = gw.dec + np.nan_to_num(gw.pmdec) * (EP - 2016) / 3.6e6
    x, y = w.all_world2pix(ra.values, de.values, 0)
    I = dict(d=d, x=x, y=y)
    it = int(np.where(gw.source_id.values == sid)[0][0])
    r0 = w.all_pix2world([[x[it], y[it]], [x[it] + 1, y[it]], [x[it], y[it] + 1]], 0)
    cc = np.cos(np.radians(r0[0, 1]))
    J = np.array([[(r0[1, 0] - r0[0, 0]) * cc * 3600, (r0[2, 0] - r0[0, 0]) * cc * 3600], [(r0[1, 1] - r0[0, 1]) * 3600, (r0[2, 1] - r0[0, 1]) * 3600]])
    pix = float(np.sqrt(abs(np.linalg.det(J)))); Jn = J / pix
    G = gw.phot_g_mean_mag.values
    cand = []
    for i in range(len(gw)):
        if i == it or G[i] > 17.3 or np.hypot(x[i] - x[it], y[i] - y[it]) * pix > 200:
            continue
        if not (H + 45 < x[i] < d.shape[1] - H - 45 and H + 45 < y[i] < d.shape[0] - H - 45):
            continue
        dd = np.hypot((gw.ra - gw.ra[i]) * c0, gw.dec - gw.dec[i]) * 3600
        if np.any((dd > 0.01) & (dd < 6) & (G - G[i] < 5)):
            continue
        st, _, _ = stamp(d, x[i], y[i])
        if st.max() > 20000 or st.max() - local_sky(d, x[i], y[i]) < 800:
            continue
        cand.append(i)
    E = build_epsf(I, cand)
    raw, _, _ = stamp(d, x[it], y[it]); sky = local_sky(d, x[it], y[it]); st = raw - sky
    sig = np.sqrt(np.clip(raw, 1, None) / GAIN + RN ** 2)
    yy, xx = np.mgrid[-H:H + 1, -H:H + 1]
    mask = (raw < SATLIM) & (np.hypot(xx, yy) <= H)
    p0, cf, chi = fit_single(st, sig, mask, E)
    SC = Scanner(E, p0, pix, rmin=0.6, rmax=4.0, step=0.1)
    b, chi0, c0f = SC.run(st, sig, mask)
    so = Jn @ np.array(b[1]); sep = float(np.hypot(*b[1])); pa = float(np.degrees(np.arctan2(so[0], so[1])) % 360)
    nsat = int((~(raw < SATLIM) & (np.hypot(xx, yy) <= H)).sum())
    out = dict(file=os.path.basename(f), star=name, sid=sid, n_psf=len(cand), peak_raw=float(raw.max()), mask_r=float(np.sqrt(nsat / np.pi) * pix),
               chi2r=chi0 / (mask.sum() - 3), dchi2=b[0], dchi2_over_chi2r=b[0] / (chi0 / (mask.sum() - 3)), sep=sep, pa=pa,
               dm=float(-2.5 * np.log10(b[2])) if b[2] > 0 else None)
    if inject:
        rng = np.random.default_rng(9); F = c0f[0]; rec = []
        for sp in (1.0, 1.5, 2.0):
            for dm in (2.0, 2.5, 3.0, 3.5, 4.0):
                ok = []
                for pa_ in np.arange(0, 360, 60.0) + rng.uniform(0, 60):
                    dpx = np.linalg.solve(Jn, sp * np.array([np.sin(np.radians(pa_)), np.cos(np.radians(pa_))]))
                    comp = F * 10 ** (-0.4 * dm) * psf(E, p0[0] + dpx[0] / pix, p0[1] + dpx[1] / pix)
                    raw2 = raw + comp; mk2 = (raw2 < SATLIM) & (np.hypot(xx, yy) <= H)
                    sg2 = np.sqrt(np.clip(raw2, 1, None) / GAIN + RN ** 2)
                    bb, _, _ = SC.run(st + comp, sg2, mk2)
                    ok.append(bool(np.hypot(bb[1][0] - dpx[0], bb[1][1] - dpx[1]) < 0.35 and bb[0] > b[0]))
                rec.append(dict(sep=sp, dm=dm, frac=float(np.mean(ok))))
        out['inj'] = rec
    print(f"{out['file'][4:30]} {name:6s} peak {out['peak_raw']:.0f} mask r {out['mask_r']:.2f}\" npsf {len(cand)} chi2r {out['chi2r']:.1f} "
          f"best: dchi2 {b[0]:.0f} (/chi2r {out['dchi2_over_chi2r']:.0f}) sep {sep:.2f}\" PA {pa:.0f} dm {out['dm']}", flush=True)
    return out


if __name__ == '__main__':
    res = []
    for f in sorted(glob.glob(V3 + '/imaging/big/vhs_*.fits.fz')):
        for name, sid in STARS.items():
            try:
                res.append(analyse(f, name, sid, inject=(name in ('N228',))))
            except Exception as e:
                print('fail', f, name, e)
    to_json(res, V3 + '/results/h02_vhs_null.json')
