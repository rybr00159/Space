"""f0: per-cadence raw-pixel 3x3-box aperture light curves around each transit (phase 0) and at phase 0.5, for the
fold figure. This re-implements the referee3_clean q02 'box3, pointing-regressor' measurement (raw TESScut pixels,
3x3 box at the M3V pixel, faint-pixel background, pointing from isolated bright stars, per-window quadratic baseline +
pointing terms fitted jointly with a transit-shaped template), because q02 saved only per-window flux losses, not the
per-cadence curves. Encircled energy (PRF units) is taken from referee3_clean/out/q02_bigap.json.
Check: the per-window flux losses reproduce q02_bigap_per_transit.csv (L_box3_pt) - printed at the end.
Output: figscripts/data/rawap_box3.npz (small)."""
import sys, glob, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd, batman
from astropy.io import fits
from astropy.wcs import WCS

R3 = f'{B}/referee3_clean'
T0, P = 3840.91488, 5.6915020          # the ephemeris q02 used (window centres); fold uses the adopted one later
K, B_, AR, U = 0.117, 0.52, 25.3, (0.125, 0.533)
IDS = dict(Mira=5925373268528684288, M3V=5925373272850784768)
HALF = 0.30


def shape(t, tc, ex):
    ns = 5 if ex * 86400 < 300 else (9 if ex * 86400 < 900 else 21)
    off = ((np.arange(ns) + 0.5) / ns - 0.5) * ex
    tt = (t[:, None] + off[None]).ravel()
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = tc, P, K, AR, np.degrees(np.arccos(B_ / AR)), 0., 90., list(U), 'quadratic'
    f = batman.TransitModel(p, tt).light_curve(p).reshape(-1, ns).mean(1)
    p0 = batman.TransitModel(p, np.array([tc])).light_curve(p)[0]
    return (1 - f) / (1 - p0)


def accept(t, tc, ex):
    tau = t - tc
    oot = np.abs(tau) > 0.07
    nexp = 2 * (HALF - 0.07) / ex
    return (np.sum(oot & (tau < 0)) >= 0.25 * nexp / 2) and (np.sum(oot & (tau > 0)) >= 0.25 * nexp / 2) and np.sum(np.abs(tau) < 0.028) >= 2


def main():
    gaia = pd.read_csv(f'{B}/phot/data/gaia/gaia_360arcsec.csv')
    var = json.load(open(f'{B}/search/out/gaia_variables.json'))['variables']
    var_ids = set(int(v['source_id']) for v in var) - {IDS['Mira']}
    addtab = np.load(f'{R3}/out/q01_addtab.npz')
    S = json.load(open(f'{R3}/out/q02_bigap.json'))
    cols = {k: [] for k in ('sector', 'epoch', 'phase', 't', 'dt', 'f', 'e', 'cad', 'L')}
    chk = []
    for sec in (12, 39, 66, 93, 103, 104):
        fn = glob.glob(f'{CUBES}/tess-s{sec:04d}-*_astrocut.fits')[0]
        with fits.open(fn) as h:
            d = h[1].data; w = WCS(h[2].header)
            t = np.array(d['TIME'], float); F = np.array(d['FLUX'], np.float64); q = np.array(d['QUALITY'], int)
        ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
        t, F = t[ok], F[ok]
        t = t + np.interp(t, addtab[f't_{sec}'], addtab[f'add_{sec}']) / 86400.0
        ex = float(np.median(np.diff(t)))
        ny, nx = F.shape[1:]
        ep = 2000.0 + (np.median(t) + 2457000 - 2451545) / 365.25
        g = gaia.copy(); dtt = ep - 2016.0
        g['ra_e'] = g.ra + np.nan_to_num(g.pmra) * dtt / 3.6e6 / np.cos(np.radians(g.dec))
        g['dec_e'] = g.dec + np.nan_to_num(g.pmdec) * dtt / 3.6e6
        g['x'], g['y'] = w.all_world2pix(g.ra_e.values, g.dec_e.values, 0)
        pos = {k: g[g.source_id == v][['x', 'y']].values[0] for k, v in IDS.items()}
        yy, xx = np.mgrid[0:ny, 0:nx]
        dmira = np.hypot(xx - pos['Mira'][0], yy - pos['Mira'][1])
        bad = dmira < 6
        for _, s in g[(g['T'] < 14)].iterrows():
            bad |= np.hypot(xx - s.x, yy - s.y) < 2.5
        for _, s in g[g.source_id.isin(var_ids)].iterrows():
            bad |= np.hypot(xx - s.x, yy - s.y) < 3.0
        med = np.median(F, 0)
        thr = np.median(med[~bad])
        bpix = (~bad) & (med <= thr)
        b = np.median(F[:, bpix], axis=1)
        offs = []
        for _, s in g[(g['T'] < 11.8)].iterrows():
            if s.source_id in var_ids or s.source_id == IDS['Mira']: continue
            if not (2.5 < s.x < nx - 3.5 and 2.5 < s.y < ny - 3.5): continue
            if np.hypot(s.x - pos['Mira'][0], s.y - pos['Mira'][1]) < 5: continue
            nb = g[(np.hypot(g.x - s.x, g.y - s.y) < 2.0) & (g['T'] < s['T'] + 2) & (g.source_id != s.source_id)]
            if len(nb): continue
            ix, iy = int(round(s.x)), int(round(s.y))
            sub = F[:, iy - 2:iy + 3, ix - 2:ix + 3] - b[:, None, None]
            tot = sub.sum((1, 2))
            cx = (sub * (np.arange(-2, 3)[None, None, :])).sum((1, 2)) / tot
            cy = (sub * (np.arange(-2, 3)[None, :, None])).sum((1, 2)) / tot
            offs.append((cx - np.median(cx), cy - np.median(cy), s['T']))
        wts = np.array([10 ** (-0.4 * (o[2] - 10)) for o in offs])
        dxp = np.sum([o[0] * wi for o, wi in zip(offs, wts)], 0) / wts.sum()
        dyp = np.sum([o[1] * wi for o, wi in zip(offs, wts)], 0) / wts.sum()
        XP = np.column_stack([dxp, dyp])
        m3 = pos['M3V']; ix, iy = int(round(m3[0])), int(round(m3[1]))
        bx = np.zeros((ny, nx), bool); bx[iy - 1:iy + 2, ix - 1:ix + 2] = True
        A = F[:, bx].sum(1) - bx.sum() * b
        EE = S[str(sec)]['EE']['box3']
        n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
        for e in range(n0, n1 + 1):
            for ph in (0.0, 0.5):
                tc = T0 + (e + ph) * P
                m = np.abs(t - tc) < HALF
                if m.sum() < 8 or not accept(t[m], tc, ex): continue
                tau = t[m] - tc
                s = shape(t[m], tc, ex)
                M = np.column_stack([-s, np.ones_like(tau), tau, tau ** 2, XP[m]])
                c, *_ = np.linalg.lstsq(M, A[m], rcond=None)
                base = M[:, 1:] @ c[1:]
                r = A[m] - M @ c
                sig = 1.4826 * np.median(np.abs(r - np.median(r)))
                norm = F_M3V * EE
                cols['sector'] += [sec] * m.sum(); cols['epoch'] += [e] * m.sum(); cols['phase'] += [ph] * m.sum()
                cols['t'] += list(t[m]); cols['dt'] += list(tau); cols['f'] += list(1 + (A[m] - base) / norm)
                cols['e'] += [sig / norm] * m.sum(); cols['cad'] += [round(ex * 86400)] * m.sum(); cols['L'] += [c[0]] * m.sum()
                if ph == 0.0: chk.append((sec, e, c[0]))
        print(sec, 'EE', round(EE, 3), 'n pointing stars', len(offs), flush=True)
    os.makedirs(f'{W}/figscripts/data', exist_ok=True)
    np.savez_compressed(f'{W}/figscripts/data/rawap_box3.npz', **{k: np.array(v) for k, v in cols.items()})
    q = pd.read_csv(f'{R3}/out/q02_bigap_per_transit.csv')
    dd = [(s, e, L, float(q[(q.sector == s) & (q.epoch == e)].L_box3_pt.values[0])) for s, e, L in chk if len(q[(q.sector == s) & (q.epoch == e)])]
    diff = np.array([a[2] - a[3] for a in dd])
    print('windows', len(chk), 'matched to q02', len(dd), 'max |L - L_q02| e/s', np.abs(diff).max())


if __name__ == '__main__':
    main()
