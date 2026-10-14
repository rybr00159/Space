"""v4 VHS: polar maps (sep 0.8-2.3", PA step 15 deg) of the joint-fit companion dchi2 (primary shift re-optimised at every
grid point, companion flux free with either sign) for T|C1 (pawprints) or T|C2 (tiles), and the null C1|C2.
Summed over frames per band, the peak location and its significance relative to the largest null peak."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h01_vhs_loop import *

SEPS = np.arange(0.8, 2.31, 0.15)
PAS = np.arange(0, 360, 15.0)
pp = [f for f in FILES if 'big' in f] + [f for f in FILES if 'big' not in f]


def jmap(S, P, Ji):
    p0, c2_0, c0, npx, solve = fit_pair(S, P)
    D = np.zeros((len(SEPS), len(PAS))); F = np.zeros_like(D)
    for i, s in enumerate(SEPS):
        for j, pa in enumerate(PAS):
            off = Ji @ np.array([s * np.sin(np.radians(pa)), s * np.cos(np.radians(pa))])
            o = minimize(lambda q: solve(q, comp_off=(off[0], off[1]))[0], p0, method='Nelder-Mead', options=dict(xatol=2e-3, fatol=1e-2, maxiter=150))
            c2, cc, _ = solve(o.x, comp_off=(off[0], off[1]))
            D[i, j] = c2_0 - c2; F[i, j] = cc[2] / cc[0]
    return D, F, c2_0 / (npx - 4)


res = []
for fn in pp:
    hd0 = fits.getheader(fn, 0)
    for ext in range(1, len(fits.open(fn))):
        d = fits.getdata(fn, ext).astype(float); hd = fits.getheader(fn, ext)
        if min(d.shape) < 200:
            continue
        w = WCS(hd); band = hd0.get('FILTER') or hd.get('FILTER')
        st = {k: stamp(d, w, k) for k in IDS}
        x0, y0 = w.all_world2pix([POS['T']], 0)[0]
        r0 = w.all_pix2world([[x0, y0], [x0 + 1, y0], [x0, y0 + 1]], 0); cd = np.cos(np.radians(r0[0][1]))
        J = np.column_stack([[(r0[1][0] - r0[0][0]) * cd * 3600, (r0[1][1] - r0[0][1]) * 3600], [(r0[2][0] - r0[0][0]) * cd * 3600, (r0[2][1] - r0[0][1]) * 3600]])
        Ji = np.linalg.inv(J)
        tm = 'C1' if st['C1']['raw'].max() > 0.5 * st['T']['raw'].max() else 'C2'
        DT, FT, crT = jmap(st['T'], st[tm], Ji)
        row = dict(file=os.path.basename(fn), ext=ext, band=band, template=tm, D_T=DT, F_T=FT, chi2r_T=crT)
        if tm == 'C1':
            DN, FN, crN = jmap(st['C1'], st['C2'], Ji)
            row.update(D_null=DN, F_null=FN, chi2r_null=crN)
        i, j = np.unravel_index(np.argmax(np.where(FT > 0, DT, -np.inf)), DT.shape)
        msg = f"{os.path.basename(fn)[4:30]} e{ext} {band} [{tm}]: T peak dchi2 {DT[i, j]:.0f} at {SEPS[i]:.2f}\" PA {PAS[j]:.0f} f {FT[i, j]:.4f}"
        if tm == 'C1':
            k, l = np.unravel_index(np.argmax(np.where(FN > 0, DN, -np.inf)), DN.shape)
            msg += f" | null C1|C2 peak dchi2 {DN[k, l]:.0f} at {SEPS[k]:.2f}\" PA {PAS[l]:.0f} f {FN[k, l]:.4f}"
        print(msg, flush=True)
        res.append(row)
for band in ('J', 'Ks'):
    rows = [r for r in res if r['band'] == band and 'big' not in r['file'] or (r['band'] == band)]
    S = sum(r['D_T'] / r['chi2r_T'] for r in res if r['band'] == band)
    i, j = np.unravel_index(np.argmax(S), S.shape)
    prof_pa = S.max(0)
    print(f"{band}: summed (dchi2/chi2r) map peak {S[i, j]:.0f} at {SEPS[i]:.2f}\" PA {PAS[j]:.0f}; PA profile (max over sep) "
          + ' '.join(f"{int(p)}:{v:.0f}" for p, v in zip(PAS, prof_pa)))
    sep_prof = S[:, j]
    print(f"    sep profile at PA {PAS[j]:.0f}: " + ' '.join(f"{s:.2f}:{v:.0f}" for s, v in zip(SEPS, sep_prof)))
to_json(dict(seps=SEPS, pas=PAS, frames=res), f'{V4}/results/h04_vhs_map.json')
