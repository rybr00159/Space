"""v4 VHS: residual images of the closed-loop template subtraction + flux-bias calibration of the 1.1" feature.
Residual R(S|P) = S - a*shift(P) - b (masked pixels shown grey) for T|C1, C1|T, T|C2, C1|C2 in each pawprint.
Flux bias: a fake companion (C2 shape) with known f is injected into T at the MIRROR position (1.1", PA 170) and its
flux recovered with template C1 exactly as for the feature: recovered/injected ratio -> corrected contrast."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h01_vhs_loop import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

pp = [f for f in FILES if 'big' in f]
fig, axs = plt.subplots(len(pp), 4, figsize=(14, 3.6 * len(pp)))
outb = []
for row, fn in enumerate(pp):
    hd0 = fits.getheader(fn, 0)
    ext = 2 if len(fits.open(fn)) > 2 else 1
    d = fits.getdata(fn, ext).astype(float); hd = fits.getheader(fn, ext); w = WCS(hd)
    band = hd0.get('FILTER')
    st = {k: stamp(d, w, k) for k in IDS}
    x0, y0 = w.all_world2pix([POS['T']], 0)[0]
    r0 = w.all_pix2world([[x0, y0], [x0 + 1, y0], [x0, y0 + 1]], 0); cd = np.cos(np.radians(r0[0][1]))
    J = np.column_stack([[(r0[1][0] - r0[0][0]) * cd * 3600, (r0[1][1] - r0[0][1]) * 3600], [(r0[2][0] - r0[0][0]) * cd * 3600, (r0[2][1] - r0[0][1]) * 3600]])
    Ji = np.linalg.inv(J)
    for col, (a, b) in enumerate([('T', 'C1'), ('C1', 'T'), ('T', 'C2'), ('C1', 'C2')]):
        p, c2, c, npx, solve = fit_pair(st[a], st[b])
        m = ndshift(st[b]['s'], (p[1], p[0]), order=3, mode='constant', cval=0.0)
        satm = ndshift((st[b]['raw'] > SAT).astype(float), (p[1], p[0]), order=1, mode='constant', cval=0.0) > 0.01
        mask = binary_dilation((st[a]['raw'] > SAT) | satm, iterations=1)
        R = st[a]['s'] - c[0] * m - c[1]
        sig = np.sqrt(st[a]['rms'] ** 2 + np.clip(st[a]['s'], 0, None) / GAIN + (FLOOR * c[0] * m) ** 2)
        Z = np.where(mask, np.nan, R / sig)
        ax = axs[row, col]
        ax.imshow(Z, origin='lower', cmap='RdBu_r', vmin=-6, vmax=6, extent=[-H - 0.5, H + 0.5, -H - 0.5, H + 0.5])
        # mark 1.12" PA 350 and N/E arrows
        e_, n_ = 1.12 * np.sin(np.radians(350)), 1.12 * np.cos(np.radians(350))
        px = Ji @ np.array([e_, n_])
        ax.plot(px[0] + st[a]['fx'], px[1] + st[a]['fy'], 'o', mfc='none', mec='k', ms=14)
        nvec = Ji @ np.array([0, 2.0]); evec = Ji @ np.array([2.0, 0])
        ax.annotate('', xy=(-10 + nvec[0], -10 + nvec[1]), xytext=(-10, -10), arrowprops=dict(arrowstyle='->'))
        ax.annotate('', xy=(-10 + evec[0], -10 + evec[1]), xytext=(-10, -10), arrowprops=dict(arrowstyle='->', color='g'))
        ax.set_xlim(-9, 9); ax.set_ylim(-9, 9)
        ax.set_title(f"{os.path.basename(fn)[15:27]} {band} {a}|{b}  (resid/sigma)")
    # flux bias at the mirror position
    p2, _, cc2, _, _ = fit_pair(st['T'], st['C2'])
    rec = []
    for f_in in (0.01, 0.02, 0.04):
        for pa in (170.0, 80.0, 260.0):
            e_, n_ = 1.12 * np.sin(np.radians(pa)), 1.12 * np.cos(np.radians(pa))
            pxo = Ji @ np.array([e_, n_])
            comp = cc2[0] * f_in * ndshift(st['C2']['s'], (p2[1] + pxo[1], p2[0] + pxo[0]), order=3, mode='constant', cval=0.0)
            p_, c0, c_, npx_, solve_ = fit_pair(st['T'], st['C1'], inj=comp)
            _, cc, _ = solve_(p_, comp_off=(pxo[0], pxo[1]))
            rec.append(dict(f_in=f_in, pa=pa, f_out=float(cc[2] / cc[0])))
    ratio = float(np.median([r['f_out'] / r['f_in'] for r in rec]))
    # the feature itself with template C1
    e_, n_ = 1.12 * np.sin(np.radians(350)), 1.12 * np.cos(np.radians(350)); pxo = Ji @ np.array([e_, n_])
    p_, c0, c_, npx_, solve_ = fit_pair(st['T'], st['C1'])
    _, cc, _ = solve_(p_, comp_off=(pxo[0], pxo[1]))
    f_obs = float(cc[2] / cc[0])
    outb.append(dict(file=os.path.basename(fn), band=band, f_obs=f_obs, recovery_ratio=ratio, f_corr=f_obs / ratio,
                     dm_obs=float(-2.5 * np.log10(max(f_obs, 1e-6))), dm_corr=float(-2.5 * np.log10(max(f_obs / ratio, 1e-6))), trials=rec))
    print(f"{os.path.basename(fn)[4:30]} {band}: feature f_obs {f_obs:.4f} (dm {-2.5 * np.log10(max(f_obs, 1e-6)):.2f}); injected-flux recovery ratio {ratio:.2f} "
          f"-> corrected dm {-2.5 * np.log10(max(f_obs / ratio, 1e-6)):.2f}; trials {[round(r['f_out'] / r['f_in'], 2) for r in rec]}", flush=True)
plt.tight_layout(); plt.savefig(f'{V4}/plots/h02_vhs_residuals.png', dpi=70)
to_json(outb, f'{V4}/results/h02_vhs_fluxbias.json')
