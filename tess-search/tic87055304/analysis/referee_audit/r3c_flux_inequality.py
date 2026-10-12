"""Zero-PSF-model inequality (S93 and S104, the sectors where the neighbour has its own pixel n != target pixel t).

For ANY PSF that is symmetric and unimodal, the pixel containing a star receives at least as much of that star's light
as an edge-adjacent pixel in the same row/column (1-D argument per axis; t and n share a row in S104 and a column in
S93). If the TARGET hosts the transit, dF(p) = delta_T * T(p), hence dF(n) <= dF(t).
Observed dF(n) > dF(t) therefore rejects a target host using only (i) which pixel contains the target (WCS; margin =
distance of the target from the t/n pixel boundary) and (ii) PSF symmetry/unimodality.
Noise: same fake-transit stacks as r3 (phases 0.20-0.80), z = (dF(n) - dF(t)) / std_null.
"""
import numpy as np
from rcommon import Cube, OUT, dump
from r3_pixel_ratio import collect, clean, stack, PHASES, NBOOT

RNG = np.random.default_rng(99)
res = {}
for s in (93, 104):
    c = Cube(s)
    tp = (int(round(c.xt)), int(round(c.yt))); npx = (int(round(c.xn)), int(round(c.yn)))
    real = collect(c); fake = [w for ph in PHASES for w in collect(c, ph)]
    real, fake, _ = clean(real, fake)          # R3_CLEAN default 'mad10' unless overridden
    for w in real + fake:
        w['wt'] = 1.0 / w['sD'][tp[1], tp[0]] ** 2
    D, O, B = stack(real, range(len(real)))
    diff = D[npx[1], npx[0]] - D[tp[1], tp[0]]
    nd = []
    for _ in range(NBOOT):
        Df, _, _ = stack(fake, RNG.choice(len(fake), len(real), replace=False))
        nd.append(Df[npx[1], npx[0]] - Df[tp[1], tp[0]])
    nd = np.array(nd)
    axis = 0 if npx[1] == tp[1] else 1
    pos = c.xt if axis == 0 else c.yt
    boundary = (tp[axis] + npx[axis]) / 2.0
    pt = []
    for k in range(len(real)):
        Dk, _, _ = stack(real, [k])
        pt.append(float(Dk[npx[1], npx[0]] - Dk[tp[1], tp[0]]))
    res[str(s)] = dict(pix_t=tp, pix_n=npx, dF_n=float(D[npx[1], npx[0]]), dF_t=float(D[tp[1], tp[0]]),
                       O_n=float(O[npx[1], npx[0]]), O_t=float(O[tp[1], tp[0]]), diff=float(diff),
                       sigma_null=float(nd.std()), z=float(diff / nd.std()), frac_null_ge=float(np.mean(nd >= diff)),
                       target_distance_from_pixel_boundary_px=float(abs(pos - boundary)), per_transit_diff=pt,
                       n_transits=len(real))
    print(f"S{s}: dF(n)={D[npx[1], npx[0]]:.2f} dF(t)={D[tp[1], tp[0]]:.2f} e/s (O_n {O[npx[1], npx[0]]:.0f}, O_t "
          f"{O[tp[1], tp[0]]:.0f}); diff {diff:.2f} +- {nd.std():.2f} -> z = {diff/nd.std():.1f}; per transit {np.round(pt,2)}; "
          f"target {abs(pos-boundary):.3f} px from the t/n boundary")
zs = [res[k]['z'] for k in res]
res['stouffer'] = float(np.sum(zs) / np.sqrt(len(zs)))
print('Stouffer', round(res['stouffer'], 2))
dump(res, f'{OUT}/r3c_flux_inequality.json')
