"""Sanity check: do the fitted difference-image amplitudes reproduce the QLP aperture depths, and
what intrinsic depth does each host hypothesis imply?

depth_ap(host) = A * sum_ap PSF_host / sum_ap M   for a circular aperture of radius 2.5 px centred on
the target (QLP '2.5:4.0:3.0' apertures; S104 used a 3x3 square), M = mean OOT image minus the
constant background fitted in the registration (the TESScut FLUX here is not background-subtracted:
FLUX_BKG = 0; floor ~180-240 e/s/px of sky + unresolved stars).  QLP's own background treatment in
this crowded field is not reproduced, so only rough agreement is expected.
Intrinsic depth on the host = A / F_host with F_host from the registration fit (free fluxes) and,
alternatively, from the Gaia-predicted T-mag ratio (target flux from the fit).
Outputs: depth_consistency.json
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pipeline import *

real = json.load(open(f'{HERE}/diff_prf_main.json'))
reg = json.load(open(f'{HERE}/registration.json'))
qlp = {13: 4101, 39: 4238, 66: 3792, 93: 3629, 104: 2939}
out = {}
for sec in SECTORS:
    ctx = SectorContext(sec, psfkind='prf', reg_variant='local_r4')
    S = ctx.S
    M = S.flux[S.good].mean(0)
    st = real[str(sec)]['stacked']
    A = st['free_amp']
    aT, aN = st['a_T'], st['a_N']
    if sec == 104:
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]
        cx, cy = np.round(ctx.pT).astype(int)
        ap = (np.abs(xx - cx) <= 1) & (np.abs(yy - cy) <= 1)
    else:
        ap = circle_mask(M.shape, ctx.pT[0], ctx.pT[1], 2.5)
    bkg = reg[str(sec)]['variants']['prf:local_r4']['bkg']   # fitted constant (sky + unresolved stars), e/s/px
    Fap = (M[ap] - bkg).sum()
    dT = A * ctx.psf(*ctx.pT)[ap].sum() / Fap
    dN = A * ctx.psf(*ctx.pN)[ap].sum() / Fap
    d2 = (aT * ctx.psf(*ctx.pT)[ap].sum() + aN * ctx.psf(*ctx.pN)[ap].sum()) / Fap
    v = reg[str(sec)]['variants']['prf:local_r4']
    FT, FN = v['F_target'], v['F_neigh']
    ratio_cat = reg[str(sec)]['neigh_target_ratio_catalog']
    out[str(sec)] = dict(A=A, F_ap=float(Fap), bkg_fit=bkg, depth_ap_if_target=float(dT), depth_ap_if_neigh=float(dN),
                         depth_ap_2source=float(d2), qlp_depth=qlp[sec] * 1e-6,
                         F_target_fit=FT, F_neigh_fit=FN, F_neigh_catalog_ratio=FT * ratio_cat,
                         intrinsic_depth_target=A / FT, intrinsic_depth_neigh_fitflux=A / FN,
                         intrinsic_depth_neigh_catflux=A / (FT * ratio_cat))
    o = out[str(sec)]
    print(f"S{sec}: A={A:.1f} e/s  aperture depth: if target {dT*1e6:.0f} ppm, if neighbour {dN*1e6:.0f} ppm, 2-source {d2*1e6:.0f} ppm; QLP {qlp[sec]} ppm | "
          f"intrinsic: target {o['intrinsic_depth_target']*1e6:.0f} ppm; neighbour {o['intrinsic_depth_neigh_fitflux']*100:.2f}% (fit flux {FN:.0f}) / {o['intrinsic_depth_neigh_catflux']*100:.2f}% (catalogue ratio)")
to_json(out, f'{HERE}/depth_consistency.json')
