"""v4 VHS: joint fit (primary shift + companion position + linear amplitudes) of the ~1.1" PA ~347 feature, with a
flux/position calibration by injecting C2-shaped companions of known f into T at the same separation and 7 other
PAs (one per 45 deg, starting 35 deg from the feature), fitted identically (template C1, equally saturated).
Also the significance: dchi2(joint companion) at the feature vs the distribution of the best joint dchi2 obtained
when starting the fit at the 7 other PAs (local maxima of noise / PSF mismatch at the same separation)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h01_vhs_loop import *

pp = [f for f in FILES if 'big' in f] + [f for f in FILES if 'big' not in f]
FREE = os.environ.get('H03_FREE', '0') == '1'
SEP0 = float(os.environ.get('H03_SEP', 1.1))


def joint(S, P, e0, n0, J, Ji, inj=None, free_pos=None):
    free_pos = FREE if free_pos is None else free_pos
    p0, c2_0, c0, npx, solve = fit_pair(S, P, inj=inj)
    def f(q):
        off = Ji @ np.array([q[2], q[3]]) if free_pos else Ji @ np.array([e0, n0])
        c2, cc, _ = solve(q[:2], comp_off=(off[0], off[1]))
        return c2
    q0 = np.r_[p0, e0, n0] if free_pos else p0
    if free_pos:
        o = minimize(f, q0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3, maxiter=600))
    else:
        o = minimize(lambda q: f(np.r_[q, e0, n0]), p0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3, maxiter=300))
        o.x = np.r_[o.x, e0, n0]
    off = Ji @ o.x[2:4]
    c2, cc, _ = solve(o.x[:2], comp_off=(off[0], off[1]))
    return dict(dchi2=float(c2_0 - c2), f=float(cc[2] / cc[0]), e=float(o.x[2]), n=float(o.x[3]), sep=float(np.hypot(*o.x[2:4])),
                pa=float(np.degrees(np.arctan2(o.x[2], o.x[3])) % 360), chi2r0=float(c2_0 / (npx - 4)))


out = []
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
        Ji = np.linalg.inv(J) * 1.0 / PIX * PIX      # arcsec -> px
        # template choice: C1 (equally saturated) when available in the stamp, else C2
        tmpl = 'C1' if st['C1']['raw'].max() > 0.5 * st['T']['raw'].max() else 'C2'
        e0, n0 = SEP0 * np.sin(np.radians(347)), SEP0 * np.cos(np.radians(347))
        feat = joint(st['T'], st[tmpl], e0, n0, J, Ji)
        ctrl = []
        for pa in 347 + 35 + 45 * np.arange(7):
            e1, n1 = SEP0 * np.sin(np.radians(pa)), SEP0 * np.cos(np.radians(pa))
            ctrl.append(joint(st['T'], st[tmpl], e1, n1, J, Ji))
        # same for comparison star C1 with template C2 (and vice versa) at the feature position: does C1 show it?
        cmp_ = joint(st['C1'], st['C2'], e0, n0, J, Ji) if tmpl == 'C1' else None
        # calibration: inject at the 7 control PAs with f_in, refit
        p2, _, cc2, _, _ = fit_pair(st['T'], st['C2'])
        cal = []
        for f_in in (0.01, 0.02, 0.04):
            for pa in 347 + 35 + 45 * np.arange(7):
                e1, n1 = SEP0 * np.sin(np.radians(pa)), SEP0 * np.cos(np.radians(pa))
                pxo = Ji @ np.array([e1, n1])
                comp = cc2[0] * f_in * ndshift(st['C2']['s'], (p2[1] + pxo[1], p2[0] + pxo[0]), order=3, mode='constant', cval=0.0)
                r = joint(st['T'], st[tmpl], e1, n1, J, Ji, inj=comp)
                cal.append(dict(f_in=f_in, pa_in=float(pa % 360), f_out=r['f'], dpos=float(np.hypot(r['e'] - e1, r['n'] - n1)), dchi2=r['dchi2']))
        rat = {f_in: float(np.median([c['f_out'] / c['f_in'] for c in cal if c['f_in'] == f_in])) for f_in in (0.01, 0.02, 0.04)}
        # invert: find f_in whose median recovered f matches the observed f
        fin = np.array(sorted(rat)); fout = np.array([rat[k] * k for k in fin])
        f_true = float(np.interp(feat['f'], fout, fin)) if fout[0] < feat['f'] < fout[-1] else float(feat['f'] / rat[0.02])
        rowc = dict(file=os.path.basename(fn), ext=ext, band=band, template=tmpl, feature=feat, controls=ctrl, comparison_C1_C2=cmp_, calib=cal,
                    median_ratio=rat, f_true_est=f_true, dm_obs=float(-2.5 * np.log10(max(feat['f'], 1e-6))), dm_true_est=float(-2.5 * np.log10(max(f_true, 1e-6))))
        out.append(rowc)
        cd2 = sorted([c['dchi2'] for c in ctrl])
        print(f"{os.path.basename(fn)[4:30]} e{ext} {band} [{tmpl}]: feature sep {feat['sep']:.2f}\" PA {feat['pa']:.0f} f {feat['f']:.4f} dchi2 {feat['dchi2']:.0f} "
              f"(chi2r0 {feat['chi2r0']:.1f}); controls dchi2 {np.round(cd2, 0)} f {np.round([c['f'] for c in ctrl], 4)}", flush=True)
        if cmp_:
            print(f"     C1|C2 at feature: sep {cmp_['sep']:.2f} PA {cmp_['pa']:.0f} f {cmp_['f']:.4f} dchi2 {cmp_['dchi2']:.0f}", flush=True)
        print(f"     calib median f_out/f_in {rat};  dpos med {np.median([c['dpos'] for c in cal]):.2f}\";  dm_obs {rowc['dm_obs']:.2f} -> dm_true ~ {rowc['dm_true_est']:.2f}", flush=True)
to_json(out, f'{V4}/results/h03_vhs_joint_sep{SEP0}_free{int(FREE)}.json')
