"""c01: rebuild the INV (phot model B / A) photometry with the field variables excluded from every field-wide
estimator, and measure the leak of the RRab and the two Gaia EBs before and after.

Modes (cleanlib.var_mask): none (= phase-1 phot, reproduced to 0.0 e-/s), rr3 (RRab + 2 EBs masked), all (every Gaia
variable with G < 18.6 masked; ADOPTED), allreg (all + aperture light curves of the RRab and both EBs as explicit
regressors in the pixel regression).

Leak statistics per sector and light curve (OOT, quality_ok, detrended as phot p06):
  RRab: first-harmonic semi-amplitude of a 3-harmonic fit at the Gaia period (phase referred to the Gaia epoch of
        maximum); projection on the RRab-maximum phase (negative = inverted); p_emp from 200 fake periods 0.45-0.95 d
        (avoiding P_RR, P_RR/2, the EB periods and their halves within 2%); 95% upper limit = amplitude whose
        fake-period distribution of (measured - injected) ... approximated as amp + 1.645 x rms(fake amplitudes)/1.25
        (Rayleigh rms -> sigma); reported both.
  EBs:  4-harmonic delta-chi2 at the EB period vs the same fake periods.
Regressors dx, dy, bkg, BG0: the same statistics (in their own units).
Outputs: out/leak_inv.json, data/lcinv_{mode}_s{sector}.npz
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *

rng = np.random.default_rng(31)
MODES = ("none", "rr3", "all", "allreg")


def fake_periods(n=200):
    bad = [P_RR, P_RR / 2, P_EB075, P_EB075 / 2, P_EB159 / 2, P_EB159 / 3]
    out = []
    while len(out) < n:
        p = rng.uniform(0.45, 0.95)
        if all(abs(p / b - 1) > 0.02 for b in bad): out.append(p)
    return np.array(out)


FAKE = fake_periods()


def leak_stats(t, y, w, Pp, nh=3, T0=TMAX_RR):
    a1, proj, c = harm_fit(t, y, w, Pp, T0, nh)
    X0 = np.column_stack([np.ones(len(t))]) * np.sqrt(w)[:, None]
    def dchi(P_):
        _, _, cc = harm_fit(t, y, w, P_, T0, nh)
        ph = 2 * np.pi * (t - T0) / P_
        X = np.column_stack([np.ones(len(t))] + [f(h * ph) for h in range(1, nh + 1) for f in (np.cos, np.sin)])
        r1 = (y - X @ cc) * np.sqrt(w); r0 = (y - np.average(y, weights=w)) * np.sqrt(w)
        return float(np.sum(r0 ** 2) - np.sum(r1 ** 2))
    d_real = dchi(Pp)
    fa = []; fd = []
    for P_ in FAKE[:120]:
        a_, _, _ = harm_fit(t, y, w, P_, T0, nh); fa.append(a_)
        fd.append(dchi(P_))
    fa = np.array(fa); fd = np.array(fd)
    sig1 = np.sqrt(np.mean(fa ** 2) / 2)          # per-quadrature sigma from Rayleigh-distributed noise amplitudes
    return dict(amp1=a1, proj=proj, dchi2=d_real, p_amp=float(np.mean(fa >= a1)), p_dchi2=float(np.mean(fd >= d_real)),
                fake_amp_median=float(np.median(fa)), fake_amp_p95=float(np.percentile(fa, 95)), sigma_quad=float(sig1),
                ul95_amp=float(np.sqrt(max(a1 ** 2 - 2 * sig1 ** 2, 0)) + 1.645 * sig1 * 1.2))


def run(sector):
    res = {}
    for mode in MODES:
        mjit = "all" if mode == "allreg" else mode
        S = setup(sector, mjit)
        extra = None
        if mode == "allreg":
            extra = {}
            for nm, sid in (("RR", RR_ID), ("EB075", EB075_ID), ("EB159", EB159_ID)):
                v = var_aperture_lc(S["C"], sector, sid)
                if v is not None and np.isfinite(v).sum() > 100:
                    v = np.where(np.isfinite(v), v, np.nanmedian(v)); extra[nm] = v
        R = photometry(sector, S, extra_regs=extra)
        L = make_lc(sector, R)
        ok = L["quality_ok"] & ~L["intransit_window"] & np.isfinite(L["det_M3V_B"])
        t = L["t_tdb"][ok]
        st = {}
        for key in ("M3V_B", "M3V_A", "Mira_B", "pair_A"):
            y = (L[f"det_{key}"][ok] - 1) * 1e3; w = 1 / (L[f"err_{key}"][ok] * 1e3) ** 2
            k = np.isfinite(y) & np.isfinite(w)
            st[key] = dict(RR=leak_stats(t[k], y[k], w[k], P_RR), EB075=leak_stats(t[k], y[k], w[k], P_EB075, nh=4),
                           EB159=leak_stats(t[k], y[k], w[k], P_EB159, nh=4), rstd_ppt=float(robust_std(y[k])))
        g = L["quality_ok"] & np.isfinite(R["dx"]) & np.isfinite(R["BG"][:, 0])
        for nm, arr in (("dx", R["dx"]), ("dy", R["dy"]), ("bkg", R["bkg"]), ("BG0", R["BG"][:, 0]), ("BGy", R["BG"][:, 2])):
            y = arr[g] - np.nanmedian(arr[g])
            # high-pass (1-d running median) so that slow drifts do not dominate
            from scipy.ndimage import median_filter
            n = max(5, int(1.0 / R["exptime"]) | 1)
            y = y - median_filter(y, size=n, mode="nearest")
            w = np.ones(g.sum()) / max(robust_std(y), 1e-9) ** 2
            st[nm] = dict(RR=leak_stats(L["t_tdb"][g], y, w, P_RR), EB075=leak_stats(L["t_tdb"][g], y, w, P_EB075, nh=4))
        res[mode] = st
        keep = {k: v for k, v in L.items() if isinstance(v, np.ndarray) and (k.startswith(("det_", "raw_", "err_", "dF_", "trend_")) or k in
                ("t_tdb", "t_tess", "q", "quality_ok", "bkg", "dx", "dy", "intransit_window"))}
        np.savez_compressed(f"{CL}/data/lcinv_{mode}_s{sector}.npz", exptime=R["exptime"], BG0=R["BG"][:, 0], Fap=R["Fap"],
                            npix_jit=S["J"].get("npix_fit", -1), **keep)
        print(f"S{sector} {mode:6s}: M3V_B RR amp1 {st['M3V_B']['RR']['amp1']:.2f} ppt (proj {st['M3V_B']['RR']['proj']:+.2f}, "
              f"p {st['M3V_B']['RR']['p_amp']:.3f}, UL95 {st['M3V_B']['RR']['ul95_amp']:.2f}); EB075 dchi2 {st['M3V_B']['EB075']['dchi2']:.1f} "
              f"(p {st['M3V_B']['EB075']['p_dchi2']:.2f}); EB159 {st['M3V_B']['EB159']['dchi2']:.1f} (p {st['M3V_B']['EB159']['p_dchi2']:.2f}); "
              f"dx RR amp {st['dx']['RR']['amp1']:.2e} (p {st['dx']['RR']['p_amp']:.2f}); BG0 RR {st['BG0']['RR']['amp1']:.3f} (p {st['BG0']['RR']['p_amp']:.2f}); "
              f"rstd {st['M3V_B']['rstd_ppt']:.2f} ppt", flush=True)
    return res


if __name__ == "__main__":
    secs = [int(a) for a in sys.argv[1:]] or list(SECTORS)
    allr = {}
    fn = f"{CL}/out/leak_inv.json"
    if os.path.exists(fn): allr = json.load(open(fn))
    for s in secs:
        allr[str(s)] = run(s)
        jd(allr, fn)
