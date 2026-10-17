"""e01 (exploratory): residual RRab signal in the cleaned INV curve of S93 (0.78 ppt, p < 0.01): which estimator still
carries it? Variants: larger variable masks (radius x1.6), no BG0 regressor term, jitter from the 'all' mask with
second-order terms... Prints the RRab first-harmonic amplitude in M3V_B and in each regressor."""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
import cleanlib as cl
from cleanlib import *
from c01_build_inv import leak_stats
s = int(sys.argv[1]) if len(sys.argv) > 1 else 93
orig = cl.var_mask
for lab, rs in (("all_x1.0", 1.0), ("all_x1.6", 1.6), ("all_x2.2", 2.2)):
    cl.var_mask = lambda sector, mode, ny=21, nx=21, radius_scale=1.0, _rs=rs: orig(sector, mode, ny, nx, _rs)
    J = jitter_clean(s, "all", save=False)
    S = setup(s, "all"); S["J"] = J; S["vm"] = cl.var_mask(s, "all")
    R = photometry(s, S)
    L = make_lc(s, R, models=(("B", "M3V"),))
    ok = L["quality_ok"] & ~L["intransit_window"] & np.isfinite(L["det_M3V_B"])
    y = (L["det_M3V_B"][ok] - 1) * 1e3; w = 1 / (L["err_M3V_B"][ok] * 1e3) ** 2
    st = leak_stats(L["t_tdb"][ok], y, w, P_RR)
    regs = {}
    for nm, arr in (("Fap", R["Fap"]), ("BG0", R["BG"][:, 0]), ("BGx", R["BG"][:, 1]), ("BGy", R["BG"][:, 2]), ("dx", R["dx"]), ("dy", R["dy"]), ("bkg", R["bkg"])):
        g = ok & np.isfinite(arr)
        from scipy.ndimage import median_filter
        yy = arr[g] - median_filter(arr[g], size=max(5, int(1 / R["exptime"]) | 1), mode="nearest")
        regs[nm] = leak_stats(L["t_tdb"][g], yy, np.ones(g.sum()) / robust_std(yy) ** 2, P_RR)["p_amp"]
    print(f"S{s} {lab}: npix masked {S['vm'].sum()}, jitter fit pixels {J['npix_fit']}; M3V_B RR amp {st['amp1']:.2f} ppt p {st['p_amp']:.3f}; regressor p-values {regs}", flush=True)
