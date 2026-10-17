"""c02: rebuild the referee ('REF') local-window photometry with the field variables excluded from B(t) and the
pointing (refclean.aux_clean), and measure the RRab / EB leak before and after on a contiguous ('tiled') REF-method
light curve (windows every 0.3 d, central +-0.15 d kept).

Modes: orig (the referee's own aux file), none (aux recomputed by this code without masks: must equal orig),
rr3, all (ADOPTED).
Per mode and sector: strict real windows at the phase-1 ephemeris, 150 fake windows (|phase| > 0.45 d from transit),
the tiled light curve, and leak statistics (same as c01: 3 harmonics at the Gaia RRab period, 4 at the EB periods;
120 fake periods). Outputs: data/refwin_{mode}_s{sector}.npz, data/reftile_{mode}_s{sector}.npz, out/leak_ref.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from refclean import *
from cleanlib import P_RR, P_EB075, P_EB159, P1, T01, F_M3V
from c01_build_inv import leak_stats
import r05_local as R5

MODES = ("orig", "none", "rr3", "all")


def run(s):
    res = {}
    rng = np.random.default_rng(2026 + s)
    for mode in MODES:
        S = SectorC(s, mode)
        if mode == "none":
            o = np.load(f"{RFR}/data/aux_s{s}.npz")
            res["none_vs_orig_maxdiff"] = {k: float(np.nanmax(np.abs(S.a[k] - o[k]))) for k in ("B", "dx", "dy")}
        # real windows at the phase-1 ephemeris
        cen, ks = R5.real_centres(S, P=P1, T0=T01)
        out = {}
        for lab, bmax in (("strict", 60.0), ("relaxed", 800.0)):
            rr = []; kk = []
            for tc, k in zip(cen, ks):
                r = S.window(tc, bmax=bmax)
                if r is not None: rr.append(r); kk.append(k)
            for key, v in R5.pack(rr).items(): out[f"{lab}_{key}"] = v
            out[f"{lab}_epoch"] = np.array(kk)
        fk = R5.fake_centres(S, 150, rng, P=P1, T0=T01)
        for key, v in R5.pack(R5.windows_for(S, fk, bmax=60.0)).items(): out[f"fake_{key}"] = v
        np.savez_compressed(f"{CL}/data/refwin_{mode}_s{s}.npz", **out, exptime=S.exptime, pos=json.dumps(S.pos))
        # tiled light curve and leak statistics
        t, a, sg, w = tiled(S)
        ph = ((t - T01) / P1 + 0.5) % 1 - 0.5
        k = np.abs(ph * P1) > 0.1
        np.savez_compressed(f"{CL}/data/reftile_{mode}_s{s}.npz", t=t, aM=a, sig=sg, win=w)
        y = a[k] / F_M3V * 1e3; ww = 1 / (sg[k] / F_M3V * 1e3) ** 2
        st = dict(RR=leak_stats(t[k], y, ww, P_RR), EB075=leak_stats(t[k], y, ww, P_EB075, nh=4),
                  EB159=leak_stats(t[k], y, ww, P_EB159, nh=4), rstd_ppt=float(rstd(y)), n=int(k.sum()))
        for nm in ("B", "dx", "dy"):
            g = S.qok & np.isfinite(S.a[nm])
            from scipy.ndimage import median_filter
            yy = S.a[nm][g] - np.nanmedian(S.a[nm][g]); n = max(5, int(1.0 / S.exptime) | 1)
            yy = yy - median_filter(yy, size=n, mode="nearest")
            st[nm] = dict(RR=leak_stats(S.ttdb[g], yy, np.ones(g.sum()) / max(rstd(yy), 1e-9) ** 2, P_RR))
        res[mode] = st
        print(f"S{s} {mode:5s}: tiled REF RR amp1 {st['RR']['amp1']:.2f} ppt (proj {st['RR']['proj']:+.2f}, p {st['RR']['p_amp']:.3f}, "
              f"UL95 {st['RR']['ul95_amp']:.2f}); EB075 dchi2 {st['EB075']['dchi2']:.1f} (p {st['EB075']['p_dchi2']:.2f}); EB159 "
              f"{st['EB159']['dchi2']:.1f} (p {st['EB159']['p_dchi2']:.2f}); dx RR p {st['dx']['RR']['p_amp']:.2f}, B RR p {st['B']['RR']['p_amp']:.2f}; "
              f"rstd {st['rstd_ppt']:.2f} ppt; strict epochs {out['strict_epoch'].tolist()}", flush=True)
    return res


if __name__ == "__main__":
    secs = [int(a) for a in sys.argv[1:]] or list(SECTORS)
    fn = f"{CL}/out/leak_ref.json"
    allr = json.load(open(fn)) if os.path.exists(fn) else {}
    for s in secs:
        allr[str(s)] = run(s)
        jd(allr, fn)
