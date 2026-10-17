"""Validation of the PSF model against the transit signal itself: in-transit difference images.

For each sector, from the p05 difference-image cube D (static scene, background and jitter response removed):
dI = mean over transits of [mean D(out of transit, 0.6 T14 < |t - tc| < 0.25 d) - mean D(flat part, |t - tc| < 0.3 T14)]
(positive = pixel dims in transit). Its per-pixel error is the scatter of the same statistic at 200 random fake
epochs. dI (pixels within 3 px of the Mira or the M3V) is fitted with a constant plus: the M3V PRF, the pair PRF,
the Mira's empirical image, and a PRF at a free position. The free-position fit gives the dimming centroid in arcsec
relative to the M3V (field-star WCS); its error from the fake-epoch scatter of the fitted position (null fits
on noise images with the real dI added would be equivalent; here we refit dI + null realisation).
Outputs: out/validate_diffimg.json, figs/diffimg_s{sector}.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p05_phot import setup, photometry
from p06_lc import BKG_EXC_MAX
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14 = TR["P"], TR["T0"], TR["T14_d"]
RNG = np.random.default_rng(3)


def stack(t, D, ok, tcs):
    imgs = []
    for tc in tcs:
        inn = ok & (np.abs(t - tc) < 0.3 * T14)
        out = ok & (np.abs(t - tc) > 0.6 * T14) & (np.abs(t - tc) < 0.25)
        if inn.sum() < 2 or out.sum() < 10: continue
        imgs.append(np.nanmean(D[out], 0) - np.nanmean(D[inn], 0))
    return (np.mean(imgs, 0), len(imgs)) if imgs else (None, 0)


def main():
    res = {}
    for s in SECTORS:
        S = setup(s); C, W, J, g, prf, rBA = S
        R = photometry(s, S=S, verbose=False, return_D=True)
        t = R["t"]; D = R["D"]; chiB = R["FIT"]["B"]["chi2"]; bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
        ok = R["qok"] & (bexc < BKG_EXC_MAX) & R["jit_ok"] & np.isfinite(chiB) & (chiB < 3 * np.nanmedian(chiB[R["q0"]]))
        tt = t[np.isfinite(t)]
        ks = np.arange(np.ceil((tt.min() - T0) / P), np.floor((tt.max() - T0) / P) + 1)
        dI, n = stack(t, D, ok, T0 + ks * P)
        if dI is None: continue
        # null realisations
        nulls = []
        real = np.abs((((np.nan_to_num(t) - T0) / P + 0.5) % 1 - 0.5) * P) < 0.15
        for k in range(200):
            ph = RNG.uniform(0.1, 0.9)
            nI, nn = stack(t, D, ok & ~real, T0 + (ks + ph) * P)
            if nI is not None and nn >= max(1, n - 1): nulls.append(nI)
        nulls = np.array(nulls); err = np.std(nulls, 0)
        yy, xx = np.mgrid[0:C.ny, 0:C.nx]
        m3 = R["pos"]["M3V"]; mi = R["pos"]["Mira"]; pa = R["pos"]["pair"]
        fm = (np.hypot(xx - mi[0], yy - mi[1]) <= 3) | (np.hypot(xx - m3[0], yy - m3[1]) <= 3)
        b = R["templates"]["Mira"]
        tm = dict(M3V=prf.image(*m3), pair=prf.image(*pa), Mira=b)
        fits_ = {}
        w = 1 / err[fm]
        for nm, T in tm.items():
            A = np.column_stack([T[fm], np.ones(fm.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(A, dI[fm] * w, rcond=None)
            fits_[nm] = dict(amp=float(c[0]), chi2=float(np.sum((A @ c - dI[fm] * w) ** 2)))
        def rr(p, img=dI):
            A = np.column_stack([prf.image(p[0], p[1])[fm], np.ones(fm.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(A, img[fm] * w, rcond=None); return A @ c - img[fm] * w
        best = None
        for gx in np.linspace(-1.5, 1.5, 7):
            for gy in np.linspace(-1.5, 1.5, 7):
                p0 = [m3[0] + gx, m3[1] + gy]; c2 = np.sum(rr(p0) ** 2)
                if best is None or c2 < best[0]: best = (c2, p0)
        q = least_squares(rr, best[1], diff_step=1e-3, x_scale=0.05)
        Jm = C.jacobian(*m3)
        off = Jm @ (q.x - np.array(m3))
        # position error: refit dI + null realisations (first 60)
        pos = []
        for nI in nulls[:60]:
            q2 = least_squares(lambda p: rr(p, dI + (nI - nulls.mean(0))), q.x, diff_step=1e-3, x_scale=0.05)
            pos.append(Jm @ (q2.x - np.array(m3)))
        pos = np.array(pos)
        rel = {k: (Jm @ (np.array(v) - np.array(m3))).tolist() for k, v in (("Mira", mi), ("pair", pa))}
        fits_["free"] = dict(chi2=float(np.sum(q.fun ** 2)), offset_from_M3V_EN=off.tolist(), offset_err=pos.std(0).tolist() if len(pos) > 5 else None)
        res[s] = dict(n_transits=n, npix=int(fm.sum()), fits=fits_, rel_positions_from_M3V=rel)
        print(f"S{s}: {n} transits, {fm.sum()} px; chi2 M3V {fits_['M3V']['chi2']:.1f} (amp {fits_['M3V']['amp']:.1f}), pair {fits_['pair']['chi2']:.1f}, "
              f"Mira {fits_['Mira']['chi2']:.1f}, free {fits_['free']['chi2']:.1f}; free offset from M3V (E,N) {np.round(off,1).tolist()} +- "
              f"{np.round(pos.std(0),1).tolist()}\"; pair at {np.round(rel['pair'],1).tolist()}, Mira at {np.round(rel['Mira'],1).tolist()}", flush=True)
        fig, ax = plt.subplots(1, 4, figsize=(16, 4))
        vm = np.nanmax(np.abs(dI[fm]))
        ax[0].imshow(np.where(fm, dI, np.nan), origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm); ax[0].set_title(f"S{s} in-transit dimming (e/s), {n} transits")
        for k, nm in enumerate(("M3V", "pair", "Mira")):
            T = tm[nm]; A = np.column_stack([T[fm], np.ones(fm.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(A, dI[fm] * w, rcond=None)
            mod = np.full(dI.shape, np.nan); mod[fm] = c[0] * T[fm] + c[1]
            ax[k + 1].imshow(np.where(fm, (dI - mod) / err, np.nan), origin="lower", cmap="RdBu_r", vmin=-4, vmax=4)
            ax[k + 1].set_title(f"resid/err, host {nm}: chi2 {fits_[nm]['chi2']:.1f}/{fm.sum()-2}")
        for a in ax:
            for (x_, y_), c_ in ((mi, "k"), (m3, "y"), (pa, "m")): a.plot(x_, y_, "+", color=c_, ms=10, mew=2)
            a.plot(q.x[0], q.x[1], "gx", ms=9, mew=2)
            a.set_xlim(mi[0] - 4.5, mi[0] + 4.5); a.set_ylim(mi[1] - 4.5, mi[1] + 4.5)
        fig.tight_layout(); fig.savefig(f"{ROOT}/figs/diffimg_s{s}.png", dpi=70); plt.close(fig)
    tot = {k: sum(v["fits"][k]["chi2"] for v in res.values()) for k in ("M3V", "pair", "Mira", "free")}
    dof = sum(v["npix"] - 2 for v in res.values())
    res["total_chi2"] = tot; res["total_dof"] = dof
    print("total chi2:", {k: round(v, 1) for k, v in tot.items()}, "dof", dof)
    save_json({str(k): v for k, v in res.items()}, f"{ROOT}/out/validate_diffimg.json")


if __name__ == "__main__":
    main()
