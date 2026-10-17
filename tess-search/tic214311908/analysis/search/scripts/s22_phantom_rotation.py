"""s22: is the recurring ~9-10-d (and ~4.5-7-d) power in the deblended M3V light curve stellar rotation or a systematic
shared by every PSF fit next to the Mira?

Per sector, on the cached difference images (s04 crops): per-cadence weighted fits of
  model B'  : Mira empirical image + M3V PRF + plane                        -> 'M3V'
  phantom k : Mira empirical image + M3V PRF + PRF at an empty position k + plane -> 'phantom_k' (10 positions as s03:
              5 on the M3V-Mira ring, 5 at 1.8 px from the M3V)
Undetrended fluxes (e-/s), quality_ok cadences outside the 5.69-d windows, regressed on background excess and pointing,
RRab leak (0.648528 d, 3 harmonics) removed, 1-h bins. For every series: LS periodogram 1-14 d (top peak) and the
sinusoid (+ quadratic trend) semi-amplitude and phase at the M3V's top period of that sector.
If the M3V's peak were stellar it should stand out against the phantoms' amplitudes at the same period.
Output: out/phantom_rotation.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/search/scripts")
sys.path.insert(1, "tic214311908/analysis/phot/scripts")
from slib import *
from s03_controls import phantoms
from s07_periodic_model import harm_X, harm_fit
from s14_rotation import ls_top, sin_vec
from astropy.timeseries import LombScargle


def main():
    import photlib
    out = {}
    for s in SECTORS:
        z = np.load(f"{ROOT}/data/pix_s{s}.npz")
        t = z["t_tdb"]; ok = z["ok"].astype(bool) & ~z["intr"].astype(bool) & np.isfinite(t)
        D = z["D"].astype(float); ok &= np.all(np.isfinite(D.reshape(len(t), -1)), 1)
        names = list(z["prf_names"]); b = z["b"]; T3 = z["prf_imgs"][names.index("M3V")]
        x0, y0 = int(z["x0"]), int(z["y0"]); ny, nx = D.shape[1:]
        xym = z["pos_xy"][list(z["pos_names"]).index("Mira")]; xy3 = z["pos_xy"][list(z["pos_names"]).index("M3V")]
        yy, xx = np.mgrid[0:ny, 0:nx]
        fm = ((np.hypot(xx - xym[0], yy - xym[1]) <= 3) | (np.hypot(xx - xy3[0], yy - xy3[1]) <= 3)).ravel()
        Dg = D[ok].reshape(ok.sum(), -1)[:, fm]
        sig = 1.4826 * np.nanmedian(np.abs(np.diff(Dg, axis=0)), axis=0) / np.sqrt(2) + 0.5
        u = ((xx - xym[0]) / 3).ravel()[fm]; v = ((yy - xym[1]) / 3).ravel()[fm]
        Y = (D[ok].reshape(ok.sum(), -1)[:, fm] / sig).T
        C = photlib.Cutout(s); prf = photlib.PRFModel(C); prf.blur = json.load(open(f"{PHOT}/out/wcs_s{s}.json"))["blur"]
        series = {}
        base = [b.ravel()[fm], T3.ravel()[fm]]
        def fit(cols):
            A = np.column_stack(cols + [np.ones(fm.sum()), u, v]) / sig[:, None]
            return (np.linalg.pinv(A.T @ A) @ A.T @ Y)
        c = fit(base); series["M3V"] = c[1]
        for k, p in enumerate(phantoms(s)):
            Tp = prf.image(p[0] + x0, p[1] + y0)[y0:y0 + ny, x0:x0 + nx]
            c = fit(base + [Tp.ravel()[fm]]); series[f"phantom_{k}"] = c[2]
        # regressors
        J = np.load(f"{PHOT}/data/jitter_s{s}.npz")
        bk = z["bkg"][ok]; bexc = bk - np.median(bk)
        X = np.column_stack([np.ones(ok.sum()), bexc, z["dx"][ok] - np.nanmean(z["dx"][ok]), z["dy"][ok] - np.nanmean(z["dy"][ok])])
        tt = t[ok]
        kb = np.floor((tt - tt[0]) * 24).astype(int); ub, inv = np.unique(kb, return_inverse=True); nbin = np.bincount(inv)
        tb = np.bincount(inv, tt) / nbin
        res = {}
        binned = {}
        for nm, F in series.items():
            cc, *_ = np.linalg.lstsq(X, F, rcond=None); r = F - X[:, 1:] @ cc[1:]
            ch, _, _ = harm_fit(tt, r, np.ones(len(r)), 0.648528); r = r - (harm_X(tt, 0.648528) @ ch - ch[0])
            binned[nm] = np.bincount(inv, r) / nbin
        top3, _, _ = ls_top(tb, binned["M3V"], 1.0, 14.0)
        P3 = top3[0][0]
        v3 = sin_vec(tb, binned["M3V"], P3)
        for nm, y in binned.items():
            top, _, _ = ls_top(tb, y, 1.0, 14.0)
            vv = sin_vec(tb, y, P3)
            res[nm] = dict(top=top[:2], amp_at_P3_e_s=float(np.hypot(*vv)),
                           phase_diff_rad=float((np.arctan2(vv[1], vv[0]) - np.arctan2(v3[1], v3[0]) + np.pi) % (2 * np.pi) - np.pi),
                           rms_1h_e_s=float(np.std(y)))
        ph_amps = np.array([v["amp_at_P3_e_s"] for k, v in res.items() if k.startswith("phantom")])
        out[s] = dict(P_M3V_top=P3, M3V_amp_e_s=res["M3V"]["amp_at_P3_e_s"], phantom_amps_e_s=ph_amps.tolist(),
                      frac_phantoms_ge_M3V=float(np.mean(ph_amps >= res["M3V"]["amp_at_P3_e_s"])), series=res)
        print(f"S{s}: M3V top {P3:.2f} d amp {res['M3V']['amp_at_P3_e_s']:.1f} e/s; phantoms at that P: " +
              ", ".join(f"{a:.1f}" for a in ph_amps) + f" (frac >= M3V {out[s]['frac_phantoms_ge_M3V']:.2f}); phantom top periods " +
              ", ".join(f"{v['top'][0][0]:.1f}" for k, v in res.items() if k.startswith("phantom")), flush=True)
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/phantom_rotation.json")


if __name__ == "__main__":
    main()
