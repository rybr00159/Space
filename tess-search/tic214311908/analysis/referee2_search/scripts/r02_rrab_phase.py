"""r02: is the 0.6485-d signal in the M3V light curve the (inverted) light curve of the Gaia DR3 RRab
5925372615693467904 (200" S)?

1. Gaia DR3 vari_rrlyrae (queried by referee 2, data/gaia_rrab.csv): pf = 0.648528 +- 0.000003 d, epoch_g of maximum
   1707.312788 (BJD_TCB - 2455197.5) -> predicted maximum-light times in BTJD (TCB-TDB difference < 1e-5 d ignored),
   propagated phase uncertainty.
2. Simple aperture photometry of the RRab directly from the TESScut cubes (3x3 px at its WCS position, clipped to the
   cutout; local background = median of the 5x5 annulus pixels below the 30th percentile), detrended with a 1-d
   running median, folded at pf.
3. Folded waveform of the M3V-B leak (phot det, transits masked, 30 bins) and its correlation with the RRab waveform
   at zero lag (sign and r); lag of the best cross-correlation in phase.
Outputs: out/r02_rrab_phase.json, figs/r02_rrab_fold.png
"""
import sys, os, glob
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CUBES = "/home/user/work/tesscut5/214311906"
NB = 30


def fold_bins(t, y, Pp, Tmax):
    ph = ((t - Tmax) / Pp) % 1
    k = np.floor(ph * NB).astype(int)
    med = np.array([np.nanmedian(y[k == j]) if np.any(k == j) else np.nan for j in range(NB)])
    return med


def main():
    g = pd.read_csv(f"{ME}/data/gaia_rrab.csv").iloc[0]
    Pp = float(g.pf); Pe = float(g.pf_error)
    Tmax0 = float(g.epoch_g) + 2455197.5 - 2457000.0
    out = dict(gaia=dict(pf=Pp, pf_err=Pe, Tmax_btjd=Tmax0, ra=float(g.ra), dec=float(g.dec), G=float(g.phot_g_mean_mag)))
    d = load_m3v("B")
    oot = ~d.intr.values
    res = {}
    fig, axs = plt.subplots(2, 3, figsize=(14, 7)); axs = axs.ravel()
    for i, sec in enumerate(SECTORS):
        fn = glob.glob(f"{CUBES}/tess-s{sec:04d}-*_astrocut.fits")[0]
        with fits.open(fn, memmap=True) as h:
            tab = h[1].data; w = WCS(h[2].header)
            tt = np.array(tab["TIME"], float); fl = np.array(tab["FLUX"], float); q = np.array(tab["QUALITY"], int)
        x, y = w.all_world2pix([[g.ra, g.dec]], 0)[0]
        ny, nx = fl.shape[1:]
        ix, iy = int(round(x)), int(round(y))
        xs = range(max(ix - 1, 0), min(ix + 2, nx)); ys = range(max(iy - 1, 0), min(iy + 2, ny))
        ap = np.zeros((ny, nx), bool)
        for a in ys:
            for b in xs: ap[a, b] = True
        ann = np.zeros((ny, nx), bool)
        for a in range(max(iy - 3, 0), min(iy + 4, ny)):
            for b in range(max(ix - 3, 0), min(ix + 4, nx)):
                if not ap[a, b]: ann[a, b] = True
        good = (q == 0) & np.isfinite(tt) & np.all(np.isfinite(fl[:, ap]), axis=1)
        tt, fl = tt[good], fl[good]
        bk = []
        for f_ in fl:
            v = f_[ann]; v = v[np.isfinite(v)]
            bk.append(np.mean(v[v <= np.percentile(v, 30)]))
        bk = np.array(bk)
        lc = fl[:, ap].sum(1) - ap.sum() * bk
        # 1-d running median detrend
        trend = np.array([np.median(lc[np.abs(tt - tc) < 0.5]) for tc in tt])
        rr = lc - trend
        # TESScut TIME is FFI-barycentred (target offset <~ 31 s: irrelevant here)
        Tmax_pred_phase_err = abs((np.median(tt) - Tmax0) / Pp) * Pe / Pp   # in phase units
        wf_rr = fold_bins(tt, rr, Pp, Tmax0)
        m = d.s.values == sec
        mm = m & oot
        wf_m3v = fold_bins(d.t.values[mm], (d.det.values[mm] - 1) * FREF, Pp, Tmax0)
        ok = np.isfinite(wf_rr) & np.isfinite(wf_m3v)
        r0 = float(np.corrcoef(wf_rr[ok], wf_m3v[ok])[0, 1])
        # cross-correlation over lags
        lags = np.arange(-NB // 2, NB // 2)
        cc = [np.corrcoef(np.roll(wf_rr, L)[ok], wf_m3v[ok])[0, 1] for L in lags]
        Lbest = int(lags[np.argmin(cc)])
        # phase of RRab maximum in the TESS data
        phmax = (np.argmax(wf_rr) + 0.5) / NB
        res[sec] = dict(rr_xy=[float(x), float(y)], n_ap=int(ap.sum()), rr_amp_ptp_e_s=float(np.nanmax(wf_rr) - np.nanmin(wf_rr)),
                        m3v_leak_ptp_e_s=float(np.nanmax(wf_m3v) - np.nanmin(wf_m3v)), corr_zero_lag=r0,
                        most_negative_corr=float(np.min(cc)), lag_of_most_negative_phase=Lbest / NB,
                        rr_max_phase_vs_gaia=float(phmax if phmax < 0.5 else phmax - 1), gaia_phase_err=float(Tmax_pred_phase_err),
                        ratio_leak_over_rr=float(-np.polyfit(wf_rr[ok], wf_m3v[ok], 1)[0]))
        print(sec, res[sec], flush=True)
        ph = (np.arange(NB) + 0.5) / NB
        axs[i].plot(ph, wf_rr / np.nanmax(np.abs(wf_rr)), "o-", label="RRab aperture (norm)")
        axs[i].plot(ph, wf_m3v / np.nanmax(np.abs(wf_m3v)), "s-", label="M3V-B leak (norm)")
        axs[i].set_title(f"S{sec} r={r0:+.2f}"); axs[i].set_xlabel("phase from Gaia max")
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{ME}/figs/r02_rrab_fold.png", dpi=80)
    out["sectors"] = res
    jdump(out, f"{ME}/out/r02_rrab_phase.json")


if __name__ == "__main__":
    main()
