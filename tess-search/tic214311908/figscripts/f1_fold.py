"""f1: fold.png - phase-folded transit of the cleaned INV and REF light curves (clean/) and of the raw-pixel 3x3-box
aperture (referee3_clean q02 method, recomputed per cadence by f0_rawap.py), with the adopted transit model,
residuals, and the phase-0.5 window (no secondary eclipse).
Each +-0.30-d window gets its own linear (INV, REF) or quadratic + pointing (aperture) baseline, fitted jointly with a
free-amplitude copy of the transit shape; only the baseline is removed. Fold on the adopted ephemeris
(T0 = 3840.9150 BTJD_TDB, P = 5.691506 d). Bins of 15 min, inverse-variance weighted (errors: max of formal and
scatter/sqrt(n))."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd

style()
H = 24.0
EDGES = np.arange(-4.0, 4.0001, 0.25) / H        # 15-min bins, in days


def load_inv():
    d = pd.read_csv(f'{B}/clean/lc_M3V_INV_clean.csv', usecols=['time_btjd_tdb', 'sector', 'cadence_s', 'quality_ok', 'flux_raw_norm',
                                                            'flux_err', 'beta_transit', 'flare_flag'])
    d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].reset_index(drop=True)
    return d


def windows_inv(d, phase):
    """windows of half-width 0.30 d at the given orbital phase; acceptance as referee3_bd lc3.load_inv"""
    t = d.time_btjd_tdb.values
    ep = np.round((t - T0_ADOPT) / P_ADOPT - phase).astype(int)
    tc = T0_ADOPT + (ep + phase) * P_ADOPT
    dt = t - tc
    out = []
    for s in np.unique(d.sector):
        for e in np.unique(ep[(d.sector.values == s) & (np.abs(dt) < 0.30)]):
            m = (d.sector.values == s) & (ep == e) & (np.abs(dt) < 0.30)
            if np.sum(m & (np.abs(dt) < 0.03)) >= 2 and np.sum(m & (dt < -0.07)) >= 10 and np.sum(m & (dt > 0.07)) >= 10:
                out.append(np.where(m)[0])
    return out, dt


def detrend(dt, y, err, cad, idx_list, order=1, shape_fn=None):
    """per-window baseline removal; returns concatenated dt, detrended y, weights, cadence"""
    X, Y, W, C = [], [], [], []
    for idx in idx_list:
        s = 1 - shape_fn(dt[idx], cad[idx])
        w = 1 / err[idx] ** 2
        base, a = fit_window_baseline(dt[idx], y[idx], w, -s, order=order)
        X.append(dt[idx]); Y.append(y[idx] - base + 1); W.append(w); C.append(cad[idx])
    return np.concatenate(X), np.concatenate(Y), np.concatenate(W), np.concatenate(C)


def shape_transit(dt, cad):
    return model_by_cadence(dt, cad)


def shape_occ(dt, cad):
    """occultation-shaped template (same duration, uniform disk): used only to profile the phase-0.5 baselines"""
    return model_by_cadence(dt, cad, u1=0.0, u2=0.0)


def main():
    res = {}
    # ---------------- INV
    d = load_inv()
    err = d.flux_err.values * d.beta_transit.values
    idx, dt = windows_inv(d, 0.0)
    res['INV'] = detrend(dt, d.flux_raw_norm.values, err, d.cadence_s.values.astype(float), idx, 1, shape_transit) + (len(idx),)
    idx5, dt5 = windows_inv(d, 0.5)
    sec_inv = detrend(dt5, d.flux_raw_norm.values, err, d.cadence_s.values.astype(float), idx5, 1, shape_occ) + (len(idx5),)
    # ---------------- REF
    r = pd.read_csv(f'{B}/clean/lc_M3V_REF_clean.csv')
    r = r[np.isfinite(r.flux_norm) & np.isfinite(r.flux_err)].reset_index(drop=True)
    tc = T0_ADOPT + r.epoch.values * P_ADOPT
    dtr = r.time_btjd_tdb.values - tc
    win = r.sector.values * 1000 + r.epoch.values
    idxr = [np.where(win == wv)[0] for wv in np.unique(win)]
    res['REF'] = detrend(dtr, r.flux_norm.values, r.flux_err.values * r.beta_transit.values, r.cadence_s.values.astype(float), idxr, 1,
                         shape_transit) + (len(idxr),)
    # ---------------- raw-pixel aperture (already detrended by f0 with quadratic + pointing baseline)
    z = np.load(f'{W}/figscripts/data/rawap_box3.npz')
    m0 = z['phase'] == 0.0
    tca = T0_ADOPT + z['epoch'] * P_ADOPT
    dta = z['t'] - tca
    nwin_ap = len(np.unique(z['sector'][m0] * 1000 + z['epoch'][m0]))
    res['AP'] = (dta[m0], z['f'][m0], 1 / z['e'][m0] ** 2, z['cad'][m0].astype(float), nwin_ap)
    m5 = z['phase'] == 0.5
    tca5 = T0_ADOPT + (z['epoch'] + 0.5) * P_ADOPT
    sec_ap = (z['t'][m5] - tca5[m5], z['f'][m5], 1 / z['e'][m5] ** 2, z['cad'][m5].astype(float),
              len(np.unique(z['sector'][m5] * 1000 + z['epoch'][m5])))

    # ---------------- figure
    fig = plt.figure(figsize=(10.4, 6.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.35, 1], height_ratios=[2.6, 1], hspace=0.08, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0]); axr = fig.add_subplot(gs[1, 0], sharex=ax); ax5 = fig.add_subplot(gs[:, 1])
    labels = {'INV': ('INV PRF photometry (cleaned)', C_INV, 'o'), 'REF': ('REF PRF photometry (cleaned)', C_REF, 's'),
              'AP': ('Raw pixels, 3x3 box (no scene model)', AQUA_DARK, '^')}
    offs = {'INV': 0.0, 'REF': -32.0, 'AP': -64.0}
    roffs = {'INV': 0.0, 'REF': -9.0, 'AP': -18.0}
    tm = np.linspace(-4, 4, 1601) / H
    summary = {}
    for key in ('INV', 'REF', 'AP'):
        x, y, w, cad, nw = res[key]
        lab, col, mk = labels[key]
        mod_pts = model_by_cadence(x, cad)
        # smooth model: cadence mix weighted like the data
        wc = {c: w[cad == c].sum() for c in np.unique(cad)}
        mcurve = sum(wc[c] * transit_model(tm, exp_s=c) for c in wc) / sum(wc.values())
        xb, yb, eb = bin_curve(x, (y - 1) * 1e3, w, EDGES)
        ax.errorbar(xb * H, yb + offs[key], eb, fmt=mk, ms=4.2, color=col, mec='white', mew=0.5, elinewidth=1.0, zorder=3)
        ax.plot(tm * H, (mcurve - 1) * 1e3 + offs[key], '-', color=INK, lw=1.4, zorder=4)
        xr, yr, er = bin_curve(x, (y - mod_pts) * 1e3, w, EDGES)
        axr.axhline(roffs[key], color=MUTED, lw=0.8, zorder=1)
        ax.axhline(offs[key], color=GRID, lw=0.6, zorder=0)
        axr.errorbar(xr * H, yr + roffs[key], er, fmt=mk, ms=3.6, color=col, mec='white', mew=0.4, elinewidth=0.9, zorder=3)
        ax.text(-3.92, offs[key] + 5.5, f'{lab}  ({nw} transits)', color=col if key != 'AP' else AQUA_DARK, fontsize=8.6,
                ha='left', va='bottom', fontweight='bold')
        chi2 = float(np.sum(((yr) / er) ** 2)); summary[key] = dict(n_windows=int(nw), n_points=int(len(x)), resid_chi2_bins=chi2, n_bins=int(len(xr)))
    d0 = (1 - transit_model(np.array([0.0]))[0]) * 1e3
    aR = aR_of_rho(RHO_ADOPT)
    ax.text(3.95, -99, f'Black: adopted transit model, depth {d0:.1f} ppt\n(Rp/R* {K_ADOPT:.3f}, b {B_ADOPT:.2f}, a/R* {aR:.1f}, '
            f'T14 1.70 h); 15-min bins', ha='right', va='bottom', fontsize=8.0, color=INK2)
    ax.set_ylabel('Relative flux (ppt of the M3V flux, offset)')
    ax.set_ylim(-100, 12); ax.set_yticks([-80, -60, -40, -20, 0])
    ax.set_xlim(-4, 4)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title('TIC 214311908: folded TESS transit, six sectors (S12-S104)', loc='left')
    axr.set_ylabel('Residual (ppt, offset)')
    axr.text(0.01, 0.97, 'Data minus model; same order and symbols as panel a', transform=axr.transAxes, fontsize=7.8, color=INK2, va='top')
    axr.set_xlabel('Time from mid-transit (h)')
    axr.set_ylim(-27, 8)
    # ---------------- phase 0.5 panel
    o5 = {'INV': 0.0, 'AP': -14.0}
    occ_curve = (1 - (1 - transit_model(tm, u1=0., u2=0., exp_s=200)) / (1 - transit_model(np.array([0.]), u1=0., u2=0.)[0]) * 3.0e-3)
    for key, data, lab, col, mk in (('INV', sec_inv, 'INV (cleaned)', C_INV, 'o'), ('AP', sec_ap, 'Raw pixels, 3x3 box', AQUA_DARK, '^')):
        x, y, w, cad, nw = data
        xb, yb, eb = bin_curve(x, (y - 1) * 1e3, w, EDGES)
        ax5.errorbar(xb * H, yb + o5[key], eb, fmt=mk, ms=4.0, color=col, mec='white', mew=0.5, elinewidth=1.0, zorder=3)
        ax5.axhline(o5[key], color=MUTED, lw=0.8)
        ax5.plot(tm * H, (occ_curve - 1) * 1e3 + o5[key], '--', color=INK2, lw=1.0)
        ax5.text(-3.9, o5[key] + 6.0, f'{lab} ({nw} windows)', color=col, fontsize=8.2, fontweight='bold')
        summary[f'phase05_{key}'] = dict(n_windows=int(nw), n_points=int(len(x)))
    ph = json.load(open(f'{B}/clean/out/phase_scan.json'))
    q12 = json.load(open(f'{B}/referee3_clean/out/q12_secondary_ap.json'))
    si = ph['inv']['secondary_phase05']; sr = ph['ref']['secondary_phase05']; sa = q12['L_box3_pt_ph0.5']
    ax5.text(0.03, 0.035, f"Fitted occultation depth at phase 0.5:\n"
             f"INV {si['depth_ppt']:+.1f} $\\pm$ {si['err_cal_ppt']:.1f} ppt, REF (not shown) {sr['depth_ppt']:+.1f} $\\pm$ {sr['err_cal_ppt']:.1f} ppt\n"
             f"raw 3x3 box {sa['depth_ppt']:+.1f} $\\pm$ {sa['err']:.1f} ppt\n"
             f"Dashed: a 3-ppt occultation (~3$\\sigma$ limit)", transform=ax5.transAxes, fontsize=7.8, color=INK2, va='bottom')
    ax5.set_xlim(-4, 4); ax5.set_ylim(-27, 9)
    ax5.set_xlabel('Time from phase 0.5 (h)')
    ax5.set_ylabel('Relative flux (ppt, offset)')
    ax5.set_title('Phase 0.5: no secondary eclipse', loc='left')
    panel_label(ax, 'a'); panel_label(axr, 'b'); panel_label(ax5, 'c')
    save(fig, 'fold.png')
    summary['model_depth_ppt'] = d0
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
