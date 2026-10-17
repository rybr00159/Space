"""f5: scenarios.png - the folded cleaned REF light curve (clean/lc_M3V_REF_clean.csv, 21 windows) with four
best-fit models, recomputed from saved parameters:
  planet        referee3_bd/out/s03_fit_real.json  ref_beta.planet   (cleaned REF, clean/ beta noise)
  grazing BD    referee3_bd/out/s03_fit_real.json  ref_beta.grazBD   (k >= 0.19, undiluted, on the M3V)
  V-shaped      referee2_shape/out/scen_ref_fake.json vshape         (phase-1 REF curve; b >= 1 - k)
  WD occulted   referee2_shape/out/scen_ref_fake.json wd_cool        (phase-1 REF curve; WD totally occulted by the M3V)
The referee2_shape models were fitted to the phase-1 REF photometry, not to the cleaned curve; they are placed at
their mid-time offset relative to the planet fit of the same run. Per-window linear baselines are those of the
planet fit. Delta-chi2 values quoted in the legend come from the files named above (not refitted here); a check of
each model on the cleaned REF data with its own profiled baselines is printed to stdout."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import pandas as pd

style()
H = 24.0
RSUN = 6.957e10


def overlap_frac(d, k):
    d = np.asarray(d, float); out = np.zeros_like(d)
    out[d <= k - 1] = 1.0
    inside = d <= 1 - k
    out[inside] = k ** 2
    part = (d > k - 1) & (~inside) & (d < 1 + k)
    dd = d[part]
    a1 = np.arccos(np.clip((dd ** 2 + 1 - k ** 2) / (2 * dd), -1, 1))
    a2 = np.arccos(np.clip((dd ** 2 + k ** 2 - 1) / (2 * dd * k), -1, 1))
    area = a1 + k ** 2 * a2 - 0.5 * np.sqrt(np.clip((-dd + 1 + k) * (dd + 1 - k) * (dd - 1 + k) * (dd + 1 + k), 0, None))
    out[part] = area / np.pi
    return out


def wd_model(dt, cad, w, P):
    """uniform WD of radius RA (fraction fA of the light) occulted by the M3V (radius RB), circular orbit"""
    a, RA, RB, fA = w['a_rsun'], w['R_WD'], w['R_occ'], w['fA']
    cosi = w['b_occ'] * RB / a
    out = np.ones(len(dt))
    for c in np.unique(cad):
        m = cad == c
        ns = 5 if c < 300 else (11 if c < 900 else 31)
        off = ((np.arange(ns) + 0.5) / ns - 0.5) * c / 86400.
        tt = (dt[m][:, None] + off[None, :]).ravel()
        ph = 2 * np.pi * tt / P
        x = a * np.sin(ph); y = a * cosi * np.cos(ph)
        dd = np.sqrt(x ** 2 + y ** 2) / RA
        fr = np.where(np.cos(ph) > 0, overlap_frac(dd, RB / RA), 0.)
        out[m] = (1 - fA * fr).reshape(-1, ns).mean(1)
    return out


def main():
    s3 = json.load(open(f'{B}/referee3_bd/out/s03_fit_real.json'))['fits']['ref_beta']
    sc = json.load(open(f'{B}/referee2_shape/out/scen_ref_fake.json'))
    r = pd.read_csv(f'{B}/clean/lc_M3V_REF_clean.csv')
    r = r[np.isfinite(r.flux_norm) & np.isfinite(r.flux_err)].reset_index(drop=True)
    dt = r.time_btjd_tdb.values - (T0_CLEAN + r.epoch.values * P_CLEAN)
    cad = r.cadence_s.values.astype(float)
    y = r.flux_norm.values; e = r.flux_err.values * r.beta_transit.values; w = 1 / e ** 2
    win = r.sector.values * 1000 + r.epoch.values
    idx = [np.where(win == v)[0] for v in np.unique(win)]
    pl, bd = s3['planet'], s3['grazBD']
    tp = sc['planet_free']['x'][0]
    vx = sc['vshape']['x']; wdp = sc['wd_cool']
    models = {}
    models['planet'] = lambda t, cc: model_by_cadence(t - pl['dT0_min'] / 1440, cc, k=pl['k'], b=pl['b'], rho=pl['rho_c'], u1=pl['u1'], u2=pl['u2'], per=P_CLEAN)
    models['bd'] = lambda t, cc: model_by_cadence(t - bd['dT0_min'] / 1440, cc, k=bd['k'], b=bd['b'], rho=bd['rho_c'], u1=bd['u1'], u2=bd['u2'], per=P_CLEAN)
    vb = sc['vshape']['b']; vrho = 10 ** vx[4]
    dTv = pl['dT0_min'] / 1440 + (vx[0] - tp)
    models['v'] = lambda t, cc: model_by_cadence(t - dTv, cc, k=vx[2], b=vb, rho=vrho, u1=vx[5], u2=vx[6], per=P_CLEAN)
    dTw = pl['dT0_min'] / 1440 + (wdp['x'][0] - tp)
    models['wd'] = lambda t, cc: wd_model(t - dTw, cc, wdp, P_CLEAN)
    # evaluate at data points and check chi2 with own baselines
    mp = {k: f(dt, cad) for k, f in models.items()}
    chi = {}
    for k, m in mp.items():
        c2 = 0.0
        for ii in idx:
            X = np.column_stack([np.ones(len(ii)), dt[ii]])
            sw = np.sqrt(w[ii])
            c, *_ = np.linalg.lstsq(X * sw[:, None], (y[ii] - m[ii]) * sw, rcond=None)
            c2 += np.sum(w[ii] * (y[ii] - m[ii] - X @ c) ** 2)
        chi[k] = c2
    print('chi2 on cleaned REF with own linear baselines (fixed shapes):', {k: round(v - chi['planet'], 2) for k, v in chi.items()})
    print('file values: grazBD Lambda', round(s3['Lambda_grazBD'], 2), '| referee2_shape (phase-1 REF): vshape', round(sc['vshape']['dchi2_vs_planet_free'], 2),
          'wd_cool', round(sc['wd_cool']['dchi2_vs_planet_free'], 2))
    # detrend with the planet baselines
    yd = np.empty_like(y)
    for ii in idx:
        X = np.column_stack([np.ones(len(ii)), dt[ii]]); sw = np.sqrt(w[ii])
        c, *_ = np.linalg.lstsq(X * sw[:, None], (y[ii] - mp['planet'][ii]) * sw, rcond=None)
        yd[ii] = y[ii] - (X @ c)          # planet-fit baseline removed (flux_norm is 1 + dF/979)
    EDG = np.arange(-3.0, 3.0001, 0.2) / H          # 12-min bins
    xb, ybn, eb = bin_curve(dt, (yd - 1) * 1e3, w, EDG)
    # smooth model curves (cadence mix of the data)
    tm = np.linspace(-3, 3, 1201) / H
    wc = {c: w[cad == c].sum() for c in np.unique(cad)}
    def smooth(key):
        tot = 0
        for c, wt in wc.items():
            tot = tot + wt * models[key](tm, np.full(len(tm), c))
        return tot / sum(wc.values())
    curves = {k: smooth(k) for k in models}
    fig = plt.figure(figsize=(10.2, 6.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.3, 1], hspace=0.07)
    ax = fig.add_subplot(gs[0]); axr = fig.add_subplot(gs[1], sharex=ax)
    ax.errorbar(xb * H, ybn, eb, fmt='o', ms=4.2, color=MUTED, mec='white', mew=0.4, elinewidth=1.0, zorder=2,
                label=f'REF cleaned photometry, 21 transits (12-min bins)')
    sty = {'planet': (INK, '-', 2.0, f"Planet (Rp/R* {pl['k']:.3f}, b {pl['b']:.2f}, flat-bottomed): reference, $\\Delta\\chi^2$ = 0"),
           'wd': (ORANGE, '--', 1.7, f"White dwarf totally occulted by the M3V ({wdp['M_WD']:.2f} Msun, ingress {wdp['ingress_min']:.0f} min): "
                                      f"$\\Delta\\chi^2$ +{sc['wd_cool']['dchi2_vs_planet_free']:.1f} light curve only"),
           'bd': (BLUE, '-.', 1.7, f"Brown dwarf grazing the M3V (k {bd['k']:.2f}, b {bd['b']:.2f}, needs e $\\geq$ 0.5): $\\Delta\\chi^2$ +{s3['Lambda_grazBD']:.1f}"),
           'v': (AQUA_DARK, ':', 2.0, f"V-shaped (grazing, b {vb:.2f}, k {vx[2]:.2f}): $\\Delta\\chi^2$ +{sc['vshape']['dchi2_vs_planet_free']:.1f}")}
    for k in ('wd', 'bd', 'v', 'planet'):
        col, ls, lw, lab = sty[k]
        ax.plot(tm * H, (curves[k] - 1) * 1e3, ls=ls, color=col, lw=lw, zorder=4 if k == 'planet' else 3, label=lab)
        axr.plot(tm * H, (curves[k] - curves['planet']) * 1e3, ls=ls, color=col, lw=lw, zorder=3)
    xr, yr, er = bin_curve(dt, (yd - mp['planet']) * 1e3, w, EDG)
    ax.set_ylabel('Relative flux (ppt)')
    ax.set_ylim(-33, 4)
    plt.setp(ax.get_xticklabels(), visible=False)
    h, l = ax.get_legend_handles_labels()
    order = [4, 3, 0, 1, 2] if len(h) == 5 else list(range(len(h)))
    ax.legend([h[i] for i in order], [l[i] for i in order], loc='lower center', fontsize=7.7, handlelength=2.6, borderaxespad=0.2, ncol=1)
    ax.set_title('$\\bf{a}$   False-positive shapes against the cleaned REF light curve', loc='left')
    axr.errorbar(xr * H, yr, er, fmt='o', ms=3.8, color=MUTED, mec='white', mew=0.4, elinewidth=0.9, zorder=2)
    axr.axhline(0, color=INK, lw=1.4, zorder=1)
    axr.set_ylim(-5.5, 5.5)
    axr.set_ylabel('Data or model\nminus planet (ppt)')
    axr.set_xlabel('Time from mid-transit (h)')
    axr.set_xlim(-3, 3)
    axr.text(0.01, 0.95, '$\\bf{b}$   Residuals from the planet model, with each alternative model minus the planet',
             transform=axr.transAxes, fontsize=8.6, va='top', color=INK)
    save(fig, 'scenarios.png')
    print('binned residual chi2 vs planet', round(float(np.sum((yr / er) ** 2)), 1), 'n bins', len(yr))
    for k in curves:
        print(k, 'depth at centre (ppt)', round((1 - curves[k][600]) * 1e3, 2))


if __name__ == '__main__':
    main()
