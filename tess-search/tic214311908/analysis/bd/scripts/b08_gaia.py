"""b08: Gaia DR3 constraints on the BD channels.
(1) RVS radial-velocity error. The target has rv_method_used = 2 (grvs_mag 12.55): the DR3 RV and its error come from the
    combined (summed) CCF of rv_nb_transits = 21 transits (12 of them deblended), and rv_amplitude_robust,
    rv_chisq_pvalue and rv_renormalised_gof are NOT provided for such stars (NaN in the archive; checked). The only
    handle is radial_velocity_error = 2.217 km/s. Model (generalises ../referee2_shape/scripts/r13_rv_smear.py): each
    transit's CCF is a Gaussian of sigma w (12 or 16 km/s: RVS R = 11500 plus template matching) with peak C0 (0.3 or
    0.6); an orbit with semi-amplitude K and eccentricity e, sampled at 21 random mean anomalies, shifts each CCF by
    v_i = K [cos(f_i + w) + e cos w]; the error of the summed CCF follows Zucker (2003), sigma^2 = -[N C''/C C^2/(1 - C^2)]^-1
    at the peak. Ratio r = sigma(K)/sigma(0), distribution over phases (1000 draws). Intrinsic error distribution p0 from
    Gaia DR3 M dwarfs with rv_method_used = 2 (b08a: parallax > 8, 2.55 < BP-RP < 2.95, 12.2 < grvs_mag < 12.9),
    restricted to rv_expected_sig_to_noise within +-15 % of the target's 14.70 (KDE in log error). Likelihood ratio
    LR(K) = < p0(err_obs/r)/r >_phases / p0(err_obs).
    Companion-hosted case: the companion supplies a fraction f_c of the RVS light and moves with K2; the summed CCF is
    (1 - f_c) G(0) + f_c <G(v_i)>.
(2) RUWE. A BD on the M3V moves the photocentre (BD dark) on an edge-on orbit of angular semi-major axis
    alpha = a M_bd/(M* + M_bd) x parallax; with P = 5.69 d every CCD observation sees an independent displacement, adding
    <delta^2> = alpha^2/2 x <cos^2 theta> (= alpha^2/4 averaged over the node orientation; 0.15-0.35 alpha^2 for the
    anisotropic scan of this source) per CCD; UWE^2 -> UWE^2 + <delta^2>/sigma_CCD^2 (sigma_CCD 0.19-0.25 mas,
    ../referee2_shape/out/sigma_ccd_check.json). LR from the RUWE distribution of comparison M dwarfs of similar G.
Outputs: out/b08_gaia.json, out/b08_rv_lr.npz"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
from scipy.stats import gaussian_kde
rng = np.random.default_rng(8)
ERR_OBS, SN_OBS, RUWE_OBS, GMAG = 2.217367, 14.695826, 1.076204, 14.281
cmp_ = pd.read_csv(f'{BD}/data/gaia/rv_comparison.csv')
sel = cmp_[(np.abs(cmp_.rv_expected_sig_to_noise/SN_OBS - 1) < 0.15) & (cmp_.source_id != 5925373272850784768)]
le = np.log(sel.radial_velocity_error.values); kde = gaussian_kde(le)
p0 = lambda err: kde(np.log(err))/err      # density in err
out = dict(n_comparison=int(len(sel)), err_pct={p: float(np.percentile(sel.radial_velocity_error, p)) for p in (1, 5, 16, 50, 84, 95)},
           target_err_percentile=float(np.mean(sel.radial_velocity_error < ERR_OBS)*100),
           deblended_effect=dict(median_err_deblended_frac_gt_0p4=float(sel[sel.rv_nb_deblended_transits/sel.rv_nb_transits > 0.4].radial_velocity_error.median()) if np.any(sel.rv_nb_deblended_transits/sel.rv_nb_transits > 0.4) else None,
                                 n_deblended_frac_gt_0p4=int(np.sum(sel.rv_nb_deblended_transits/sel.rv_nb_transits > 0.4)),
                                 median_err_all=float(sel.radial_velocity_error.median())))
print(out, flush=True)
v = np.linspace(-160, 160, 1601); dv = v[1] - v[0]
def zucker(C):
    i = int(np.argmax(C)); i = min(max(i, 2), len(v) - 3)
    c = C[i]; c2 = (C[i+1] - 2*C[i] + C[i-1])/dv**2
    return np.sqrt(-1/(c2/c*c**2/(1 - c**2)))
def kepler_v(K, e, w, M):
    E = M.copy()
    for _ in range(30): E = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
    f = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
    return K*(np.cos(f + w) + e*np.cos(w))
Ks = np.array([0., 1., 2., 3., 4., 5., 6., 8., 10., 12., 14., 16., 18., 20., 25., 30., 40.])
es = np.array([0., 0.3, 0.5, 0.7])
NDR = 250
LR = {}
for (w_, C0) in ((12., 0.3), (16., 0.3), (20., 0.3), (12., 0.6), (16., 0.6), (20., 0.6)):
    base = zucker(C0*np.exp(-0.5*(v/w_)**2))
    tab = np.zeros((len(es), len(Ks))); rmed = np.zeros_like(tab)
    for ie, e in enumerate(es):
        for ik, K in enumerate(Ks):
            rr = []
            for j in range(NDR if K > 0 else 1):
                vi = kepler_v(K, e, rng.uniform(0, 2*np.pi), rng.uniform(0, 2*np.pi, 21)) if K > 0 else np.zeros(21)
                C = C0*np.mean(np.exp(-0.5*((v[:, None] - vi[None, :])/w_)**2), axis=1)
                rr.append(zucker(C)/base)
            rr = np.array(rr)
            tab[ie, ik] = float(np.mean(p0(ERR_OBS/rr)/rr)/p0(ERR_OBS)[0])
            rmed[ie, ik] = float(np.median(rr))
    LR[f'w{int(w_)}_C{C0}'] = dict(LR=tab, ratio_median=rmed)
    print(w_, C0, 'LR e=0:', np.round(tab[0], 3), flush=True)
# companion-hosted: companion RVS light fraction f_c, K2
fcs = np.array([0.05, 0.1, 0.2, 0.3]); K2s = np.array([5., 10., 20., 30., 40.])
LRc = np.zeros((len(fcs), len(K2s)))
w_, C0 = 14., 0.45
base = zucker(C0*np.exp(-0.5*(v/w_)**2))
for i, fc in enumerate(fcs):
    for j, K2 in enumerate(K2s):
        rr = []
        for _ in range(200):
            vi = K2*np.sin(rng.uniform(0, 2*np.pi, 21))
            C = C0*((1 - fc)*np.exp(-0.5*(v/w_)**2) + fc*np.mean(np.exp(-0.5*((v[:, None] - vi[None, :])/w_)**2), axis=1))
            rr.append(zucker(C)/base)
        rr = np.array(rr)
        LRc[i, j] = float(np.mean(p0(ERR_OBS/rr)/rr)/p0(ERR_OBS)[0])
print('companion LR', np.round(LRc, 3), flush=True)
np.savez(f'{BD}/out/b08_rv_lr.npz', Ks=Ks, es=es, **{f'LR_{k}': v_['LR'] for k, v_ in LR.items()}, **{f'rmed_{k}': v_['ratio_median'] for k, v_ in LR.items()},
         fcs=fcs, K2s=K2s, LR_comp=LRc)
out['rv_LR'] = {k: dict(Ks=Ks, es=es, LR=v_['LR'], ratio_median=v_['ratio_median']) for k, v_ in LR.items()}
out['rv_LR_companion'] = dict(fcs=fcs, K2s=K2s, LR=LRc, w=w_, C0=C0)
# ---------------- RUWE
m = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
rs = m[(np.abs(m.phot_g_mean_mag - GMAG) < 0.6) & np.isfinite(m.ruwe)].ruwe.values
rs = rs[rs < 3]
kr = gaussian_kde(np.log(rs))
pr = lambda x: kr(np.log(x))/x
plx = 15.497
Mbds = np.array([13, 20, 30, 40, 50, 60, 70, 78.5])
ru = {}
for sccd in (0.19, 0.22, 0.25):
    row = []
    for Mb in Mbds:
        q = Mb/1047.57
        a1 = 0.0448*q/(0.370 + q)*plx                 # mas
        lrs = []
        for cth in (0.15, 0.25, 0.35):
            dU2 = a1**2/2*cth*2/sccd**2               # <cos^2 theta> = cth*2 -> 0.3..0.7 of alpha^2/2
            x0 = RUWE_OBS**2 - dU2
            lrs.append(float(pr(np.sqrt(x0))[0]*RUWE_OBS/np.sqrt(x0)/pr(RUWE_OBS)[0]) if x0 > 0.3 else 0.0)
        row.append(dict(M_bd=float(Mb), alpha_mas=float(a1), LR_by_projection=lrs, LR_mean=float(np.mean(lrs))))
    ru[f'sigma_ccd_{sccd}'] = row
    print('RUWE', sccd, [(r['M_bd'], round(r['alpha_mas'], 3), round(r['LR_mean'], 3)) for r in row], flush=True)
out['ruwe'] = dict(n_comparison=int(len(rs)), ruwe_pct={p: float(np.percentile(rs, p)) for p in (16, 50, 84, 95)}, target=RUWE_OBS, table=ru)
jdump(dict(method=__doc__, **out), f'{BD}/out/b08_gaia.json')
