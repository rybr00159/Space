"""Referee r09: the M3V's TESS flux in the same (PRF-fit) units as the transit flux loss, from Gaia XP spectra.
Synthetic TESS photon rate N = int F_lambda * lambda * R(lambda) dlambda (SVO TESS/TESS.Red response; XP sampled
spectra 336-1020 nm; the part of the band above 1020 nm is added by extending F_lambda flat from its 1000-1020 nm mean,
and as a variant ignored). For the r01b field stars (isolated, T < 14.5) the ratio PRF-fit flux / N is formed per sector;
its median times N(M3V) is the M3V's flux in PRF units. Also the implied TESS magnitude using the TIC-Tmag zero point of
the same stars. Output: out/xpflux.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
X = pd.read_csv(f"{R}/data/ext/xp_sampled.csv")
rr = np.loadtxt(f"{R}/data/ext/tess_red.dat"); wlr = rr[:, 0] / 10; Rr = rr[:, 1]
FS = pd.read_csv(f"{R}/out/fluxscale_stars.csv")
g = pd.read_csv(f"{R}/data/gaia330.csv")


def synth(sid, ext=True):
    d = X[X.source_id == sid]
    if not len(d): return np.nan
    wl = d.wl.values; fl = d.flux.values
    R_ = np.interp(wl, wlr, Rr, left=0, right=0)
    N = np.trapezoid(fl * wl * R_, wl)
    if ext:
        tail = np.arange(1022, 1130, 2.0); Rt = np.interp(tail, wlr, Rr, left=0, right=0)
        fl_t = np.mean(fl[(wl >= 1000)])
        N += np.trapezoid(fl_t * tail * Rt, tail)
    return N


out = {}
for ext in (True, False):
    Nm = synth(GID["M3V"], ext); Nb = synth(GID["bg3"], ext)
    res = {}
    allr = []
    for s in SECTORS:
        d = FS[(FS.sector == s) & FS.iso & FS.inside & (FS.Ttic < 14.5)].copy()
        d["N"] = [synth(i, ext) for i in d.source_id]
        d = d[np.isfinite(d.N)]
        d["k"] = d.fit / d.N
        d["zp"] = d.Ttic + 2.5 * np.log10(d.N)        # T = zp - 2.5 log N
        kmed = float(np.median(d.k)); kerr = float(1.2533 * rstd(d.k) / np.sqrt(len(d)))
        res[s] = dict(n=len(d), k_median=kmed, k_rstd=float(rstd(d.k)), k_err=kerr, F_M3V_prf=kmed * Nm, F_M3V_prf_err=kerr * Nm,
                      T_M3V_xp=float(np.median(d.zp) - 2.5 * np.log10(Nm)), T_zp_rstd=float(rstd(d.zp)),
                      colour_slope=float(np.polyfit(d.bp_rp, d.k / kmed, 1)[0]) if len(d) > 4 else np.nan, F_bg3_prf=kmed * Nb)
        allr.append(d.assign(sector=s))
    A = pd.concat(allr)
    out["ext" if ext else "noext"] = dict(per_sector=res, F_M3V_prf_mean=float(np.mean([v["F_M3V_prf"] for v in res.values()])),
                                          T_M3V_xp_median=float(np.median([v["T_M3V_xp"] for v in res.values()])),
                                          F_bg3_prf_mean=float(np.mean([v["F_bg3_prf"] for v in res.values()])),
                                          T_zp_scatter_all=float(rstd(A.zp)), n_stars=int(A.source_id.nunique()),
                                          bp_rp_range=[float(A.bp_rp.min()), float(A.bp_rp.max())])
    print(f"[{'with' if ext else 'without'} >1020 nm tail] " + "; ".join(f"S{s}: k {v['k_median']:.4g} (rstd {v['k_rstd']/v['k_median']*100:.1f}%, n {v['n']}), F_M3V {v['F_M3V_prf']:.0f}+-{v['F_M3V_prf_err']:.0f}, T_xp {v['T_M3V_xp']:.3f}, colour slope {v['colour_slope']:+.3f}/mag"
                                                                    for s, v in res.items()))
    print(f"   mean F_M3V(PRF units) {out['ext' if ext else 'noext']['F_M3V_prf_mean']:.0f} e/s; T_M3V(XP, TIC zero point) {out['ext' if ext else 'noext']['T_M3V_xp_median']:.3f}; "
          f"F(bg3) {out['ext' if ext else 'noext']['F_bg3_prf_mean']:.1f} e/s; zp scatter {rstd(A.zp):.3f} mag over {A.source_id.nunique()} stars BP-RP {A.bp_rp.min():.2f}-{A.bp_rp.max():.2f}")
save_json(out, f"{R}/out/xpflux.json")
