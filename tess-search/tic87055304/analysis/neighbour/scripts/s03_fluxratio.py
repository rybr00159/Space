"""TESS-band flux ratio r = F_nb / F_target and the dilution factors for each light curve.

A. r from catalogue photometry:
   A1 Stassun+19 (TIC-8) colour relation applied to Gaia DR3 G, BP, RP of both stars (Monte Carlo over photometric
      errors; the relation's 0.006-mag scatter is added in quadrature per star as an upper bound);
   A2 TIC v8.2 Tmag values (DR2-based);
   A3 Gaia DR3 XP sampled spectra (336-1020 nm) integrated (photon counting) over the TESS response; validated by
      integrating over G and comparing with the G-band flux ratio.
   (A4, SED synthetic TESS magnitudes, is in s04_stellar.py.)
B. Fraction of each star's PRF inside the QLP apertures (TESS_PRF at the CCD position, Gaia DR3 positions propagated
   to each sector) -> dilution D for the QLP SAP target light curve under two readings of QLP's flux normalisation:
     raw aperture flux:          D_raw  = (F_ap - F_n f_n) / (F_n f_n)
     target-dilution-corrected:  D_corr = (F_t f_t) / (F_n f_n) - 1   (first-round MCMC used 4.35)
C. Empirical calibration: the same template-amplitude estimator on QLP SAP and on the PSF pair light curve, same
   epochs.  rho_s = depth_QLP / depth_pair; the neighbour's undiluted depth is depth_pair (1 + 1/r), so the effective
   QLP dilution is 1 + D_eff = (1 + 1/r) / rho_s.
Output: results/fluxratio.json, fig/fluxratio_dilution.png
"""
import os, sys, json, time, glob, warnings
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbtools as nt

ROOT = nt.ROOT
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
CAT = os.path.join(ROOT, "data/cat")
NB, TG = 6720986208895208960, 6720986204578421888
g = pd.read_csv("/home/user/work/tic87/data/gaia_dr3_90arcsec.csv").set_index("source_id")
res = {}
rng = np.random.default_rng(1)


def stassun(G, c):
    return G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473


# ---------------------------------------------------------------- A1
def mags(sid, n):
    r = g.loc[sid]
    out = {}
    for b in ("g", "bp", "rp"):
        fl, fe = r[f"phot_{b}_mean_flux"], r[f"phot_{b}_mean_flux_error"]
        out[b] = r[f"phot_{b}_mean_mag"] - 2.5 * np.log10(1 + fe / fl * rng.standard_normal(n))
    return out


n = 20000
mn, mt = mags(NB, n), mags(TG, n)
Tn = stassun(mn["g"], mn["bp"] - mn["rp"]); Tt = stassun(mt["g"], mt["bp"] - mt["rp"])
r_A1 = 10 ** (-0.4 * (Tn - Tt))
# relation scatter: 0.006 mag per star (Stassun+19, for this colour range); for two stars of almost identical colour
# (BP-RP 0.928 vs 0.861) most of it is common-mode, so treat 0.006 mag on the DIFFERENCE as a conservative bound
r_A1_err = np.hypot(np.std(r_A1), np.median(r_A1) * 0.4 * np.log(10) * 0.006)
res["A1_stassun_gaiadr3"] = dict(T_nb=float(np.median(Tn)), T_tg=float(np.median(Tt)), r=float(np.median(r_A1)),
                                 r_err=float(r_A1_err), r_err_phot_only=float(np.std(r_A1)),
                                 G_ratio=float(g.loc[NB, "phot_g_mean_flux"] / g.loc[TG, "phot_g_mean_flux"]),
                                 RP_ratio=float(g.loc[NB, "phot_rp_mean_flux"] / g.loc[TG, "phot_rp_mean_flux"]))
# ---------------------------------------------------------------- A2
try:
    tn = pd.read_csv(os.path.join(CAT, "tic82_87055315.csv")).iloc[0]
    tt = pd.read_csv(os.path.join(CAT, "tic82_87055304.csv")).iloc[0]
    rr = 10 ** (-0.4 * (tn.Tmag - tt.Tmag))
    res["A2_tic82"] = dict(T_nb=float(tn.Tmag), eT_nb=float(tn.e_Tmag), T_tg=float(tt.Tmag), eT_tg=float(tt.e_Tmag),
                           r=float(rr), r_err=float(rr * 0.4 * np.log(10) * np.hypot(tn.e_Tmag, tt.e_Tmag)))
except Exception as e:
    res["A2_tic82"] = f"failed: {e}"
# ---------------------------------------------------------------- A3 XP sampled spectra
import sedlib  # local copy; filter curves from stellar/cache (read-only)
fxp = os.path.join(ROOT, "data/xp_sampled.npz")
if not os.path.exists(fxp):
    from astroquery.gaia import Gaia
    got = None
    for i in range(4):
        try:
            got = Gaia.load_data(ids=[NB, TG], data_release="Gaia DR3", retrieval_type="XP_SAMPLED",
                                 data_structure="INDIVIDUAL", format="csv", verbose=False)
            break
        except Exception as e:
            print("XP retry", e); time.sleep(20)
    if got is not None:
        spec = {}
        for key, tabs in got.items():
            for t in (tabs if isinstance(tabs, list) else [tabs]):
                df = t.to_pandas()
                sid = int(df["source_id"].iloc[0]) if "source_id" in df else int(str(key).split("_")[-1].split(".")[0])
                def arr(x):
                    if isinstance(x, str): return np.array([float(v) for v in x.strip("()[]").split(",")])
                    return np.array(x, float)
                if "flux" in df and len(df) == 1:
                    spec[sid] = (arr(df["flux"].iloc[0]), arr(df["flux_error"].iloc[0]))
                elif "flux" in df:
                    spec[sid] = (df["flux"].values.astype(float), df["flux_error"].values.astype(float))
        if NB in spec and TG in spec:
            np.savez(fxp, wl=np.arange(336, 1021, 2.0), fnb=spec[NB][0], enb=spec[NB][1], ftg=spec[TG][0], etg=spec[TG][1])
if os.path.exists(fxp):
    z = np.load(fxp)
    wl = z["wl"] * 10.0  # Angstrom
    def band(fl, b):
        f = sedlib.FILT[b]
        T = np.interp(wl, f["w"], f["t"], left=0, right=0)
        return np.trapezoid(fl * T * wl, wl)
    out = {}
    for b in ("TESS", "G", "RP"):
        # MC over XP flux errors (uncorrelated approximation; the sampled-spectrum errors are correlated, so this is
        # indicative only)
        rs = []
        for k in range(500):
            a = z["fnb"] + z["enb"] * rng.standard_normal(len(wl)); c = z["ftg"] + z["etg"] * rng.standard_normal(len(wl))
            rs.append(band(a, b) / band(c, b))
        out[b] = dict(r=float(band(z["fnb"], b) / band(z["ftg"], b)), r_err_mc=float(np.std(rs)))
        tf = sedlib.FILT[b]
        out[b]["frac_response_beyond_1020nm"] = float(np.trapezoid(tf["t"][tf["w"] > 10200] * tf["w"][tf["w"] > 10200], tf["w"][tf["w"] > 10200]) /
                                                      np.trapezoid(tf["t"] * tf["w"], tf["w"]))
    out["G_check"] = dict(xp_G_ratio=out["G"]["r"], gaia_G_ratio=res["A1_stassun_gaiadr3"]["G_ratio"])
    res["A3_xp_sampled"] = out
else:
    res["A3_xp_sampled"] = "XP sampled spectra not retrieved"
print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float), flush=True)

# ---------------------------------------------------------------- B: QLP aperture fractions
import psfcommon as pc
raw = np.load("/home/user/work/tic87/data/qlp_raw_columns.npy", allow_pickle=True).item()
ap_info = {}


def circ_weights(xc, yc, rad, shape=(21, 21), ss=20):
    ny, nx = shape
    w = np.zeros(shape)
    o = (np.arange(ss) + 0.5) / ss - 0.5
    ox, oy = np.meshgrid(o, o)
    for j in range(ny):
        for i in range(nx):
            if abs(i - xc) > rad + 1 or abs(j - yc) > rad + 1: continue
            w[j, i] = np.mean(np.hypot(i + ox - xc, j + oy - yc) <= rad)
    return w


fig_rows = []
for s in nt.SECTORS:
    sec = pc.load_sector(s)
    st = pc.gaia_stars(sec, gmax=19.0)
    prf = pc.PRFModel(sec)
    it = int(np.where(st.Source.values == TG)[0][0]); inb = int(np.where(st.Source.values == NB)[0][0])
    xt, yt = st.x.values[it], st.y.values[it]
    # QLP aperture centre: SAP_X/SAP_Y are CCD coordinates; cutout pixel (0,0) is CCD (col0, row0) per psfcommon
    qx = np.nanmedian(np.asarray(raw[s]["SAP_X"], float)); qy = np.nanmedian(np.asarray(raw[s]["SAP_Y"], float))
    cand = {"qlp_minus_col0": (qx - sec["col0"], qy - sec["row0"]), "target_gaia": (xt, yt)}
    flux = 10 ** (-0.4 * (st["T"].values - 10.0))
    imgs = {}
    def frac_all(wts):
        return np.array([np.sum(prf.img(x, y) * wts) for x, y in zip(st.x.values, st.y.values)])
    out = {}
    for cname, (cx, cy) in cand.items():
        if s == 104:
            # QLP '3x3 square': the 3x3 block of pixels centred on the pixel containing the aperture centre
            w = np.zeros((21, 21)); ic, jc = int(np.round(cx)), int(np.round(cy)); w[jc - 1:jc + 2, ic - 1:ic + 2] = 1
            wts = {"3x3": w}
        else:
            wts = {"circ2.5_area": circ_weights(cx, cy, 2.5), "circ2.5_centres": (circ_weights(cx, cy, 2.5, ss=1) > 0).astype(float)}
        for wn, w in wts.items():
            fr = frac_all(w)
            Fap = np.sum(flux * fr)
            Ftap, Fnap = flux[it] * fr[it], flux[inb] * fr[inb]
            out[f"{cname}|{wn}"] = dict(centre=[float(cx), float(cy)], frac_t=float(fr[it]), frac_n=float(fr[inb]),
                                        F_ap_over_Ft=float(Fap / flux[it]), others_over_Ft=float((Fap - Ftap - Fnap) / flux[it]),
                                        D_raw=float((Fap - Fnap) / Fnap), D_corr=float(Ftap / Fnap - 1),
                                        D_corr_gaia_r=float((fr[it] / fr[inb]) / res["A1_stassun_gaiadr3"]["r"] - 1),
                                        D_target_host_raw=float((Fap - Ftap) / Ftap))
    ap_info[s] = dict(target_xy=[float(xt), float(yt)], nb_xy=[float(st.x.values[inb]), float(st.y.values[inb])],
                      sep_px=float(np.hypot(st.x.values[inb] - xt, st.y.values[inb] - yt)), apertures=out)
    print(s, json.dumps(ap_info[s], default=float)[:900], flush=True)
res["B_qlp_apertures"] = ap_info

# ---------------------------------------------------------------- C: empirical QLP / pair depth ratio
pair = nt.load_lc("pair"); qlp = nt.load_lc("qlp")
emp = {}
r_ad = res["A1_stassun_gaiadr3"]["r"]
for s in nt.SECTORS:
    a, q = pair[s], qlp[s]
    # common epochs
    ra = nt.depth_fit(a["t"], a["f_raw"], a["exptime"]); rq = nt.depth_fit(q["t"], q["f_raw"], q["exptime"])
    com = sorted(set(ra["epochs"]) & set(rq["epochs"]))
    ra = nt.depth_with_errors(a["t"], a["f_raw"], a["exptime"], epochs=com)
    rq = nt.depth_with_errors(q["t"], q["f_raw"], q["exptime"], epochs=com)
    # per-epoch ratio scatter (the two series share the photon noise of the same pixels, so the ratio error is
    # estimated from the per-epoch differences, not by adding the errors)
    pa, pq = np.array(ra["per_epoch"]), np.array(rq["per_epoch"])
    rho = rq["A"] / ra["A"]
    # error on rho from the scatter of per-epoch differences q - rho*a
    dd = pq - rho * pa
    rho_err = np.std(dd, ddof=1) / np.sqrt(len(dd)) / ra["A"] if len(dd) > 2 else np.nan
    emp[s] = dict(n_common=len(com), depth_pair=ra["A"], depth_pair_err=ra["err"], depth_qlp=rq["A"], depth_qlp_err=rq["err"],
                  rho=float(rho), rho_err=float(rho_err),
                  rho_pred_corr=float((1 + r_ad) * ap_info[s]["apertures"][[k for k in ap_info[s]["apertures"] if k.startswith("qlp")][0]]["frac_n"] /
                                      ap_info[s]["apertures"][[k for k in ap_info[s]["apertures"] if k.startswith("qlp")][0]]["frac_t"]),
                  D_eff=float((1 + 1 / r_ad) / rho - 1), D_eff_err=float((1 + 1 / r_ad) / rho ** 2 * rho_err))
    print(s, json.dumps(emp[s], default=float), flush=True)
rv = [emp[s]["rho"] for s in nt.SECTORS]; re_ = [emp[s]["rho_err"] for s in nt.SECTORS]
m4, e4, c4 = nt.wmean(rv[:4], re_[:4])
m5, e5, c5 = nt.wmean(rv, re_)
res["C_empirical"] = dict(per_sector=emp, rho_S13_S93=dict(mean=m4, err=e4, chi2=c4, dof=3),
                          rho_all=dict(mean=m5, err=e5, chi2=c5, dof=4),
                          D_eff_S13_S93=float((1 + 1 / r_ad) / m4 - 1), D_eff_S13_S93_err=float((1 + 1 / r_ad) / m4 ** 2 * e4),
                          D_eff_S104=emp[104]["D_eff"], D_eff_S104_err=emp[104]["D_eff_err"])
json.dump(res, open(os.path.join(ROOT, "results/fluxratio.json"), "w"), indent=1, default=float)

# ---------------------------------------------------------------- figure
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
labs = ["Stassun+19\n(Gaia DR3)", "TIC v8.2", "XP sampled\n(TESS band)"]
vals = [res["A1_stassun_gaiadr3"]["r"], res["A2_tic82"]["r"] if isinstance(res["A2_tic82"], dict) else np.nan,
        res["A3_xp_sampled"]["TESS"]["r"] if isinstance(res["A3_xp_sampled"], dict) else np.nan]
errs = [res["A1_stassun_gaiadr3"]["r_err"], res["A2_tic82"]["r_err"] if isinstance(res["A2_tic82"], dict) else np.nan,
        res["A3_xp_sampled"]["TESS"]["r_err_mc"] if isinstance(res["A3_xp_sampled"], dict) else np.nan]
ax[0].errorbar(range(3), vals, yerr=errs, fmt="o")
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labs, fontsize=8); ax[0].set_ylabel("F_nb / F_target (TESS band)")
ax[0].set_title("TESS-band flux ratio")
ax[1].errorbar(nt.SECTORS, [emp[s]["D_eff"] for s in nt.SECTORS], yerr=[emp[s]["D_eff_err"] for s in nt.SECTORS], fmt="o", label="empirical (QLP vs PSF pair)")
k0 = lambda s: [k for k in ap_info[s]["apertures"] if k.startswith("qlp")][0]
ax[1].plot(nt.SECTORS, [ap_info[s]["apertures"][k0(s)]["D_corr_gaia_r"] for s in nt.SECTORS], "s", mfc="none", label="PRF: QLP = target-normalised")
ax[1].plot(nt.SECTORS, [ap_info[s]["apertures"][k0(s)]["D_raw"] for s in nt.SECTORS], "^", mfc="none", label="PRF: QLP = raw aperture")
ax[1].axhline(4.35, color="k", ls=":", label="first-round 4.35")
ax[1].set_xlabel("sector"); ax[1].set_ylabel("dilution D for QLP SAP (neighbour host)"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/fluxratio_dilution.png"), dpi=110)
print("DONE")
