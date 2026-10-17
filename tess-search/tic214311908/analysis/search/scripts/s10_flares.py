"""s10: flares on the deblended M3V light curve with an injection-calibrated detector, attribution of candidates
(M3V vs the M-dwarf pair, which leaks ~40% into model B), and the flare rate against M dwarfs of known activity/age.

Light curve: model B, quality_ok, detrended, rescaled to 979 e-/s, with the RRab leak (s07) removed; per-cadence
errors = phot flux_err x 0.955 (sector robust sigma x local chi2 factor). z = (f - 1) / err.
Detector: runs of consecutive cadences with z > 2.5 (gaps < 1.5 cadences); statistic Zrun = sum(z)/sqrt(n) and peak z.
Candidate if peak z >= 3.5 and Zrun >= Z_thr. Z_thr is calibrated on the inverted light curve (same detector on -f):
the smallest threshold giving <= 1 false event in the whole data set.
Injections: Davenport (2014) flare template (t_1/2 log-uniform 2-40 min, amplitude log-uniform 0.3-60 %), integrated
over each exposure, added to the undetrended flux; the chain (biweight 0.75 d with 5.69-d mask, RRab removal) is
re-run; ~45 flares per realisation at random times, 30 realisations. Recovery if a detected run overlaps
[t_peak - 1 cadence, t_peak + 2 t_1/2 + 1 cadence]. Completeness vs equivalent duration (ED).
Energies: E_TESS = ED x L_TESS, L_TESS = 1.4e31 erg/s (Gaia XP spectrum x TESS response, peak-normalised; d = 64.31 pc;
the Sullivan+15 zero point gives 1.16e31, a 3350-K blackbody 1.73e31: +-20 %).
Attribution: for every candidate, amplitude in the model-A M3V and pair light curves, and a whole-cutout PRF fit of the
flare difference image (mean D over the event minus the mean of +-0.15 d around it), as in s05b.
Comparison: Medina et al. 2020 (ApJ 905, 107): 0.1-0.3 M_sun within 15 pc, alpha = 1.98 +- 0.02;
active (H-alpha emission, Ro < 0.5) log R31.5 = -1.30 +- 0.08 flares/day above E_TESS = 3.16e31 erg; inactive < -3.86.
Output: out/flares.json, out/flare_candidates.csv, out/flare_injections.csv, figs/flares.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s07_periodic_model import harm_X, harm_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RNG = np.random.default_rng(41)
P_RR = 0.648528
L_TESS = 1.4e31
ZLO = 2.5


def davenport(x):
    """Davenport et al. 2014 template; x = (t - t_peak) / t_half."""
    y = np.zeros_like(x)
    r = (x >= -1) & (x <= 0)
    y[r] = 1 + 1.941 * x[r] - 0.175 * x[r] ** 2 - 2.246 * x[r] ** 3 - 1.125 * x[r] ** 4
    d = x > 0
    y[d] = 0.6890 * np.exp(-1.600 * x[d]) + 0.3030 * np.exp(-0.2783 * x[d])
    return np.clip(y, 0, None)


def flare_model(t, cad_s, tp, th, amp):
    out = np.zeros(len(t))
    near = np.abs(t - tp) < (2 * th + 30 * th + 0.05)
    if not near.any(): return out
    tt = t[near]; cc = cad_s[near] / 86400
    ss = np.linspace(-0.5, 0.5, 11)
    acc = np.zeros(near.sum())
    for s_ in ss:
        acc += davenport((tt + s_ * cc - tp) / th)
    out[near] = amp * acc / len(ss)
    return out


def load():
    d = load_star("M3V", "B")
    A = load_star("M3V", "A")
    pa = load_star("pair", "A")
    return d, A, pa


def clean_chain(d, raw):
    t = d.t.values; s = d.sector.values; intr = d.intr.values
    y = np.empty(len(t))
    for sec in SECTORS:
        m = s == sec
        y[m] = raw[m] / detrend(t[m], raw[m], intr[m]) - 1
    for sec in SECTORS:
        m = (s == sec) & ~intr
        c, _, _ = harm_fit(t[m], y[m], np.ones(m.sum()), P_RR)
        mf = s == sec
        y[mf] -= harm_X(t[mf], P_RR) @ c - c[0]
    return y


def runs(t, z, cad_s, zlo=ZLO):
    """Runs of consecutive cadences with z > zlo; returns list of (i0, i1) inclusive."""
    hi = z > zlo
    out = []; i = 0; n = len(z)
    while i < n:
        if hi[i]:
            j = i
            while j + 1 < n and hi[j + 1] and (t[j + 1] - t[j]) < 1.5 * cad_s[j] / 86400: j += 1
            out.append((i, j)); i = j + 1
        else: i += 1
    return out


def detect(t, y, err, cad_s, zthr, zpeak=3.5):
    z = y / err
    ev = []
    for i0, i1 in runs(t, z, cad_s):
        zz = z[i0:i1 + 1]; n = len(zz)
        Zr = zz.sum() / np.sqrt(n)
        if zz.max() >= zpeak and Zr >= zthr:
            ed = float(np.sum(y[max(i0 - 1, 0):i1 + 2] * cad_s[max(i0 - 1, 0):i1 + 2]))    # seconds
            ev.append(dict(i0=int(i0), i1=int(i1), t0=float(t[i0]), t1=float(t[i1]), n=int(n), zpeak=float(zz.max()), Zrun=float(Zr),
                           amp=float(y[i0:i1 + 1].max()), ED_s=ed))
    return ev


def run_stats(t, y, err, cad_s):
    z = y / err
    return [(float(z[i0:i1 + 1].sum() / np.sqrt(i1 - i0 + 1)), float(z[i0:i1 + 1].max())) for i0, i1 in runs(t, z, cad_s)]


def pixel_attribution(ev, sector):
    """Whole-cutout single-PRF fit of the flare difference image (s05b machinery)."""
    from s05b_pixfull_periodic import PRFs, search
    z = np.load(f"{ROOT}/data/pixfull_s{sector}.npz", allow_pickle=True)
    t = z["t_tdb"]; ok = z["ok"].astype(bool) & np.isfinite(t)
    D = z["D"].astype(float)
    inn = ok & (t >= ev["t0"] - 1e-4) & (t <= ev["t1"] + 1e-4)
    out_ = ok & (np.abs(t - 0.5 * (ev["t0"] + ev["t1"])) < 0.15) & ~((t >= ev["t0"] - 0.02) & (t <= ev["t1"] + 0.05))
    if inn.sum() == 0 or out_.sum() < 5: return None
    img = np.nanmean(D[inn], 0) - np.nanmean(D[out_], 0)
    # per-pixel error from the scatter of single-cadence differences in the out-of-event set, scaled to n_in
    E = np.nanstd(D[out_] - np.nanmean(D[out_], 0), 0) * np.sqrt(1 / inn.sum() + 1 / out_.sum()) + 1e-3
    prfs = PRFs(sector)
    m = np.ones(img.shape, bool)
    b0, grid = search(img[None], E[None], prfs, m, (0.0, img.shape[1] - 1.0), (0.0, img.shape[0] - 1.0), 0.25)
    x, y = b0[1]
    b1, _ = search(img[None], E[None], prfs, m, (x - 0.3, x + 0.3), (y - 0.3, y + 0.3), 0.03)
    m3 = z["m3v_xy"]; pa = z["pair_xy"]; jac = z["jac"]
    from s05b_pixfull_periodic import fit_all
    c3, F3 = fit_all(img[None], E[None], prfs.img(*m3), m)
    cp, Fp = fit_all(img[None], E[None], prfs.img(*pa), m)
    return dict(best_dEN_from_M3V=(jac @ (np.array(b1[1]) - m3)).tolist(), flux_best_e_s=float(b1[2][0]), chi2_best=b1[0],
                dchi2_M3V=c3 - b1[0], dchi2_pair=cp - b1[0], flux_M3V_e_s=float(F3[0]), flux_pair_e_s=float(Fp[0]),
                pair_dEN_from_M3V=(jac @ (np.array(pa) - m3)).tolist())


def main():
    d, A, pa = load()
    t = d.t.values; cad = d.cad.values.astype(float); err = d.err.values; raw = d.raw.values
    y = clean_chain(d, raw)
    # ---- threshold from the inverted light curve
    st_inv = np.array(run_stats(t, -y, err, cad))
    st_real = np.array(run_stats(t, y, err, cad))
    grid = np.arange(3.0, 12.01, 0.1)
    n_inv = np.array([np.sum((st_inv[:, 0] >= g) & (st_inv[:, 1] >= 3.5)) for g in grid])
    n_real = np.array([np.sum((st_real[:, 0] >= g) & (st_real[:, 1] >= 3.5)) for g in grid])
    zthr = float(grid[np.argmax(n_inv <= 1)])
    print("threshold Zrun", zthr, "inverted events at threshold", n_inv[np.argmax(n_inv <= 1)], flush=True)
    ev = detect(t, y, err, cad, zthr)
    # ---- injections
    rows = []
    for k in range(30):
        tps, ths, amps, rawi = [], [], [], raw.copy()
        sel = RNG.choice(len(t), 45, replace=False)
        for i in sel:
            tp = t[i] + RNG.uniform(-0.01, 0.01); th = np.exp(RNG.uniform(np.log(2), np.log(40))) / 1440
            amp = np.exp(RNG.uniform(np.log(0.003), np.log(0.6)))
            rawi += flare_model(t, cad, tp, th, amp); tps.append(tp); ths.append(th); amps.append(amp)
        yi = clean_chain(d, rawi)
        evi = detect(t, yi, err, cad, zthr)
        for tp, th, amp in zip(tps, ths, amps):
            i = int(np.argmin(np.abs(t - tp))); c_ = cad[i] / 86400
            if abs(t[i] - tp) > 2 * c_: continue          # flare fell in a gap
            # injected ED (s): exposure-integrated template summed over cadences
            mdl = flare_model(t, cad, tp, th, amp)
            ed_in = float(np.sum(mdl * cad))
            hit = [e for e in evi if e["t1"] >= tp - c_ and e["t0"] <= tp + 2 * th + c_]
            rows.append(dict(tp=tp, t_half_min=th * 1440, amp=amp, ED_in_s=ed_in, sector=int(d.sector.values[i]), cad_s=cad[i],
                             rec=bool(hit), ED_rec_s=float(hit[0]["ED_s"]) if hit else np.nan))
        if k % 5 == 0: print("inj", k, len(rows), np.mean([r["rec"] for r in rows]), flush=True)
    R = pd.DataFrame(rows); R.to_csv(f"{ROOT}/out/flare_injections.csv", index=False)
    # completeness vs ED (all, and per cadence class)
    bins = np.logspace(-0.5, 3, 15)
    comp = {}
    for nm, m in (("all", np.ones(len(R), bool)), ("200s", R.cad_s < 300), ("600s", (R.cad_s > 300) & (R.cad_s < 1000)), ("1800s", R.cad_s > 1000)):
        r_ = R[m]; c = []
        for a, b in zip(bins[:-1], bins[1:]):
            k = (r_.ED_in_s >= a) & (r_.ED_in_s < b)
            c.append((float(np.sqrt(a * b)), int(k.sum()), float(r_.rec[k].mean()) if k.sum() else np.nan))
        comp[nm] = c
    def ed_at(frac, nm="all"):
        c = [x for x in comp[nm] if np.isfinite(x[2])]
        for e, n, f in c:
            if f >= frac: return e
        return np.nan
    ED50, ED90 = ed_at(0.5), ed_at(0.9)
    # false-positive check with injections: detections not matching an injection are not tracked; inverted rate used
    # ---- candidates: attribution
    MA = dict(zip(np.round(A.t.values, 6), A.det.values - 1)); PA = dict(zip(np.round(pa.t.values, 6), pa.det.values - 1))
    cand = []
    for e in ev:
        tt = t[e["i0"]:e["i1"] + 1]
        aA = np.nanmean([MA.get(round(x, 6), np.nan) for x in tt]) * 979
        aP = np.nanmean([PA.get(round(x, 6), np.nan) for x in tt]) * 895
        aB = float(np.mean(y[e["i0"]:e["i1"] + 1]) * 979)
        sec = int(d.sector.values[e["i0"]])
        pix = pixel_attribution(e, sec)
        c = dict(e, sector=sec, mean_dF_B_e_s=aB, mean_dF_M3V_A_e_s=float(aA), mean_dF_pair_A_e_s=float(aP),
                 E_TESS_erg=e["ED_s"] * L_TESS, **({f"pix_{k}": v for k, v in pix.items()} if pix else {}))
        cand.append(c)
        print(f"cand S{sec} t {e['t0']:.4f} n {e['n']} zpeak {e['zpeak']:.1f} Zrun {e['Zrun']:.1f} amp {e['amp']*1e3:.0f} ppt ED {e['ED_s']:.1f} s; "
              f"B {aB:.1f} A:M3V {aA:.1f} A:pair {aP:.1f} e/s; pix " + (f"dEN {pix['best_dEN_from_M3V'][0]:+.1f},{pix['best_dEN_from_M3V'][1]:+.1f}\" dchi2 M3V {pix['dchi2_M3V']:.1f} pair {pix['dchi2_pair']:.1f}" if pix else "n/a"), flush=True)
    pd.DataFrame(cand).to_csv(f"{ROOT}/out/flare_candidates.csv", index=False)
    # ---- rates
    days = float(np.sum(cad) / 86400)
    on_m3v = [c for c in cand if c.get("pix_dchi2_M3V", 99) < 6.2 and c.get("pix_dchi2_pair", 0) > 0]
    n_above = sum(1 for c in on_m3v if c["ED_s"] >= ED50) if np.isfinite(ED50) else 0
    E50 = ED50 * L_TESS
    # Medina+20 expectation above E50 (alpha = 1.98): R(>E) = R31.5 (E / 3.16e31)^(1 - alpha)
    exp_active = 10 ** -1.30 * (E50 / 3.16e31) ** (1 - 1.98) * days
    exp_active_lo = 10 ** -1.38 * (E50 / 3.16e31) ** (1 - 1.98) * days
    exp_inactive_max = 10 ** -3.86 * (E50 / 3.16e31) ** (1 - 1.98) * days
    from scipy.stats import poisson
    res = dict(days=days, Zrun_threshold=zthr, n_inverted_at_threshold=int(n_inv[np.argmax(n_inv <= 1)]),
               n_candidates=len(cand), n_on_M3V_pixels=len(on_m3v), completeness=comp, ED50_s=ED50, ED90_s=ED90,
               E50_TESS_erg=E50, E90_TESS_erg=ED90 * L_TESS, L_TESS_erg_s=L_TESS,
               n_M3V_flares_above_ED50=n_above,
               rate_above_E50_per_day=n_above / days / 0.75 if n_above else 0.0,
               rate_above_E50_upper95_per_day=float(poisson.ppf(0.95, max(n_above, 0) + 0.5) if False else (-np.log(0.05) if n_above == 0 else np.nan)) / days / 0.75,
               medina20_active_expected_N_above_E50=float(exp_active), medina20_active_expected_N_lo=float(exp_active_lo),
               medina20_inactive_expected_N_above_E50_max=float(exp_inactive_max),
               P_le_observed_if_active=float(poisson.cdf(n_above, exp_active)),
               note="completeness factor 0.75 = mean recovery above ED50 (approx.); rates in the TESS band")
    print(json.dumps({k: v for k, v in res.items() if k != "completeness"}, indent=0), flush=True)
    save_json(res, f"{ROOT}/out/flares.json")
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for nm in ("200s", "600s", "1800s", "all"):
        c = np.array([x for x in comp[nm] if np.isfinite(x[2])])
        if len(c): ax[0].semilogx(c[:, 0], c[:, 2], "o-", label=nm)
    ax[0].set_xlabel("injected ED (s)"); ax[0].set_ylabel("recovery fraction"); ax[0].legend()
    ax[1].semilogx(grid, n_inv, label="inverted (false)"); ax[1].semilogx(grid, n_real, label="real"); ax[1].axvline(zthr, color="k", ls=":")
    ax[1].set_xlabel("Zrun threshold"); ax[1].set_ylabel("N events (peak z >= 3.5)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/flares.png", dpi=80); plt.close(fig)


if __name__ == "__main__":
    main()
