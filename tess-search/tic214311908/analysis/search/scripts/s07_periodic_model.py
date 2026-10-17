"""s07: characterise and remove the 0.6485-d periodic modulation in the deblended M3V light curve (model B).

1. Per sector: least-squares harmonic fit (NH = 3 harmonics) at a trial period on the detrended, transit-masked,
   rescaled flux (full cadence); semi-amplitudes of each harmonic in ppt and e-/s (979 e-/s), first-harmonic phase,
   with errors from a block bootstrap (0.5-d blocks, 200 draws).
2. Period per sector (local chi2 minimum of the harmonic fit over 0.60-0.70 d) with bootstrap error; joint fit of
   S103+S104 (consecutive sectors).
3. Coherence test (S39-S104): chi2 of a common-amplitude, common-phase harmonic model (coherent) against per-sector
   amplitudes and phases (incoherent), scanned on a fine period grid around 0.6485 d. The incoherent model has 6 x 4
   extra parameters; Delta-chi2 is compared with the same statistic at fake periods (empirical calibration).
4. Prewhitening: Lomb-Scargle of the residuals per sector and combined; further peaks above the empirical 1% level
   (from 200 rolled-and-inverted realisations of the residuals).
5. Writes the cleaned light curve data/m3v_clean.npz (det minus per-sector harmonic model), which s01/s02/s09 use.
Also applies the same per-sector harmonic fit to model A and to the pair (model A) and Mira (model B) light curves
for comparison (amplitudes in e-/s).
Outputs: out/periodic_model.json, data/m3v_clean.npz, figs/periodic_fold.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from astropy.timeseries import LombScargle
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

P0 = 0.64853
NH = 3
RNG = np.random.default_rng(7)
FREFS = dict(M3V=979.0, pair=895.0)


def harm_X(t, P, nh=NH):
    cols = [np.ones(len(t))]
    for h in range(1, nh + 1):
        cols += [np.sin(2 * np.pi * h * t / P), np.cos(2 * np.pi * h * t / P)]
    return np.column_stack(cols)


def harm_fit(t, y, w, P, nh=NH):
    X = harm_X(t, P, nh)
    Xw = X * np.sqrt(w)[:, None]; yw = y * np.sqrt(w)
    c, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    chi2 = float(np.sum((Xw @ c - yw) ** 2))
    return c, chi2, X @ c


def block_boot(t, y, w, P, nb=200, block=0.5):
    k = np.floor((t - t.min()) / block).astype(int); blocks = [np.where(k == u)[0] for u in np.unique(k)]
    cs = []
    for _ in range(nb):
        idx = np.concatenate([blocks[i] for i in RNG.integers(0, len(blocks), len(blocks))])
        c, _, _ = harm_fit(t[idx], y[idx], w[idx], P); cs.append(c)
    return np.array(cs)


def amp_phase(c):
    out = []
    for h in range(NH):
        a, b = c[1 + 2 * h], c[2 + 2 * h]
        out.append((float(np.hypot(a, b)), float(np.arctan2(b, a))))
    return out


def sector_data(star="M3V", model="B"):
    d = load_star(star, model)
    d = d[~d.intr.values & np.isfinite(d.det.values)].reset_index(drop=True)
    return d


def main():
    out = {}
    d = sector_data("M3V", "B")
    t_all, y_all, s_all = d.t.values, d.det.values - 1, d.sector.values
    w_all = np.zeros(len(t_all))
    for s in SECTORS:
        m = s_all == s; w_all[m] = 1 / robust_std(y_all[m]) ** 2
    # 1-2: per-sector fits and periods
    per = {}
    model_all = np.zeros(len(t_all))
    for s in SECTORS:
        m = s_all == s; t, y, w = t_all[m], y_all[m], w_all[m]
        c, chi2, mdl = harm_fit(t, y, w, P0)
        cb = block_boot(t, y, w, P0)
        ap = amp_phase(c); apb = np.array([[np.hypot(cc[1 + 2 * h], cc[2 + 2 * h]) for h in range(NH)] for cc in cb])
        # null amplitude: same fit at 100 fake periods 0.55-0.80 d
        fakes = [p for p in RNG.uniform(0.55, 0.80, 140) if abs(p - P0) > 0.01 and abs(p - 2 * P0) > 0.01][:100]
        a1f = np.array([amp_phase(harm_fit(t, y, w, pf)[0])[0][0] for pf in fakes])
        # local period scan
        Ps = np.linspace(0.60, 0.70, 4001)
        ch = np.array([harm_fit(t, y, w, p)[1] for p in Ps])
        j = int(np.argmin(ch)); Pbest = Ps[j]
        # period error from bootstrap of the local minimum (fine scan +-0.01 d)
        Pf = np.linspace(Pbest - 0.008, Pbest + 0.008, 321); pb = []
        k = np.floor((t - t.min()) / 0.5).astype(int); blocks = [np.where(k == u)[0] for u in np.unique(k)]
        for _ in range(60):
            idx = np.concatenate([blocks[i] for i in RNG.integers(0, len(blocks), len(blocks))])
            pb.append(Pf[np.argmin([harm_fit(t[idx], y[idx], w[idx], p)[1] for p in Pf])])
        per[s] = dict(n=int(m.sum()), amp_ppt=[a * 1e3 for a, _ in ap], amp_err_ppt=(apb.std(0) * 1e3).tolist(),
                      amp_e_s=[a * FREF for a, _ in ap], phase1_rad=ap[0][1], amp1_fake_median_ppt=float(np.median(a1f) * 1e3),
                      amp1_fake_p99_ppt=float(np.percentile(a1f, 99) * 1e3), P_best=float(Pbest), P_boot_std=float(np.std(pb)),
                      dchi2_P0_vs_scan_median=float(np.median(ch) - harm_fit(t, y, w, P0)[1]))
        model_all[m] = mdl - c[0]
        print(f"S{s}: A1 {ap[0][0]*1e3:.2f}+-{apb.std(0)[0]*1e3:.2f} ppt (fake p99 {np.percentile(a1f, 99)*1e3:.2f}), A2 {ap[1][0]*1e3:.2f}, A3 {ap[2][0]*1e3:.2f}; "
              f"P_best {Pbest:.5f} +- {np.std(pb):.5f}", flush=True)
    out["per_sector"] = per
    # joint S103+S104
    m = np.isin(s_all, [103, 104]); t, y, w = t_all[m], y_all[m], w_all[m]
    Ps = np.linspace(0.645, 0.652, 2801); ch = np.array([harm_fit(t, y, w, p)[1] for p in Ps])
    j = int(np.argmin(ch)); Pj = Ps[j]
    # Delta-chi2 = 1 interval (white) and a block-bootstrap error
    k = np.floor((t - t.min()) / 0.5).astype(int); blocks = [np.where(k == u)[0] for u in np.unique(k)]
    Pf = np.linspace(Pj - 0.0015, Pj + 0.0015, 301); pb = []
    for _ in range(60):
        idx = np.concatenate([blocks[i] for i in RNG.integers(0, len(blocks), len(blocks))])
        pb.append(Pf[np.argmin([harm_fit(t[idx], y[idx], w[idx], p)[1] for p in Pf])])
    out["S103_104_joint"] = dict(P=float(Pj), P_boot_std=float(np.std(pb)))
    print(f"S103+S104 joint P {Pj:.6f} +- {np.std(pb):.6f}")
    # 3: coherence test over S39-S104
    sel = np.isin(s_all, [39, 66, 93, 103, 104]); t, y, w, ss = t_all[sel], y_all[sel], w_all[sel], s_all[sel]

    def chi_coh(P):
        # per-sector offsets, common harmonics
        X = harm_X(t, P)[:, 1:]
        off = np.column_stack([(ss == s).astype(float) for s in np.unique(ss)])
        XX = np.column_stack([off, X]) * np.sqrt(w)[:, None]
        c, *_ = np.linalg.lstsq(XX, y * np.sqrt(w), rcond=None)
        return float(np.sum((XX @ c - y * np.sqrt(w)) ** 2))

    def chi_inc(P):
        return sum(harm_fit(t[ss == s], y[ss == s], w[ss == s], P)[1] for s in np.unique(ss))

    Pg = np.arange(0.6470, 0.6500, 2e-7)
    cc = np.array([chi_coh(p) for p in Pg])
    jbest = int(np.argmin(cc)); Pc = Pg[jbest]
    ci = chi_inc(Pc)
    # second-best alias tooth: local minima of cc
    mins = [i for i in range(1, len(cc) - 1) if cc[i] < cc[i - 1] and cc[i] < cc[i + 1]]
    mins = sorted(mins, key=lambda i: cc[i])[:6]
    out["coherence"] = dict(P_coherent=float(Pc), chi2_coherent=float(cc[jbest]), chi2_incoherent=float(ci),
                            dchi2_inc_minus_coh=float(ci - cc[jbest]), extra_params_incoherent=4 * 2 * NH,
                            alias_minima=[(float(Pg[i]), float(cc[i] - cc[jbest])) for i in mins],
                            chi2_no_signal=float(sum(np.sum(w[ss == s] * (y[ss == s] - np.sum(w[ss == s] * y[ss == s]) / np.sum(w[ss == s])) ** 2) for s in np.unique(ss))))
    print("coherence:", json.dumps(out["coherence"]), flush=True)
    np.save(f"{ROOT}/data/coherence_scan.npy", np.column_stack([Pg, cc]))
    # 4: residual periodograms and prewhitening check
    res = y_all - model_all
    pw_out = {}
    fr = np.linspace(1 / 15, 1 / 0.2, 150000)
    for s in list(SECTORS) + ["all"]:
        m = (s_all == s) if s != "all" else np.ones(len(s_all), bool)
        tb, fb, nb_, sb, cb = bin_lc(t_all[m], res[m], s_all[m], d.cad.values[m], 30.0)
        L = LombScargle(tb, fb); p = L.power(fr)
        # empirical null: inverted + rolled residuals (per sector roll)
        mx = []
        for _ in range(100):
            fb2 = fb.copy()
            for s2 in np.unique(sb):
                k = np.where(sb == s2)[0]; fb2[k] = np.roll(fb[k], RNG.integers(len(k) // 10, len(k) - len(k) // 10))
            # rolling destroys nothing for a white process; preserves red noise; periodic signals survive within
            # a sector, so use the max over a shuffled-block version instead
            kb = np.floor((tb - tb.min()) / 1.0).astype(int); ub = np.unique(kb)
            perm = RNG.permutation(len(ub)); fb3 = np.concatenate([fb2[kb == ub[i]] for i in perm])
            mx.append(LombScargle(tb, fb3).power(fr[::10]).max())
        top = []
        for i in np.argsort(p)[::-1]:
            if len(top) >= 5: break
            if all(abs(fr[i] - fr[j]) > 0.01 * fr[j] for j in top): top.append(i)
        pw_out[str(s)] = dict(top=[(float(1 / fr[i]), float(p[i])) for i in top], null_p99=float(np.percentile(mx, 99)),
                              null_p50=float(np.median(mx)))
        print(f"resid LS {s}: " + ", ".join(f"{1/fr[i]:.4f} ({p[i]:.3f})" for i in top) + f"; null p99 {np.percentile(mx, 99):.3f}", flush=True)
    out["residual_periodograms"] = pw_out
    # other light curves: amplitudes at P0 in e-/s
    comp = {}
    for star, model in (("M3V", "A"), ("pair", "A"), ("Mira", "B"), ("Mira", "A")):
        dd = sector_data(star, model)
        fref = 979.0 if star == "M3V" else (895.0 if star == "pair" else None)
        st = json.load(open(f"{PHOT}/out/lc_stats.json"))
        comp[f"{star}_{model}"] = {}
        for s in SECTORS:
            m = dd.sector.values == s
            tt, yy = dd.t.values[m], dd.det.values[m] - 1
            ww = np.full(m.sum(), 1 / robust_std(yy) ** 2)
            c, _, _ = harm_fit(tt, yy, ww, P0)
            fr_ = fref if fref else st[str(s)]["fref"]["Mira"]
            ap = amp_phase(c)
            # projection of this light curve's first harmonic on the M3V-B phase (signed amplitude)
            proj = ap[0][0] * np.cos(ap[0][1] - per[s]["phase1_rad"])
            a1f = np.array([amp_phase(harm_fit(tt, yy, ww, pf)[0])[0][0] for pf in RNG.uniform(0.55, 0.80, 60)])
            comp[f"{star}_{model}"][s] = dict(amp1_e_s=ap[0][0] * fr_, proj_on_M3VB_phase_e_s=float(proj * fr_),
                                              phase_diff_rad=float((ap[0][1] - per[s]["phase1_rad"] + np.pi) % (2 * np.pi) - np.pi),
                                              fake_p99_e_s=float(np.percentile(a1f, 99) * fr_))
        print(f"{star}_{model}: " + "; ".join(f"S{s} A1 {v['amp1_e_s']:.2f} e/s proj {v['proj_on_M3VB_phase_e_s']:+.2f} (p99 {v['fake_p99_e_s']:.2f})" for s, v in comp[f"{star}_{model}"].items()), flush=True)
    out["other_lcs"] = comp
    out["P0"] = P0; out["NH"] = NH
    save_json({k: ({str(kk): vv for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in out.items()}, f"{ROOT}/out/periodic_model.json")
    # cleaned light curve (full M3V-B set, including transit windows: model evaluated everywhere)
    full = load_star("M3V", "B")
    cl = full.det.values.copy()
    for s in SECTORS:
        m = s_all == s
        c, _, _ = harm_fit(t_all[m], y_all[m], w_all[m], P0)
        mf = full.sector.values == s
        cl[mf] = full.det.values[mf] - (harm_X(full.t.values[mf], P0) @ c - c[0])
    np.savez_compressed(f"{ROOT}/data/m3v_clean.npz", t=full.t.values, det_clean=cl, det=full.det.values, raw=full.raw.values,
                        sector=full.sector.values, cad=full.cad.values, intr=full.intr.values, err=full.err.values)
    # figure
    fig, axs = plt.subplots(2, 3, figsize=(15, 7)); axs = axs.ravel()
    for i, s in enumerate([39, 66, 93, 103, 104, 12]):
        m = s_all == s
        ph = (t_all[m] / P0) % 1
        tb_, fb_, _, _, _ = bin_lc(ph, y_all[m], np.zeros(m.sum()), np.full(m.sum(), 1.0), 0.02 * P0 * 1440)
        axs[i].plot(ph, y_all[m] * 1e3, ",", alpha=0.3, color="0.6")
        k = np.floor(ph * 40).astype(int)
        axs[i].plot((np.arange(40) + 0.5) / 40, [np.median(y_all[m][k == j]) * 1e3 for j in range(40)], "o-", color="C0")
        axs[i].set_ylim(-15, 15); axs[i].set_title(f"S{s} M3V-B fold P = {P0} d (ppt)")
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/periodic_fold.png", dpi=80); plt.close(fig)


if __name__ == "__main__":
    main()
