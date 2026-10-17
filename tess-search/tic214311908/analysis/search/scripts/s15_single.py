"""s15: single-transit (mono-transit) search of the cleaned M3V light curve, with an empirically calibrated threshold,
per-event diagnostics, and injection-recovery.

Data: data/m3v_clean.npz (model B, RRab leak removed), 10-min bins (S12 30 min); 5.69-d windows removed (and, for the
sanity check, kept). Matched filter: box of duration D in {0.75, 1, 1.5, 2, 3, 4, 6} h slid in steps of one bin; local
baseline = weighted mean of the data within [D, 3D] on both sides (removes residual trends); depth and S/N with the
per-sector white weights; box coverage >= 60 % and >= 2 baseline bins on each side required.
Threshold: the highest S/N of the same search on (a) the inverted light curve, (b) the 10 phantom-position and
4 control light curves (s03; RRab leak removed; both signs). Events are merged within +-0.5 d (best S/N kept).
Diagnostics for the top events: model-A M3V and pair (A) depths, Mira (B), background, and the whole-cutout PRF
fit of the event difference image (as for flares).
Injections: single trapezoids (T14 1-6 h, depth 3-60 ppt) added to the cleaned flux at random times -> recovery.
Output: out/single_transits.json, out/single_events.csv, out/single_inject.csv
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s07_periodic_model import harm_X, harm_fit
import numba as nb

DURS_H = np.array([0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
RNG = np.random.default_rng(15)


@nb.njit(cache=True)
def box_scan(t, y, w, cad_d, durs):
    n = len(t); nd = len(durs)
    snr = np.zeros((nd, n)); dep = np.zeros((nd, n))
    for k in range(nd):
        D = durs[k]
        for i in range(n):
            tc = t[i]
            sw = 0.0; swy = 0.0; nin = 0
            bw = 0.0; bwy = 0.0; nb1 = 0; nb2 = 0
            # scan neighbours
            j = i
            while j >= 0 and t[j] > tc - 2 * D:
                j -= 1
            j += 1
            while j < n and t[j] < tc + 2 * D:
                dt = t[j] - tc
                if abs(dt) <= D / 2:
                    sw += w[j]; swy += w[j] * y[j]; nin += 1
                elif abs(dt) >= D and abs(dt) <= 2 * D:
                    bw += w[j]; bwy += w[j] * y[j]
                    if dt < 0: nb1 += 1
                    else: nb2 += 1
                j += 1
            expected = D / cad_d[i]
            if nin < 0.6 * expected or nb1 < 2 or nb2 < 2 or sw == 0 or bw == 0: continue
            d = bwy / bw - swy / sw
            e = np.sqrt(1 / sw + 1 / bw)
            dep[k, i] = d; snr[k, i] = d / e
    return snr, dep


def prep(y_det, t, s, cad, intr, keep_transits=False):
    m = np.isfinite(y_det) & (keep_transits | ~intr)
    tb, fb, nb_, sb, cb = bin_lc(t[m], y_det[m], s[m], cad[m], 10.0)
    yy, ww = prep_weights(fb - 1, sb)
    return tb, yy, ww, sb, cb


def events(tb, snr, dep, durs, thr, merge=0.5):
    best = snr.max(0); kbest = snr.argmax(0)
    idx = np.where(best > thr)[0]
    idx = idx[np.argsort(best[idx])[::-1]]
    ev = []
    for i in idx:
        if any(abs(tb[i] - e["t"]) < merge for e in ev): continue
        ev.append(dict(t=float(tb[i]), snr=float(best[i]), depth_ppt=float(dep[kbest[i], i] * 1e3), dur_h=float(durs[kbest[i]] * 24)))
    return ev


def scan(tb, yy, ww, cb):
    durs = DURS_H / 24
    return box_scan(tb, yy, ww, cb / 86400, durs)


def clean_series(t, s, det):
    y = det - 1
    intr = transit_mask(t)
    for sec in np.unique(s):
        m = (s == sec) & ~intr & np.isfinite(y)
        c, _, _ = harm_fit(t[m], y[m], np.ones(m.sum()), 0.648528)
        mf = s == sec
        y[mf] -= harm_X(t[mf], 0.648528) @ c - c[0]
    return 1 + y


def main():
    z = np.load(f"{ROOT}/data/m3v_clean.npz")
    t, f, s, cad, intr = z["t"], z["det_clean"], z["sector"], z["cad"].astype(float), z["intr"]
    tb, yy, ww, sb, cb = prep(f, t, s, cad, intr)
    snr, dep = scan(tb, yy, ww, cb)
    snr_i, dep_i = scan(tb, -yy, ww, cb)
    max_inv = float(snr_i.max())
    # phantom / control light curves (null references)
    Z = np.load(f"{ROOT}/data/ctrl_lcs.npz")
    names = sorted({k.split("__")[0] for k in Z.files})
    ref = {}
    for nm in names:
        tt, dd, ss_, cc_ = Z[f"{nm}__t"], Z[f"{nm}__det"], Z[f"{nm}__sector"], Z[f"{nm}__cad"]
        if len(np.unique(ss_)) < 5: continue
        ff = clean_series(tt, ss_, dd.copy())
        a, b_, c_, d_, e_ = prep(ff, tt, ss_, cc_, transit_mask(tt))
        s1, _ = scan(a, b_, c_, e_); s2, _ = scan(a, -b_, c_, e_)
        ref[nm] = dict(max_pos=float(s1.max()), max_neg=float(s2.max()))
    ref_all = [max_inv] + [v["max_pos"] for v in ref.values()] + [v["max_neg"] for v in ref.values()]
    thr_inv = max_inv
    thr_med = float(np.median(ref_all))
    ev = events(tb, snr, dep, DURS_H / 24, min(thr_inv, thr_med) * 0.8)
    print(f"max S/N real {snr.max():.2f}; inverted {max_inv:.2f}; phantom/control maxima median {thr_med:.2f}, range {min(ref_all):.2f}-{max(ref_all):.2f}", flush=True)
    # sanity: the 5.69-d transits individually (unmasked)
    tb2, yy2, ww2, sb2, cb2 = prep(f, t, s, cad, intr, keep_transits=True)
    snr2, dep2 = scan(tb2, yy2, ww2, cb2)
    best2 = snr2.max(0)
    tr_snr = []
    for e in range(-400, 70):
        tc = T0_B + e * P_B
        k = np.abs(tb2 - tc) < 0.05
        if k.sum() >= 3: tr_snr.append((int(e), float(best2[k].max())))
    # diagnostics of the top events
    A = load_star("M3V", "A"); pa = load_star("pair", "A"); mi = load_star("Mira", "B")
    from s10_flares import pixel_attribution
    rows = []
    for e in ev[:12]:
        D = e["dur_h"] / 24
        row = dict(e)
        for nmx, dd, fref in (("M3V_A", A, 979.0), ("pair_A", pa, 895.0), ("Mira_B", mi, 1.0)):
            inn = np.abs(dd.t.values - e["t"]) <= D / 2; bb = (np.abs(dd.t.values - e["t"]) >= D) & (np.abs(dd.t.values - e["t"]) <= 2 * D)
            row[f"dF_{nmx}"] = float((np.nanmean(dd.det.values[bb]) - np.nanmean(dd.det.values[inn])) * fref) if inn.any() and bb.any() else np.nan
        sec = int(sb[np.argmin(np.abs(tb - e["t"]))]); row["sector"] = sec
        # pixel: treat the dip as a negative 'flare' (attribution function averages D over the event window)
        pix = pixel_attribution(dict(t0=e["t"] - D / 2, t1=e["t"] + D / 2), sec)
        if pix:
            row.update({f"pix_{k}": v for k, v in pix.items()})
        rows.append(row)
        print(f"event t {e['t']:.3f} S{sec} S/N {e['snr']:.2f} depth {e['depth_ppt']:.1f} ppt dur {e['dur_h']:.1f} h; A:M3V {row['dF_M3V_A']:.1f} A:pair {row['dF_pair_A']:.1f} e/s; "
              + (f"pix best dEN {pix['best_dEN_from_M3V'][0]:+.1f},{pix['best_dEN_from_M3V'][1]:+.1f}\" flux {pix['flux_best_e_s']:.1f} dchi2(M3V) {pix['dchi2_M3V']:.1f}" if pix else ""), flush=True)
    pd.DataFrame(rows).to_csv(f"{ROOT}/out/single_events.csv", index=False)
    # injections (cleaned flux, no re-detrending: the 0.75-d biweight would partly absorb the longest boxes, noted)
    inj = []
    good_t = tb[np.isfinite(yy)]
    for k in range(400):
        T14 = np.exp(RNG.uniform(np.log(1.0), np.log(6.0))) / 24
        depth = np.exp(RNG.uniform(np.log(3e-3), np.log(6e-2)))
        tc = float(RNG.choice(good_t))
        x = np.abs(tb - tc); ing = 0.15 * T14
        mdl = np.where(x <= T14 / 2 - ing, 1.0, np.where(x < T14 / 2, (T14 / 2 - x) / ing, 0.0)) * depth
        y2 = yy - mdl
        lo = np.abs(tb - tc) < 0.6
        s2, d2 = scan(tb[lo], y2[lo], ww[lo], cb[lo])
        inj.append(dict(T14_h=T14 * 24, depth_ppt=depth * 1e3, sector=int(sb[np.argmin(np.abs(tb - tc))]), snr=float(s2.max())))
    I = pd.DataFrame(inj); I.to_csv(f"{ROOT}/out/single_inject.csv", index=False)
    thr = max(thr_inv, thr_med)
    def rec(dmin, dmax, Tmin, Tmax):
        k = (I.depth_ppt >= dmin) & (I.depth_ppt < dmax) & (I.T14_h >= Tmin) & (I.T14_h < Tmax)
        return (int(k.sum()), float((I.snr[k] > thr).mean()) if k.sum() else np.nan)
    rec_tab = {f"{a}-{b} ppt, {c}-{d} h": rec(a, b, c, d) for a, b in ((3, 6), (6, 10), (10, 15), (15, 25), (25, 60)) for c, d in ((1, 2.5), (2.5, 6))}
    out = dict(max_snr_real=float(snr.max()), max_snr_inverted=max_inv, ref_maxima=ref, threshold_inverted=thr_inv,
               threshold_ref_median=thr_med, threshold_adopted=thr, events_above_0p8thr=ev[:20],
               n_events_above_adopted=int(sum(e["snr"] > thr for e in ev)), known_transit_single_snr=tr_snr,
               n_known_transits_above_adopted=int(sum(v > thr for _, v in tr_snr)), injection_recovery=rec_tab)
    save_json(out, f"{ROOT}/out/single_transits.json")
    print(json.dumps({k: v for k, v in out.items() if k not in ("ref_maxima", "events_above_0p8thr")}, indent=0))


if __name__ == "__main__":
    main()
