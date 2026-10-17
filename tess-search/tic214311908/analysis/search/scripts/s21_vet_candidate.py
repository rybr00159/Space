"""s21: per-event vetting of a candidate periodic signal (default: the S103+S104 TLS/BLS peak P = 16.79 d).

For each predicted event inside the data: depth (box of the given duration, local linear baseline within +-3 durations)
in the cleaned model-B M3V light curve, the model-A M3V, model-A pair, model-B Mira and the background; the TESS
quality flags and background level in the window; and the whole-cutout PRF fit of the event difference image
(as s10/s15). Empirical per-event error: the same box at 200 random times in the same sector.
Usage: s21_vet_candidate.py P T0 dur_h tag     Output: out/vet_<tag>.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from s10_flares import pixel_attribution
RNG = np.random.default_rng(21)


def box_depth(t, y, tc, D):
    inn = np.abs(t - tc) <= D / 2
    bb = (np.abs(t - tc) > D) & (np.abs(t - tc) <= 3 * D)
    if inn.sum() < 2 or bb.sum() < 6: return np.nan
    X = np.column_stack([np.ones(bb.sum()), t[bb] - tc]); c, *_ = np.linalg.lstsq(X, y[bb], rcond=None)
    return float(c[0] - np.mean(y[inn] - c[1] * (t[inn] - tc)))


def main():
    P, T0, dur_h, tag = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    D = dur_h / 24
    z = np.load(f"{ROOT}/data/m3v_clean.npz")
    t, f, s = z["t"], z["det_clean"], z["sector"]
    full = pd.read_csv(f"{PHOT}/lc_M3V_all.csv")
    A = load_star("M3V", "A"); pa = load_star("pair", "A"); mi = load_star("Mira", "B")
    ev = []
    for e in range(-1000, 1000):
        tc = T0 + e * P
        if not np.any(np.abs(t - tc) < D / 2): continue
        sec = int(s[np.argmin(np.abs(t - tc))])
        m = s == sec
        dB = box_depth(t[m], f[m] - 1, tc, D) * 979
        fk = [box_depth(t[m], f[m] - 1, tt, D) * 979 for tt in RNG.uniform(t[m].min() + 3 * D, t[m].max() - 3 * D, 200)]
        err = float(np.nanstd(fk))
        row = dict(epoch=e, tc=tc, sector=sec, dF_M3V_B_clean_e_s=dB, err_emp_e_s=err, snr=dB / err)
        for nm, dd, fr in (("M3V_A", A, 979.0), ("pair_A", pa, 895.0)):
            mm = dd.sector.values == sec
            row[f"dF_{nm}_e_s"] = box_depth(dd.t.values[mm], dd.det.values[mm] - 1, tc, D) * fr
        mm = mi.sector.values == sec
        row["dF_Mira_B_rel_ppt"] = box_depth(mi.t.values[mm], mi.det.values[mm] - 1, tc, D) * 1e3
        w = (full.time_btjd_tdb.values > tc - 2 * D) & (full.time_btjd_tdb.values < tc + 2 * D)
        row["tess_quality_flags"] = sorted(set(int(q) for q in full.tess_quality.values[w]))
        row["bkg_mean_e_s_px"] = float(np.nanmean(full.bkg_e_s_px.values[w])) if w.any() else np.nan
        row["bkg_sector_median"] = float(np.nanmedian(full.bkg_e_s_px.values[full.sector.values == sec]))
        row["frac_quality_ok"] = float(np.mean(full.quality_ok.values[w])) if w.any() else np.nan
        pix = pixel_attribution(dict(t0=tc - D / 2, t1=tc + D / 2), sec)
        if pix:
            # flip sign: a dip gives a positive 'flux' in the out-minus-in sense used by pixel_attribution? (it uses in - out)
            row.update({f"pix_{k}": v for k, v in pix.items()})
        ev.append(row)
        print(f"E{e} tc {tc:.3f} S{sec}: B {dB:.1f}+-{err:.1f} e/s; A:M3V {row['dF_M3V_A_e_s']:.1f} A:pair {row['dF_pair_A_e_s']:.1f}; Mira {row['dF_Mira_B_rel_ppt']:.2f} ppt; "
              f"bkg {row['bkg_mean_e_s_px']:.0f} (med {row['bkg_sector_median']:.0f}); q {row['tess_quality_flags']} ok {row['frac_quality_ok']:.2f}; "
              + (f"pix best dEN {pix['best_dEN_from_M3V'][0]:+.1f},{pix['best_dEN_from_M3V'][1]:+.1f}\" flux {pix['flux_best_e_s']:.1f} (M3V {pix['flux_M3V_e_s']:.1f}) dchi2(M3V) {pix['dchi2_M3V']:.1f}" if pix else ""), flush=True)
    save_json(dict(P=P, T0=T0, dur_h=dur_h, events=ev), f"{ROOT}/out/vet_{tag}.json")


if __name__ == "__main__":
    main()
