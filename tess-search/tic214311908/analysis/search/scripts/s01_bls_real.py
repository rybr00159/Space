"""s01: coherent BLS of the deblended M3V light curve (all six sectors, 5.69-d transits masked) and per-sector BLS.

Data: phot/lc_M3V_all.csv model B, quality_ok, detrended flux rescaled to 979 e-/s; 5.69-d windows +-0.1 d removed;
binned to 10 min within each sector (S12 stays at 30 min); per-sector white weights from the robust std.
Search: box (BLS) Delta-chi2 over 0.3-30 d on a duration-aware frequency grid (phase drift over the 2574-d baseline
< 0.225 x central T14), five durations 0.3-1.25 x the central-transit T14 for rho = 9.4 g/cm3 (min 0.35 h).
Per sector (and S103+S104 together): 0.3 d to half the span.
Also: model A light curve (all-free deblending) for comparison.
Outputs: out/bls_real.json (top peaks with per-sector / per-event diagnostics), data/bls_real_{tag}.npz (periodograms).
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *


def search_data(model="B", flux="det", invert=False):
    """flux='det': phot detrended flux; flux='clean': det minus the per-sector 0.6485-d harmonic model (s07)."""
    if flux == "clean":
        z = np.load(f"{ROOT}/data/m3v_clean.npz")
        d = pd.DataFrame(dict(t=z["t"], clean=z["det_clean"], sector=z["sector"], cad=z["cad"], intr=z["intr"]))
    else:
        d = load_star("M3V", model)
    m = ~d.intr.values
    tb, fb, nb_, sb, cb = bin_lc(d.t.values[m], d[flux].values[m], d.sector.values[m], d.cad.values[m], 10.0)
    y, w = prep_weights(fb - 1, sb)
    if invert: y = -y
    return tb, y, w, sb, cb


def top_peaks(r, n=15, excl_frac=0.002):
    """Distinct peaks: greedy, excluding +-excl_frac in frequency around each accepted peak and its harmonics."""
    f = r["freq"]; p = r["dchi2"].copy(); out = []
    order = np.argsort(p)[::-1]
    used = []
    for i in order:
        if len(out) >= n: break
        ok = True
        for fu in used:
            for h in (1/3, 0.5, 1, 2, 3):
                if abs(f[i] - h * fu) < excl_frac * fu * max(h, 1): ok = False; break
            if not ok: break
        if ok:
            used.append(f[i]); out.append(int(i))
    return out


def per_event(t, y, w, sb, P, T0, dur):
    """Delta-chi2 contribution of each individual event and each sector for a box at (P, T0, dur)."""
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    inb = np.abs(ph * P) < dur / 2
    ep = np.round((t - T0) / P).astype(int)
    W = w.sum(); a = w[inb].sum(); c = (w * y)[inb].sum()
    dep = -c / a / (1 - a / W)
    ev = {}
    for e in np.unique(ep[inb]):
        k = inb & (ep == e)
        # contribution of this event to the global box fit: chi2 change with the global depth
        ev[int(e)] = dict(sector=int(sb[k][0]), n=int(k.sum()), depth=float(-(w[k] * y[k]).sum() / w[k].sum()),
                          err=float(1 / np.sqrt(w[k].sum())), dchi2_own=float(((w[k] * y[k]).sum()) ** 2 / w[k].sum()))
    sec = {}
    for s in np.unique(sb):
        k = inb & (sb == s)
        if k.sum() == 0: continue
        sec[int(s)] = dict(n_events=int(len(np.unique(ep[k]))), depth=float(-(w[k] * y[k]).sum() / w[k].sum()),
                           err=float(1 / np.sqrt(w[k].sum())))
    own = np.array([v["dchi2_own"] * np.sign(v["depth"]) for v in ev.values()])
    return dict(depth=float(dep), depth_err=float(1 / np.sqrt(a * (1 - a / W))), n_events=len(ev),
                max_event_frac=float(own.max() / max(own[own > 0].sum(), 1e-9)) if len(own) else np.nan,
                n_events_pos=int((own > 0).sum()), events=ev, sectors=sec)


def run(tag, t, y, w, sb, pmin, pmax, B=None):
    B = B or (t.max() - t.min())
    fr = freq_grid(pmin, pmax, B)
    r = run_bls(t, y, w, fr)
    np.savez_compressed(f"{ROOT}/data/bls_real_{tag}.npz", **r)
    pk = top_peaks(r)
    res = []
    for i in pk:
        P = 1 / r["freq"][i]; dur = DUR_FACS[r["kdur"][i]] * t14_central(P); dur = max(dur, 0.35 / 24)
        pe = per_event(t, y, w, sb, P, r["t0"][i], dur)
        res.append(dict(P=float(P), T0=float(r["t0"][i]), dur_h=float(dur * 24), dchi2=float(r["dchi2"][i]), snr=float(np.sqrt(r["dchi2"][i])),
                        depth_ppt=float(r["depth"][i] * 1e3), sde=float((r["dchi2"][i] - r["dchi2"].mean()) / r["dchi2"].std()),
                        n_events=pe["n_events"], n_events_pos=pe["n_events_pos"], max_event_frac=pe["max_event_frac"],
                        sectors={k: dict(depth_ppt=v["depth"] * 1e3, err_ppt=v["err"] * 1e3, n=v["n_events"]) for k, v in pe["sectors"].items()}))
    print(f"[{tag}] nfreq {len(fr)}; top: " + "; ".join(f"P {x['P']:.4f} d2 {x['dchi2']:.1f} dep {x['depth_ppt']:.2f} ppt dur {x['dur_h']:.2f} h nev {x['n_events']} maxfrac {x['max_event_frac']:.2f}" for x in res[:6]), flush=True)
    return dict(nfreq=len(fr), baseline=B, pmin=pmin, pmax=pmax, peaks=res,
                dchi2_mean=float(r["dchi2"].mean()), dchi2_std=float(r["dchi2"].std()))


def main_clean():
    """Search of the cleaned light curve (0.6485-d modulation removed): coherent + per sector, real and inverted."""
    out = {}
    t, y, w, sb, cb = search_data("B", "clean")
    out["all_clean"] = run("all_clean", t, y, w, sb, 0.3, 30.0)
    out["all_clean_inverted"] = run("all_clean_inv", t, -y, w, sb, 0.3, 30.0)
    groups = {"S12": [12], "S39": [39], "S66": [66], "S93": [93], "S103": [103], "S104": [104], "S103_104": [103, 104]}
    for g, ss in groups.items():
        m = np.isin(sb, ss)
        span = t[m].max() - t[m].min()
        yy = y[m] - np.sum(w[m] * y[m]) / w[m].sum()
        out[g] = run(g + "_clean", t[m], yy, w[m], sb[m], 0.3, min(30.0, span / 2))
    save_json(out, f"{ROOT}/out/bls_real_clean.json")


def main():
    out = {}
    t, y, w, sb, cb = search_data("B")
    out["all_B"] = run("all_B", t, y, w, sb, 0.3, 30.0)
    tA, yA, wA, sA, _ = search_data("A")
    out["all_A"] = run("all_A", tA, yA, wA, sA, 0.3, 30.0)
    out["all_B_inverted"] = run("all_B_inv", t, -y, w, sb, 0.3, 30.0)
    groups = {"S12": [12], "S39": [39], "S66": [66], "S93": [93], "S103": [103], "S104": [104], "S103_104": [103, 104]}
    for g, ss in groups.items():
        m = np.isin(sb, ss)
        span = t[m].max() - t[m].min()
        out[g] = run(g, t[m], y[m] - np.sum(w[m] * y[m]) / w[m].sum(), w[m], sb[m], 0.3, min(30.0, span / 2))
        out[g + "_inv"] = run(g + "_inv", t[m], -(y[m] - np.sum(w[m] * y[m]) / w[m].sum()), w[m], sb[m], 0.3, min(30.0, span / 2))
    save_json(out, f"{ROOT}/out/bls_real.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "clean": main_clean()
    else: main()
