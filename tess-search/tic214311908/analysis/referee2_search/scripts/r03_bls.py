"""r03: independent coherent box search for additional transiting planets (0.3-30 d, all six sectors) on the
referee-2 RRab-cleaned M3V light curve (r01: data/m3v_clean_r2.npz), 5.69-d windows (+-0.1 d) removed, 10-min bins
(S12 30 min), per-sector robust white weights; bls2 statistic and grid.

Usage:
  r03_bls.py real                 -> out/r03_real.json (top distinct peaks with per-event fractions, band maxima;
                                     also the inverted light curve), data/r03_real_pgram.npz (decimated)
  r03_bls.py null <kind> <n> <seed> -> out/r03_null_<kind>_<seed>.json
     kind = blk  : flux in 1-d blocks shuffled within each sector (time stamps, gaps, weights fixed);
                   destroys periodicity > ~1 d and inter-sector coherence, keeps single events and red noise
     kind = blki : same, then inverted
     kind = shift: each sector's flux kept intact but its time stamps shifted by a random offset 0-60 d
                   (no wrap junction; keeps within-sector periodicity, destroys inter-sector coherence)
Bands for maxima: 0.3-1, 1-3, 3-10, 10-30 d.
"""
import sys, os, time
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
import bls2

BANDS = [(0.3, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 30.0)]


def data():
    z = np.load(f"{ME}/data/m3v_clean_r2.npz")
    m = ~z["intr"]
    tb, yb, nb, sb = binlc(z["t"][m], z["clean"][m] - 1, z["s"][m], z["cad"][m], 10.0)
    y, w = bls2.weights(yb, sb)
    return tb, y, w, sb


def band_max(P, d2):
    return {f"{a}-{b}": float(d2[(P >= a) & (P < b)].max()) for a, b in BANDS}


def peaks(r, n=10):
    f = r["freq"]; d2 = r["dchi2"]; used = []; out = []
    for i in np.argsort(d2)[::-1]:
        if len(out) >= n: break
        if any(abs(f[i] - h * u) < 0.003 * u * max(h, 1) for u in used for h in (1 / 3, 0.5, 1, 2, 3)):
            continue
        used.append(f[i]); out.append(int(i))
    return out


def per_event(t, y, w, s, Pp, T, dur):
    x = phase_dist(t, Pp, T)
    inb = np.abs(x) < dur / 2
    ep = np.round((t - T) / Pp).astype(int)
    con = []
    for e in np.unique(ep[inb]):
        k = inb & (ep == e)
        S = -(w[k] * y[k]).sum()
        con.append((int(e), int(s[k][0]), float(S), float(w[k].sum())))
    Stot = sum(c[2] for c in con)
    pos = [c[2] for c in con if c[2] > 0]
    return dict(n_events=len(con), n_pos=len(pos), max_frac=float(max(pos) / sum(pos)) if pos else np.nan,
                events=[dict(epoch=c[0], sector=c[1], depth_ppt=c[2] / c[3] * 1e3, err_ppt=1e3 / np.sqrt(c[3])) for c in con])


def main_real():
    t, y, w, s = data()
    B = t.max() - t.min()
    fr = bls2.freq_grid(0.3, 30.0, B)
    out = dict(nfreq=len(fr), baseline=B, npts=len(t))
    for tag, sign in (("real", 1), ("inverted", -1)):
        t0 = time.time()
        r = bls2.run(t, sign * y, w, fr)
        Pg = 1 / r["freq"]
        res = []
        for i in peaks(r):
            dur = bls2.durations(Pg[i])[r["k"][i]]
            pe = per_event(t, sign * y, w, s, Pg[i], r["t0"][i], dur)
            res.append(dict(P=float(Pg[i]), T0=float(r["t0"][i]), dur_h=float(dur * 24), dchi2=float(r["dchi2"][i]),
                            depth_ppt=float(r["depth"][i] * 1e3), **pe))
        out[tag] = dict(peaks=res, bands=band_max(Pg, r["dchi2"]), runtime_s=time.time() - t0)
        print(tag, out[tag]["bands"], [(round(p["P"], 4), round(p["dchi2"], 1), p["n_events"], round(p["max_frac"], 2)) for p in res[:6]], flush=True)
        if tag == "real":
            # decimated periodogram (max in blocks of 50) for figures
            nbk = len(Pg) // 50
            np.savez_compressed(f"{ME}/data/r03_real_pgram.npz", P=Pg[:nbk * 50].reshape(nbk, 50).mean(1),
                                d2=r["dchi2"][:nbk * 50].reshape(nbk, 50).max(1))
    jdump(out, f"{ME}/out/r03_real.json")


def make_null(t, y, w, s, kind, rng):
    y2 = y.copy(); t2 = t.copy()
    if kind in ("blk", "blki"):
        for q in np.unique(s):
            idx = np.where(s == q)[0]
            kb = np.floor((t[idx] - t[idx].min()) / 1.0).astype(int)
            ub = np.unique(kb)
            perm = rng.permutation(len(ub))
            vals = np.concatenate([y[idx][kb == ub[i]] for i in perm])
            # assign shuffled values to the original time stamps in order (block lengths differ slightly: values keep order)
            y2[idx] = vals
        if kind == "blki": y2 = -y2
    elif kind == "shift":
        for q in np.unique(s):
            m = s == q
            t2[m] = t[m] + rng.uniform(0, 60.0)
        o = np.argsort(t2); t2, y2 = t2[o], y2[o]
        return t2, y2, w[o], s[o]
    return t2, y2, w, s


def main_null(kind, n, seed):
    rng = np.random.default_rng(seed)
    t, y, w, s = data()
    B = t.max() - t.min()
    res = []
    for j in range(n):
        t2, y2, w2, s2 = make_null(t, y, w, s, kind, rng)
        B2 = t2.max() - t2.min()
        fr = bls2.freq_grid(0.3, 30.0, B2)
        r = bls2.run(t2, y2, w2, fr)
        Pg = 1 / r["freq"]; i = int(np.argmax(r["dchi2"]))
        dur = bls2.durations(Pg[i])[r["k"][i]]
        pe = per_event(t2, y2, w2, s2, Pg[i], r["t0"][i], dur)
        res.append(dict(max=float(r["dchi2"][i]), P=float(Pg[i]), depth_ppt=float(r["depth"][i] * 1e3),
                        n_events=pe["n_events"], max_frac=pe["max_frac"], bands=band_max(Pg, r["dchi2"])))
        print(j, round(res[-1]["max"], 1), round(res[-1]["P"], 3), res[-1]["bands"], flush=True)
        jdump(dict(kind=kind, seed=seed, res=res), f"{ME}/out/r03_null_{kind}_{seed}.json")


if __name__ == "__main__":
    if sys.argv[1] == "real": main_real()
    else: main_null(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
