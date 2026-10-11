"""Stage 2: test whether a single-sector signal recurs in every other TESS sector of the same star.

For each candidate: fetch all sectors (SPOC 2-min > TESS-SPOC FFI > QLP), detrend, then
  * per sector: narrow BLS around the candidate period, compared with the same search at random periods
    (empirical false-alarm probability);
  * all sectors jointly: fine BLS around the period to get a multi-year ephemeris and a joint S/N.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("NUMBA_NUM_THREADS", "1")
import sys, io, time, json, warnings
import numpy as np, pandas as pd, requests
from astropy.io import fits

warnings.filterwarnings("ignore")
from bls import _bls, detrend, clean, diagnostics, robust_std

DL = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
CACHE = os.path.join(os.path.dirname(__file__), "..", "lc_cache")
session = requests.Session()


def products(tic):
    from astroquery.mast import Observations
    obs = Observations.query_criteria(target_name=[str(tic)], provenance_name=["QLP", "SPOC", "TESS-SPOC"])
    rows = []
    for r in obs:
        url = str(r["dataURL"])
        if not url.endswith(("_lc.fits", "_llc.fits")) or "fast-lc" in url: continue
        rows.append(dict(sector=int(r["sequence_number"]), prov=str(r["provenance_name"]), exptime=float(r["t_exptime"]), uri=url))
    return pd.DataFrame(rows)


def fetch(uri):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, uri.split("/")[-1])
    if not os.path.exists(path):
        for a in range(4):
            try:
                r = session.get(DL + uri, timeout=120); r.raise_for_status()
                open(path, "wb").write(r.content); break
            except Exception:
                time.sleep(2 ** a)
        else:
            return None
    return path


def read_lc(path, prov):
    with fits.open(path) as f:
        d = f[1].data
        t = np.asarray(d["TIME"], float)
        if prov == "QLP":
            col = next(c for c in ("DET_FLUX", "KSPSAP_FLUX", "SAP_FLUX") if c in d.columns.names)
            fl = np.asarray(d[col], float); q = np.asarray(d["QUALITY"])
        else:
            fl = np.asarray(d["PDCSAP_FLUX"], float); q = np.asarray(d["QUALITY"])
    ok = (q == 0) & np.isfinite(t) & np.isfinite(fl) & (fl > 0)
    t, fl = t[ok], fl[ok]
    if len(t) < 100: return None
    # bin to 10 min (no-op for 30-min data)
    dt = np.median(np.diff(t))
    if dt < 500 / 86400:
        k = np.floor((t - t[0]) / (600 / 86400)).astype(int)
        _, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
        t = np.bincount(inv, t) / cnt; fl = np.bincount(inv, fl) / cnt
    return t, fl / np.median(fl)


def load_all(tic, prefer=("SPOC", "TESS-SPOC", "QLP")):
    p = products(tic)
    if p.empty: return []
    out = []
    for s, g in p.groupby("sector"):
        for prov in prefer:
            h = g[g.prov == prov]
            if len(h):
                r = h.sort_values("exptime").iloc[0]
                path = fetch(r.uri)
                lc = read_lc(path, prov) if path else None
                if lc is not None:
                    t, f = detrend(*lc)
                    t, f, s_ = clean(t, f)
                    out.append(dict(sector=int(s), prov=prov, exptime=r.exptime, t=t, f=f, sigma=s_))
                    break
    return out


def narrow_snr(t, f, sigma, P, dur, rel=0.004, nper=41):
    """Best box S/N within P*(1 +- rel) and durations around `dur`."""
    y = f - 1.0; w = np.full(len(t), 1.0 / sigma ** 2); y = y - np.sum(w * y) / np.sum(w)
    fr = 1.0 / np.linspace(P * (1 - rel), P * (1 + rel), nper)
    durs = np.array([dur / 1.5, dur, dur * 1.5])
    pw, dp, t0, du = _bls(t - t.min(), y, w, fr, durs, min(dur / 6, 10 / 1440))
    i = int(np.argmax(pw))
    return np.sqrt(pw[i]), 1.0 / fr[i], t0[i] + t.min(), du[i], dp[i]


def sector_test(sec, P, dur, nnull=150, rng=None):
    rng = rng or np.random.default_rng(0)
    t, f, s = sec["t"], sec["f"], sec["sigma"]
    if t.max() - t.min() < 1.2 * P:
        return dict(snr=np.nan, fap=np.nan)
    snr, Pb, t0, du, dp = narrow_snr(t, f, s, P, dur)
    null = []
    while len(null) < nnull:
        Pr = 10 ** rng.uniform(np.log10(max(0.5, P / 3)), np.log10(min(P * 3, (t.max() - t.min()) / 1.2)))
        if min(abs(Pr / P - k) for k in (0.5, 1, 2, 1 / 3, 3)) < 0.02: continue
        null.append(narrow_snr(t, f, s, Pr, dur)[0])
    null = np.array(null)
    return dict(snr=snr, P=Pb, t0=t0, dur=du, depth=dp, fap=(np.sum(null >= snr) + 1) / (len(null) + 1),
                null_med=float(np.median(null)), null_max=float(null.max()))


def joint_fit(secs, P, dur, rel=0.003):
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    w = np.concatenate([np.full(len(s["t"]), 1 / s["sigma"] ** 2) for s in secs])
    base = t.max() - t.min()
    dP = P * (dur / 6) / base  # phase drift of dur/6 across the baseline
    n = int(min(20000, 2 * rel * P / dP) + 1)
    fr = 1.0 / np.linspace(P * (1 - rel), P * (1 + rel), n)
    y = f - 1.0; y = y - np.sum(w * y) / np.sum(w)
    durs = np.array([dur / 1.5, dur, dur * 1.5])
    pw, dp, t0, du = _bls(t - t.min(), y, w, fr, durs, min(dur / 6, 10 / 1440))
    i = int(np.argmax(pw)); sp = np.sqrt(pw)
    # alias ambiguity: second-best peak more than 2 period-steps away
    far = np.abs(np.arange(n) - i) > 3
    second = sp[far].max() if far.any() else 0.0
    return dict(snr=sp[i], P=1.0 / fr[i], t0=t0[i] + t.min(), dur=du[i], depth=dp[i], second=second, n=n,
                baseline=base)


def run(cand):
    try:
        return _run(cand)
    except Exception as e:
        return dict(tic=int(cand["tic"]), error=repr(e)[:300])


def _run(cand):
    tic, P, dur = int(cand["tic"]), float(cand["P"]), float(cand["dur"])
    disc = {int(x) for x in str(cand["sector"]).split(",")}  # discovery sector(s): never counted as confirmation
    s0 = min(disc)
    secs = load_all(tic)
    if not secs: return dict(tic=tic, error="no data")
    rng = np.random.default_rng(tic % 2 ** 32)
    per = []
    for s in secs:
        r = sector_test(s, P, dur, rng=rng); r.update(sector=s["sector"], prov=s["prov"], sigma=s["sigma"])
        per.append(r)
    others = [r for r in per if r["sector"] not in disc and np.isfinite(r["snr"])]
    j = joint_fit(secs, P, dur)
    # per-sector depth at the joint ephemeris
    for s, r in zip(secs, per):
        d = diagnostics(s["t"], s["f"], j["P"], j["t0"], j["dur"], s["sigma"])
        r.update(depth_joint=d["depth"], depth_err=d["depth_err"], ntr_joint=d["n_tr"])
    tc = [s for s in secs if s["sector"] not in disc]
    jo = joint_fit(tc, P, dur) if tc else None  # ephemeris without the discovery sector
    from scipy.stats import chi2
    fisher = float(-2 * np.sum(np.log([r["fap"] for r in others]))) if others else np.nan
    # depth consistency: weighted mean of other-sector depths at the joint ephemeris vs the discovery sector
    od = [(r["depth_joint"], r["depth_err"]) for r in others if r["depth_err"] > 0]
    if od:
        wts = np.array([1 / e ** 2 for _, e in od]); dm = np.sum([d * w for (d, _), w in zip(od, wts)]) / wts.sum()
        de = 1 / np.sqrt(wts.sum())
    else:
        dm = de = np.nan
    disc_rec = next((r for r in per if r["sector"] == s0), None)
    return dict(tic=tic, n_sectors=len(secs), sectors=[r["sector"] for r in per],
                n_other=len(others), n_other_sig=int(sum(r["fap"] <= 0.01 for r in others)),
                fisher=fisher, p_comb=float(chi2.sf(fisher, 2 * len(others))) if others else np.nan,
                other_depth=dm, other_depth_err=de, other_depth_snr=dm / de if de == de and de > 0 else np.nan,
                disc_depth=disc_rec["depth_joint"] if disc_rec else np.nan,
                joint=j, others_only=jo, per=per)


if __name__ == "__main__":
    cands = pd.read_csv(sys.argv[1]); out = sys.argv[2]
    from multiprocessing import get_context
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["tic"] for l in open(out)}
    todo = [c for c in cands.to_dict("records") if int(c["tic"]) not in done]
    with get_context("spawn").Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 3) as pool, open(out, "a") as fo:
        for r in pool.imap_unordered(run, todo):
            fo.write(json.dumps(r, default=float) + "\n"); fo.flush()
            print(r["tic"], r.get("n_sectors"), r.get("n_other_sig"), round(r.get("joint", {}).get("snr", 0), 1) if "joint" in r else r.get("error"), flush=True)
