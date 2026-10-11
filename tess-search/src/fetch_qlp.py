"""Download QLP light curves for a TIC list and one sector; keep good cadences binned to 10 min.

Output: data/s{sector}/chunk_{n}.npz with concatenated arrays and per-star offsets.
"""
import sys, io, time, os
from concurrent.futures import ThreadPoolExecutor
import numpy as np, requests
from astropy.io import fits

BASE = "https://archive.stsci.edu/hlsps/qlp"
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(pool_connections=32, pool_maxsize=32)
session.mount("https://", adapter)

def url(tic, sector):
    t = f"{tic:016d}"; s = f"s{sector:04d}"
    return f"{BASE}/{s}/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/hlsp_qlp_tess_ffi_{s}-{t}_tess_v01_llc.fits"

def load(tic, sector, binsize=3):
    """Return (time, flux, meta) of good cadences binned by `binsize` cadences, or None."""
    for attempt in range(4):
        try:
            r = session.get(url(tic, sector), timeout=60)
            if r.status_code == 404: return None
            r.raise_for_status(); break
        except Exception:
            time.sleep(2 ** attempt)
    else:
        return None
    with fits.open(io.BytesIO(r.content)) as f:
        h = f[0].header; d = f[1].data
        t = np.asarray(d["TIME"], float); fl = np.asarray(d["DET_FLUX"], float)
        cad = np.asarray(d["CADENCENO"]); q = np.asarray(d["QUALITY"])
        meta = (h.get("CAMERA", 0), h.get("CCD", 0))
    ok = (q == 0) & np.isfinite(t) & np.isfinite(fl) & (fl > 0)
    if ok.sum() < 300: return None
    t, fl, cad = t[ok], fl[ok], cad[ok]
    key = cad // binsize
    u, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    tb = np.bincount(inv, t) / cnt; fb = np.bincount(inv, fl) / cnt
    good = cnt >= 2
    return tb[good], fb[good] / np.median(fb[good]), cnt[good], meta

def main(listfile, sector, outdir, chunk=5000, workers=24):
    tics = np.loadtxt(listfile, dtype=np.int64, ndmin=1)
    os.makedirs(outdir, exist_ok=True)
    for c0 in range(0, len(tics), chunk):
        out = f"{outdir}/chunk_{c0 // chunk:04d}.npz"
        if os.path.exists(out): continue
        part = tics[c0:c0 + chunk]; t0 = time.time()
        with ThreadPoolExecutor(workers) as ex:
            res = list(ex.map(lambda x: load(int(x), sector), part))
        ids, offs, T, F, N, cams = [], [0], [], [], [], []
        for tic, r in zip(part, res):
            if r is None: continue
            ids.append(tic); T.append(r[0]); F.append(r[1].astype(np.float32)); N.append(r[2].astype(np.int8))
            cams.append(r[3]); offs.append(offs[-1] + len(r[0]))
        np.savez(out, tic=np.array(ids), off=np.array(offs), time=np.concatenate(T), flux=np.concatenate(F),
                 n=np.concatenate(N), camccd=np.array(cams, dtype=np.int8))
        print(f"{out}: {len(ids)}/{len(part)} stars, {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
