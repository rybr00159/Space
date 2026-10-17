"""c03: transit-window data sets from the cleaned light curves of both pipelines, and the noise calibration at fake
epochs.

INV (mode 'all'): quality_ok cadences, flare-flagged points (> 4 robust sigma, with neighbours) removed, windows of
+-0.30 d around each predicted transit (phase-1 ephemeris). Two flux versions are compared:
   'det'  = wotan-detrended flux (phase-1 p06 recipe), baseline = per-window offset + slope;
   'raw'  = undetrended normalised flux 1 + dF/979, baseline = per-window offset + slope (or quadratic, 'raw2').
REF (mode 'all'): the referee's strict windows (+-0.35 d; S39 +-0.4, S12 +-0.6 d), f = 1 + aM/979, baseline offset
+ slope (the per-pixel quadratic baselines are already removed in the pipeline).
Window acceptance (both): >= 2 in-transit points, >= 25% of the expected out-of-transit cadences on each side.
Per-point sigma = robust std of the window's out-of-transit points.

Noise calibration: the same window construction at FAKE epochs (phase-1 ephemeris shifted by random phases, at least
0.45 d from any real transit; 150 per sector for INV; the referee pipeline's own 150 fake windows for REF); a fixed
template (trapezoid T14 1.68 h, flat 0.72, i.e. the shape/referee2_shape solution) depth is fitted in each, and
beta_s = std(fake single-window depths) / median(formal single-window error). Fake-window depths also give the
empirical single-transit error.
Outputs: data/win_{inv,ref}_{version}.npz (+ fake), out/noise_calib.json
"""
import sys, json
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import CL, P1, T01, F_M3V, SECTORS, robust_std, jd, trap_frac
from tfit import *

HALF_INV = 0.30
T14_T, FL_T = 1.68 / 24, 0.72
rng = np.random.default_rng(77)


def accept(tt, tc, half, exptime):
    oot = np.abs(tt - tc) > 0.07
    nexp = 2 * (half - 0.07) / exptime
    nl = np.sum(oot & (tt < tc)); nr = np.sum(oot & (tt > tc))
    nin = np.sum(np.abs(tt - tc) < T14_T / 2 * 0.8)
    return nl >= 0.25 * nexp / 2 and nr >= 0.25 * nexp / 2 and nin >= 2


def inv_sector(s, mode="all"):
    z = np.load(f"{CL}/data/lcinv_{mode}_s{s}.npz")
    t = z["t_tdb"]; ok = z["quality_ok"].astype(bool) & np.isfinite(z["det_M3V_B"]) & np.isfinite(t)
    det = z["det_M3V_B"]; raw = z["raw_M3V_B"]
    sd = robust_std(det[ok & ~z["intransit_window"].astype(bool)])
    hi = ok & (det > 1 + 4 * sd); fl = hi | np.roll(hi, 1) | np.roll(hi, -1)
    ok &= ~fl
    return dict(t=t[ok], det=det[ok], raw=raw[ok], ex=float(z["exptime"]), nflare=int(hi.sum()))


def build_inv(centres_by_sector, version, half=HALF_INV):
    T, F, S, Wd, TC, SEC, EP, EX = [], [], [], [], [], [], [], []
    i = 0
    for s, cen in centres_by_sector.items():
        d = INV[s]
        for tc, ep in cen:
            m = np.abs(d["t"] - tc) < half
            if m.sum() < 8 or not accept(d["t"][m], tc, half, d["ex"]): continue
            f = d["det" if version == "det" else "raw"][m]
            oot = np.abs(d["t"][m] - tc) > 0.07
            sg = robust_std(f[oot] - np.polyval(np.polyfit(d["t"][m][oot] - tc, f[oot], 1), d["t"][m][oot] - tc))
            T.append(d["t"][m]); F.append(f); S.append(np.full(m.sum(), sg)); Wd.append(np.full(m.sum(), i))
            TC.append(tc); SEC.append(s); EP.append(ep); EX.append(np.full(m.sum(), d["ex"])); i += 1
    return Windows(np.concatenate(T), np.concatenate(F), np.concatenate(S), np.concatenate(Wd), TC, SEC, EP,
                   np.concatenate(EX), order=2 if version == "raw2" else 1)


def build_ref(mode="all", kind="strict"):
    T, F, S, Wd, TC, SEC, EP, EX = [], [], [], [], [], [], [], []
    i = 0
    for s in SECTORS:
        z = np.load(f"{CL}/data/refwin_{mode}_s{s}.npz")
        if f"{kind}_t" not in z.files: continue
        ex = float(z["exptime"])
        eps = z[f"{kind}_epoch"] if kind != "fake" else np.full(len(z["fake_tc"]), -9999)
        for j, tc in enumerate(z[f"{kind}_tc"]):
            m = z[f"{kind}_win"] == j
            tt = z[f"{kind}_t"][m]; f = 1 + z[f"{kind}_aM"][m] / F_M3V; oot = z[f"{kind}_oot"][m]
            if not accept(tt, tc, {12: 0.6, 39: 0.4}.get(s, 0.35), ex): continue
            sg = robust_std(f[oot] - np.polyval(np.polyfit(tt[oot] - tc, f[oot], 1), tt[oot] - tc))
            T.append(tt); F.append(f); S.append(np.full(m.sum(), sg)); Wd.append(np.full(m.sum(), i))
            TC.append(tc); SEC.append(s); EP.append(eps[j]); EX.append(np.full(m.sum(), ex)); i += 1
    return Windows(np.concatenate(T), np.concatenate(F), np.concatenate(S), np.concatenate(Wd), TC, SEC, EP,
                   np.concatenate(EX), order=1)


def real_centres(s):
    t = INV[s]["t"]
    k = np.arange(np.ceil((t.min() - T01) / P1 - 0.1), np.floor((t.max() - T01) / P1 + 0.1) + 1).astype(int)
    return [(T01 + kk * P1, kk) for kk in k]


def fake_centres(s, n=150):
    t = INV[s]["t"]; out = []
    tries = 0
    while len(out) < n and tries < 100 * n:
        tries += 1
        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
        ph = ((tc - T01) / P1 + 0.5) % 1 - 0.5
        if abs(ph * P1) < 0.45: continue
        if out and np.min(np.abs(np.array([o[0] for o in out]) - tc)) < 0.15: continue
        out.append((tc, -9999))
    return out


def single_depths(W):
    dep = []; err = []
    sh = np.zeros(len(W.t))
    for i, ii in enumerate(W.idx):
        sh[ii] = trap_frac(W.t[ii], W.tc0[i], T14_T, FL_T, float(W.ex[ii][0]))
        X = np.column_stack([-sh[ii]] + [W.tau[ii] ** o for o in range(W.order + 1)]) / W.s[ii, None]
        c, *_ = np.linalg.lstsq(X, W.f[ii] / W.s[ii], rcond=None)
        cov = np.linalg.pinv(X.T @ X); dep.append(c[0]); err.append(np.sqrt(cov[0, 0]))
    return np.array(dep), np.array(err)


def save_windows(W, fn):
    np.savez_compressed(fn, t=W.t, f=W.f, s=W.s, w=W.w, tc0=W.tc0, sector=W.sector, epoch=W.epoch, ex=W.ex, order=W.order)


def load_windows(fn):
    z = np.load(fn)
    return Windows(z["t"], z["f"], z["s"], z["w"], z["tc0"], z["sector"], z["epoch"], z["ex"], int(z["order"]))


if __name__ == "__main__":
    INV = {s: inv_sector(s) for s in SECTORS}
    real = {s: real_centres(s) for s in SECTORS}
    fake = {s: fake_centres(s) for s in SECTORS}
    out = {"nflare_removed": {str(s): INV[s]["nflare"] for s in SECTORS}}
    for ver in ("det", "raw", "raw2"):
        Wr = build_inv(real, ver); Wf = build_inv(fake, ver)
        save_windows(Wr, f"{CL}/data/win_inv_{ver}.npz"); save_windows(Wf, f"{CL}/data/win_inv_{ver}_fake.npz")
        dr, er = single_depths(Wr); df, ef = single_depths(Wf)
        cal = {}
        for s in SECTORS:
            m = Wf.sector == s; mr = Wr.sector == s
            cal[str(s)] = dict(n_fake=int(m.sum()), fake_mean_ppt=float(np.mean(df[m]) * 1e3), fake_std_ppt=float(np.std(df[m]) * 1e3),
                               fake_rstd_ppt=float(robust_std(df[m]) * 1e3), formal_median_ppt=float(np.median(ef[m]) * 1e3),
                               beta=float(np.std(df[m]) / np.median(ef[m])), n_real=int(mr.sum()),
                               real_epochs=Wr.epoch[mr].tolist(), real_depths_ppt=np.round(dr[mr] * 1e3, 2).tolist())
        out[f"inv_{ver}"] = cal
        print(ver, {s: (c["n_real"], round(c["fake_std_ppt"], 2), round(c["beta"], 2), round(c["fake_mean_ppt"], 2)) for s, c in cal.items()}, flush=True)
    Wr = build_ref("all", "strict"); Wf = build_ref("all", "fake")
    save_windows(Wr, f"{CL}/data/win_ref.npz"); save_windows(Wf, f"{CL}/data/win_ref_fake.npz")
    Wo = build_ref("orig", "strict"); save_windows(Wo, f"{CL}/data/win_reforig.npz")
    Wof = build_ref("orig", "fake"); save_windows(Wof, f"{CL}/data/win_reforig_fake.npz")
    for lab, (A, B) in (("ref", (Wr, Wf)), ("reforig", (Wo, Wof))):
        dr, er = single_depths(A); df, ef = single_depths(B)
        cal = {}
        for s in SECTORS:
            m = B.sector == s; mr = A.sector == s
            cal[str(s)] = dict(n_fake=int(m.sum()), fake_mean_ppt=float(np.mean(df[m]) * 1e3), fake_std_ppt=float(np.std(df[m]) * 1e3),
                               fake_rstd_ppt=float(robust_std(df[m]) * 1e3), formal_median_ppt=float(np.median(ef[m]) * 1e3),
                               beta=float(np.std(df[m]) / np.median(ef[m])), n_real=int(mr.sum()),
                               real_epochs=A.epoch[mr].tolist(), real_depths_ppt=np.round(dr[mr] * 1e3, 2).tolist())
        out[lab] = cal
        print(lab, {s: (c["n_real"], round(c["fake_std_ppt"], 2), round(c["beta"], 2), round(c["fake_mean_ppt"], 2)) for s, c in cal.items()}, flush=True)
    jd(out, f"{CL}/out/noise_calib.json")
