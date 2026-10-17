"""Referee r05: local-window pixel photometry (independent of the investigator's global regression + Mira template).

For every window centred on a (real or fake) transit time tc (half-width HALF = 0.35 d; 0.6 d for 30-min S12):
 1. cadences: finite, QUALITY in {0, 2048}, background B < median + BMAX (strict 60 / relaxed 800 e/s/px),
    pointing within 6 robust sigma; both sides of the window must keep >= 25% of the out-of-transit cadences.
 2. every pixel within 3.5 px of the M3V or the Mira: baseline = [1, tau, tau^2, B, dx, dy] (+ dx^2, dy^2, dxdy if
    ORDER2) fitted to the out-of-transit cadences only (|t - tc| > MASK = 0.07 d), residuals r_p(t) for all cadences.
    (The Mira's slow variation, the background and the pointing response are absorbed per pixel and per window.)
 3. per cadence weighted fits of the residual frame (weights 1/sigma_p^2 from out-of-transit residuals):
    'M': a_M3V * PRF(M3V) + const over pixels within 2.5 px of the M3V  -> M3V light curve (flux change, e-/s)
    'T': a_M3V * PRF(M3V) + a_pair * PRF(pair A+B) + a_Mira * b_Mira + const within 3.5 px -> all three stars
 4. stores the mean residual frame in the flat part (|t - tc| < 0.02 d) and sigma_p for difference-image tests.
Pixel-level injection (inject=[(img, depth_e_s, tc, T14, flat)]) subtracts depth*shape(t)*img from the cube before
step 2, so the full method can be tested with templates that differ from the fitting ones.
Positions: Gaia DR3 propagated + r01 field-star offset; pair = A + 0.367 B (Gaia G ratio).
Functions are imported by later scripts; running this file processes the real transits at the investigator's
global ephemeris and 120 fake epochs per sector and writes data/win_s{sector}.npz."""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *

P_G, T0_G = 5.6915087, 3840.91414          # investigator's global ephemeris (used only to place windows)
MASK = 0.07


def half_for(s): return 0.6 if s == 12 else (0.4 if s == 39 else 0.35)


def trap_shape(t, tc, T14, fl, exptime, nsub=9):
    y = np.zeros(len(t))
    for o in (np.arange(nsub) + 0.5) / nsub - 0.5:
        a = np.abs(t + o * exptime - tc); T23 = fl * T14
        y += np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-9)))
    return y / nsub


class Sector:
    def __init__(self, s):
        self.s = s; C = self.C = Cut(s)
        g = C.stars(); W = json.load(open(f"{R}/out/wcs_s{s}.json")); ox, oy = W["offset_px"]
        g["x"] += ox; g["y"] += oy; self.g = g
        self.a = np.load(f"{R}/data/aux_s{s}.npz")
        z = np.load(f"{R}/data/tcorr_s{s}.npy"); self.ttdb = C.t + np.interp(C.t, z[0], z[1]) / 86400
        pos = {k: g[g.source_id == GID[k]].iloc[0] for k in ("Mira", "M3V", "pA", "pB", "bg3")}
        self.pos = {k: (float(v.x), float(v.y)) for k, v in pos.items()}
        yy, xx = np.mgrid[0:C.ny, 0:C.nx]; self.xx, self.yy = xx, yy
        mx, my = self.pos["Mira"]; m3x, m3y = self.pos["M3V"]
        self.reg = ((np.hypot(xx - mx, yy - my) <= 3.5) | (np.hypot(xx - m3x, yy - m3y) <= 3.5)).ravel()
        self.regM = (np.hypot(xx - m3x, yy - m3y) <= 2.5).ravel()
        self.P_M3V = C.psf(m3x, m3y)
        self.P_pair = (C.psf(*self.pos["pA"]) + 0.367 * C.psf(*self.pos["pB"])) / 1.367
        b = np.load(f"{R}/data/mirab_s{s}.npz")["b"]
        bm = np.where(np.hypot(xx - mx, yy - my) <= 4.0, b, 0.0); self.b_Mira = bm / bm.sum()
        self.P_Mira = C.psf(mx, my)
        q = C.q; B = self.a["B"]
        self.qok = C.fin & ((q == 0) | (q == 2048)) & np.isfinite(self.a["dx"])
        q0 = C.fin & (q == 0)
        self.Bmed = float(np.nanmedian(B[q0]))
        jit = np.ones(C.n, bool)
        for k in ("dx", "dy"):
            v = self.a[k]; jit &= np.abs(v - np.nanmedian(v[q0])) < 6 * rstd(v[q0]) + 0.02
        self.qok &= jit
        self.flat = C.f.reshape(C.n, -1)
        self.exptime = C.dt

    def window(self, tc, bmax=60.0, half=None, order2=False, inject=None, tmpl_M=None, keep_cube=False, linbase=False):
        C = self.C; half = half or half_for(self.s)
        t = self.ttdb
        w = self.qok & (np.abs(t - tc) < half) & (self.a["B"] < self.Bmed + bmax)
        oot = w & (np.abs(t - tc) > MASK)
        nl = np.sum(oot & (t < tc)); nr = np.sum(oot & (t > tc))
        nexp = 2 * (half - MASK) / self.exptime
        if nl < 0.25 * nexp / 2 or nr < 0.25 * nexp / 2 or np.sum(w & ~oot) < 2: return None
        Y = self.flat[w][:, self.reg].copy()
        if inject:
            for img, dep, tci, T14, fl in inject:
                Y -= dep * trap_shape(t[w], tci, T14, fl, self.exptime)[:, None] * img.ravel()[self.reg][None]
        tau = t[w] - tc; Bw = self.a["B"][w]; dx = self.a["dx"][w]; dy = self.a["dy"][w]
        o = oot[w]
        cols = [np.ones(w.sum()), tau] + ([] if linbase else [tau ** 2]) + [Bw - Bw[o].mean(), dx - dx[o].mean(), dy - dy[o].mean()]
        if order2:
            ddx, ddy = dx - dx[o].mean(), dy - dy[o].mean(); cols += [ddx ** 2, ddy ** 2, ddx * ddy]
        X = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(X[o], Y[o], rcond=None)
        Rz = Y - X @ coef
        dof = max(1, o.sum() - X.shape[1])
        sig = np.sqrt((Rz[o] ** 2).sum(0) / dof)
        out = dict(t=t[w], tc=tc, n_oot=int(o.sum()), oot=o, sig=sig, B=Bw)
        # 'M' fit
        TM = (tmpl_M if tmpl_M is not None else self.P_M3V).ravel()[self.reg]
        mM = self.regM[self.reg]
        A = np.column_stack([TM[mM], np.ones(mM.sum())]) / sig[mM, None]
        cov = np.linalg.pinv(A.T @ A); cM = cov @ A.T @ (Rz[:, mM] / sig[mM]).T
        out["aM"] = cM[0]; out["aM_ferr"] = float(np.sqrt(cov[0, 0]))
        # 'T' fit (three stars)
        A3 = np.column_stack([TM, self.P_pair.ravel()[self.reg], self.b_Mira.ravel()[self.reg], np.ones(self.reg.sum())]) / sig[:, None]
        cov3 = np.linalg.pinv(A3.T @ A3); c3 = cov3 @ A3.T @ (Rz / sig).T
        out["aT"] = c3[:3].T; out["aT_ferr"] = np.sqrt(np.diag(cov3)[:3])
        # mean residual frame in the flat part, and in the whole transit
        fin_ = np.abs(tau) < 0.02
        img = np.full(self.reg.shape, np.nan); img[self.reg] = Rz[fin_].mean(0) if fin_.sum() else np.nan
        out["dimg"] = img; out["n_flat"] = int(fin_.sum())
        sg = np.full(self.reg.shape, np.nan); sg[self.reg] = sig; out["sigimg"] = sg
        if keep_cube: out["Rz"] = Rz
        return out


def windows_for(S, centres, **kw):
    res = []
    for tc in centres:
        r = S.window(tc, **kw)
        if r is not None: res.append(r)
    return res


def real_centres(S, P=P_G, T0=T0_G):
    t = S.ttdb[np.isfinite(S.ttdb)]
    k = np.arange(np.ceil((t.min() - T0) / P - 0.1), np.floor((t.max() - T0) / P + 0.1) + 1)
    return T0 + k * P, k.astype(int)


def fake_centres(S, n, rng, P=P_G, T0=T0_G, avoid=0.45):
    t = S.ttdb[S.qok]
    out = []
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
        ph = ((tc - T0) / P + 0.5) % 1 - 0.5
        if abs(ph * P) < avoid: continue
        if out and np.min(np.abs(np.array(out) - tc)) < 0.15: continue
        out.append(tc)
    return np.array(out)


def pack(res):
    """Flatten window results into arrays for saving."""
    if not res: return {}
    d = dict(tc=np.array([r["tc"] for r in res]), n_oot=np.array([r["n_oot"] for r in res]), n_flat=np.array([r["n_flat"] for r in res]),
             dimg=np.array([r["dimg"] for r in res]), sigimg=np.array([r["sigimg"] for r in res]),
             aM_ferr=np.array([r["aM_ferr"] for r in res]), aT_ferr=np.array([r["aT_ferr"] for r in res]))
    idx = np.concatenate([[i] * len(r["t"]) for i, r in enumerate(res)])
    d.update(win=idx, t=np.concatenate([r["t"] for r in res]), aM=np.concatenate([r["aM"] for r in res]),
             aT=np.concatenate([r["aT"] for r in res]), oot=np.concatenate([r["oot"] for r in res]), B=np.concatenate([r["B"] for r in res]))
    return d


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        S = Sector(s)
        cen, ks = real_centres(S)
        out = {}
        for lab, bmax in (("strict", 60.0), ("relaxed", 800.0)):
            res = []; kk = []
            for tc, k in zip(cen, ks):
                r = S.window(tc, bmax=bmax)
                if r is not None: res.append(r); kk.append(k)
            for key, v in pack(res).items(): out[f"{lab}_{key}"] = v
            out[f"{lab}_epoch"] = np.array(kk)
        fk = fake_centres(S, 150, rng)
        resf = windows_for(S, fk, bmax=60.0)
        for key, v in pack(resf).items(): out[f"fake_{key}"] = v
        np.savez_compressed(f"{R}/data/win_s{s}.npz", **out, exptime=S.exptime, pos=json.dumps(S.pos))
        print(f"S{s}: real strict epochs {out['strict_epoch'].tolist()}, relaxed {out['relaxed_epoch'].tolist()}; fake windows {len(resf)}; "
              f"median per-cadence aM scatter (fake, OOT) {rstd(out['fake_aM'][out['fake_oot']]):.1f} e/s", flush=True)
