"""cleanlib: rebuild of the phase-1 'INV' (phot model B) photometry with every catalogued variable star excluded from
the field-wide estimators (pointing dx/dy and flux scale, faint-pixel background plane, faint-pixel background level),
optionally with the variables' own light curves as explicit regressors.

The code is a modified copy of ../phot/scripts/p04_jitter.py, p05_phot.py and p06_lc.py (read-only imports of
photlib.py / p03_scene.py from ../phot/scripts; nothing is written into ../phot). Changes are marked 'CLEAN:'.

Variable masks (mode):
  'none' : phase-1 behaviour (no masking; used to verify that this copy reproduces phot)
  'rr3'  : the Gaia DR3 RRab 5925372615693467904 and the two Gaia DR3 EBs 5925384852082638720 (P 0.7503 d) and
           5925378392451781120 (P 1.5953 d) masked
  'all'  : every Gaia DR3 phot_variable_flag == VARIABLE source with G < 18.6 inside the cutout (except the Mira,
           which phot already treats separately): radius 2.5 px (G >= 14), 3.0 px (12 <= G < 14), 3.5 px (G < 12)
"""
import os, sys, json
for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[k] = "1"
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *            # ROOT = ../phot (read-only inputs)
from p03_scene import build_stars
from wotan import flatten

BASE = "tic214311908/analysis"
CL = f"{BASE}/clean"
PH = f"{BASE}/phot"
RF = f"{BASE}/referee_phot"
VARPOS = json.load(open(f"{CL}/out/varpos.json"))
RR_ID, EB075_ID, EB159_ID = 5925372615693467904, 5925384852082638720, 5925378392451781120
P_RR = 0.6485275345290323; TMAX_RR = 1707.3127882019542 + 2455197.5 - 2457000.0     # Gaia DR3 vari_rrlyrae
P_EB075, P_EB159 = 1 / 1.332751955493649, 1 / 0.6268244050997452
# phase-1 ephemeris (phot/out/timing.json global fit) used to place windows
P1, T01 = 5.6915087, 3840.9141
F_M3V = 979.0                     # M3V flux in PRF units (referee_phot r09/r17)
FREF_INV = dict(M3V=F_M3V, pair=tmag_flux(13.3749) + tmag_flux(14.56))


def var_mask(sector, mode, ny=21, nx=21, radius_scale=1.0):
    """Boolean image, True = pixel excluded from field-wide estimators."""
    m = np.zeros((ny, nx), bool)
    if mode == "none": return m
    yy, xx = np.mgrid[0:ny, 0:nx]
    for v in VARPOS[str(sector)]["vars"]:
        if v["source_id"] == IDS["Mira"]: continue
        if mode == "rr3" and v["source_id"] not in (RR_ID, EB075_ID, EB159_ID): continue
        if mode == "all" and v["G"] >= 18.6: continue
        r = (2.5 if v["G"] >= 14 else (3.0 if v["G"] >= 12 else 3.5)) * radius_scale
        m |= np.hypot(xx - v["x"], yy - v["y"]) <= r
    return m


# ------------------------------------------------------------------------------------------------ pointing/background
def static_model(C, W, prf, shift=(0.0, 0.0)):
    g = C.gaia(18.5)
    g = g[(g.x > -4) & (g.x < C.nx + 3) & (g.y > -4) & (g.y < C.ny + 3)].reset_index(drop=True)
    dx, dy = W["global_shift"]; ff = W["free_fluxes"]
    M = np.zeros((C.ny, C.nx))
    for s_, T, x, y in zip(g.source_id, g["T"], g.x, g.y):
        if s_ == IDS["Mira"]: continue
        f_ = ff.get(str(int(s_)), tmag_flux(T) * W["faint_scale"])
        M += f_ * prf.image(x + dx + shift[0], y + dy + shift[1])
    return M, g


def jitter_clean(sector, mode, C=None, save=True):
    """p04_jitter with the variable pixels excluded (CLEAN: mask &= ~var_mask; faint set excludes var pixels)."""
    C = C or Cutout(sector); W = json.load(open(f"{PH}/out/wcs_s{sector}.json"))
    prf = PRFModel(C); prf.blur = W["blur"]
    M0, g = static_model(C, W, prf)
    h = 0.05
    Gx = (static_model(C, W, prf, (h, 0))[0] - static_model(C, W, prf, (-h, 0))[0]) / (2 * h)
    Gy = (static_model(C, W, prf, (0, h))[0] - static_model(C, W, prf, (0, -h))[0]) / (2 * h)
    mira = g[g.source_id == IDS["Mira"]].iloc[0]; dx0, dy0 = W["global_shift"]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    vm = var_mask(sector, mode, C.ny, C.nx)
    mask = (np.hypot(xx - mira.x - dx0, yy - mira.y - dy0) > 3.0) & ~vm          # CLEAN
    u = (xx[mask] - 10) / 10; v = (yy[mask] - 10) / 10
    A = np.column_stack([M0[mask], Gx[mask], Gy[mask], np.ones(mask.sum()), u, v, u * u, v * v, u * v])
    cube = C.flux.reshape(len(C.t), -1)[:, mask.ravel()]
    fin = C.fin
    good = fin & (C.q == 0)
    sig = 1.4826 * np.median(np.abs(np.diff(cube[good], axis=0)), axis=0) / np.sqrt(2) + 1e-3
    Aw = A / sig[:, None]
    out = np.full((len(C.t), A.shape[1]), np.nan); chi2r = np.full(len(C.t), np.nan)
    Y = (cube[fin] / sig).T
    coef, *_ = np.linalg.lstsq(Aw, Y, rcond=None)
    out[fin] = coef.T
    resid = Aw @ coef - Y
    chi2r[fin] = np.sum(resid ** 2, axis=0) / (mask.sum() - A.shape[1])
    s_ = out[:, 0]
    dx = out[:, 1] / s_; dy = out[:, 2] / s_
    allpix = C.flux.reshape(len(C.t), -1)
    med = np.nanmedian(allpix[good], 0)
    ok = ~vm.ravel()                                                             # CLEAN
    faint = ok & (med < np.percentile(med[ok], 15))
    bkg = np.nanmedian(allpix[:, faint], 1)
    J = dict(t=C.t, q=C.q, dx=dx, dy=dy, scale=s_, bkg=bkg, chi2r=chi2r, fin=fin, varmask=vm, npix_fit=int(mask.sum()))
    if save: np.savez(f"{CL}/data/jitter_{mode}_s{sector}.npz", **J)
    return J


# ------------------------------------------------------------------------------------------------ photometry
def in_windows(t, epochs_t, half):
    if len(epochs_t) == 0: return np.zeros(len(t), bool)
    e = np.sort(np.asarray(epochs_t))
    i = np.clip(np.searchsorted(e, t), 1, len(e) - 1)
    return np.minimum(np.abs(t - e[i - 1]), np.abs(t - e[i])) < half


def transit_times(t, P, T0):
    k = np.arange(np.floor((np.nanmin(t) - T0) / P) - 1, np.ceil((np.nanmax(t) - T0) / P) + 2)
    return T0 + k * P


def binned(t, Y, width):
    k = np.floor((t - t[0]) / width).astype(int)
    u, inv = np.unique(k, return_inverse=True)
    n = np.bincount(inv)
    out = np.zeros((len(u),) + Y.shape[1:])
    np.add.at(out, inv, Y)
    return out / n.reshape((-1,) + (1,) * (Y.ndim - 1)), n


def setup(sector, mode, jitter_from="clean"):
    C = Cutout(sector); W = json.load(open(f"{PH}/out/wcs_s{sector}.json"))
    if jitter_from == "phot":
        J = dict(np.load(f"{PH}/data/jitter_s{sector}.npz"))
    else:
        fn = f"{CL}/data/jitter_{mode}_s{sector}.npz"
        J = dict(np.load(fn)) if os.path.exists(fn) else jitter_clean(sector, mode, C)
    g, rBA = build_stars(C, W)
    prf = PRFModel(C); prf.blur = W["blur"]
    vm = var_mask(sector, mode, C.ny, C.nx)
    return dict(C=C, W=W, J=J, g=g, prf=prf, rBA=rBA, vm=vm, mode=mode)


def var_aperture_lc(C, sector, source_id, r=1.6):
    """Simple background-subtracted aperture flux of a field variable from the raw cube (for explicit regressors)."""
    v = [v for v in VARPOS[str(sector)]["vars"] if v["source_id"] == source_id][0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    d = np.hypot(xx - v["x"], yy - v["y"])
    ap = d <= r; ann = (d > 3) & (d <= 5)
    fl = C.flux.reshape(len(C.t), -1)
    if ap.sum() == 0: return None
    bk = np.nanmedian(fl[:, ann.ravel()], 1) if ann.sum() > 3 else np.zeros(len(C.t))
    return fl[:, ap.ravel()].sum(1) - ap.sum() * bk


def photometry(sector, S, inject=None, extra_mask_times=(), rfit=3.0, extra_regs=None, return_D=False,
               fit_models=("B", "A")):
    """Copy of phot p05.photometry (model B, Bp, A) with
    CLEAN: faint-pixel background plane computed from non-variable pixels only;
    CLEAN: pointing/background regressors from jitter_clean (variables excluded);
    CLEAN: optional extra per-cadence regressors (dict name -> array) in the pixel regression (nuisance, removed).
    inject: list of dict(img=2-D image (unit sum), amp=array e-/s per cadence) subtracted from the cube."""
    C, W, J, g, prf, vm = S["C"], S["W"], S["J"], S["g"], S["prf"], S["vm"]
    cube = C.flux.copy()
    mira = g[g.name == "Mira"].iloc[0]; m3v = g[g.name == "M3V"].iloc[0]; pair = g[g.name == "pair"].iloc[0]
    if inject:
        for inj in inject:
            cube -= inj["amp"][:, None, None] * inj["img"][None]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    fin = C.fin & np.isfinite(J["dx"]) & np.isfinite(J["dy"]) & np.isfinite(J["bkg"])
    q0 = fin & (C.q == 0)
    qok = fin & ((C.q == 0) | (C.q == 2048))
    B = J["bkg"]
    jit_ok = np.ones(len(C.t), bool)
    for k in ("dx", "dy"):
        v = J[k]; jit_ok &= np.abs(v - np.nanmedian(v[q0])) < 6 * robust_std(v[q0]) + 0.02
    win = in_windows(C.t, transit_times(C.t, P_EPH, T0_EPH), 0.1)
    if len(extra_mask_times): win |= in_windows(C.t, extra_mask_times, 0.1)
    t = C.t
    flat = cube.reshape(len(t), -1)
    med = np.nanmedian(flat[q0], 0)
    ok = ~vm.ravel()                                                                # CLEAN
    faint = ok & (med < np.percentile(med[ok], 15))
    uf = ((xx.ravel()[faint]) - 10) / 10; vf = ((yy.ravel()[faint]) - 10) / 10
    Af = np.column_stack([np.ones(faint.sum()), uf, vf])
    BG = np.full((len(t), 3), np.nan)
    cf, *_ = np.linalg.lstsq(Af, flat[fin][:, faint].T, rcond=None); BG[fin] = cf.T
    uall = (xx.ravel() - 10) / 10; vall = (yy.ravel() - 10) / 10
    flat = flat - (BG[:, 0:1] + BG[:, 1:2] * uall[None] + BG[:, 2:3] * vall[None])
    bexc = B - np.median(B[q0])
    reg = qok & (bexc < 150.0) & jit_ok & ~win
    ap = (np.hypot(xx - mira.x, yy - mira.y) <= 2.5).ravel()
    Fap = np.full(len(t), np.nan)
    Fap[fin] = flat[fin][:, ap].sum(1)
    mu = dict(F=np.mean(Fap[reg]), B=np.mean(BG[reg, 0]), dx=np.mean(J["dx"][reg]), dy=np.mean(J["dy"][reg]))
    ER = {}
    if extra_regs:
        for k_, v_ in extra_regs.items():
            v_ = np.asarray(v_, float); ER[k_] = (v_, np.nanmean(v_[reg]))
            fin &= np.isfinite(v_)
        reg &= fin
    def regressors(sel):
        ddx, ddy = J["dx"][sel] - mu["dx"], J["dy"][sel] - mu["dy"]
        cols = [np.ones(sel.sum()), Fap[sel] - mu["F"], BG[sel, 0] - mu["B"], ddx, ddy]
        for k_, (v_, m_) in ER.items(): cols.append(v_[sel] - m_)
        return np.column_stack(cols)
    Xb, n = binned(t[reg], regressors(reg), 1 / 24); Yb, _ = binned(t[reg], flat[reg], 1 / 24)
    kb = n >= max(2, 0.5 * np.median(n)); Xb, Yb = Xb[kb], Yb[kb]
    coef, *_ = np.linalg.lstsq(Xb, Yb, rcond=None)
    a_img, b_img = coef[0].reshape(C.ny, C.nx), coef[1].reshape(C.ny, C.nx)
    D = np.full(flat.shape, np.nan)
    cz = coef.copy(); cz[1] = 0.0
    D[fin] = flat[fin] - regressors(fin) @ cz
    fm = (np.hypot(xx - mira.x, yy - mira.y) <= rfit) | (np.hypot(xx - m3v.x, yy - m3v.y) <= rfit)
    fmr = fm.ravel()
    TM = b_img.copy()
    T3 = prf.image(m3v.x, m3v.y); TP = prf.image(pair.x, pair.y)
    u = (xx - mira.x) / 3; v = (yy - mira.y) / 3
    base = dict(Mira=TM, M3V=T3, pair=TP, bkg0=np.ones_like(u), bkgx=u, bkgy=v)
    MODELS = dict(B=["Mira", "M3V", "bkg0", "bkgx", "bkgy"], Bp=["Mira", "pair", "bkg0", "bkgx", "bkgy"],
                  A=["Mira", "M3V", "pair", "bkg0", "bkgx", "bkgy"])
    Dg = D[q0 & fin][:, fmr]
    sig = 1.4826 * np.nanmedian(np.abs(np.diff(Dg, axis=0)), axis=0) / np.sqrt(2) + 0.5
    Y = (D[fin][:, fmr] / sig).T
    FIT = {}
    for mk in fit_models:
        nms = MODELS[mk]
        A = np.column_stack([base[n_].ravel()[fmr] for n_ in nms])
        Aw = A / sig[:, None]
        cov = np.linalg.pinv(Aw.T @ Aw)
        F = np.full((len(t), A.shape[1]), np.nan); chi2 = np.full(len(t), np.nan)
        c_ = cov @ Aw.T @ Y
        F[fin] = c_.T
        chi2[fin] = ((Aw @ c_ - Y) ** 2).sum(0) / (fmr.sum() - A.shape[1])
        FIT[mk] = dict(F=F, names=nms, ferr=np.sqrt(np.diag(cov)), chi2=chi2)
    out = dict(t=t, q=C.q, cadn=C.cadn, timecorr=C.timecorr, fin=fin, q0=q0, qok=qok, reg=reg, bkg=B, jit_ok=jit_ok,
               dx=J["dx"], dy=J["dy"], Fap=Fap, BG=BG, FIT=FIT, b_img=b_img, a_img=a_img, sig=sig, mu=mu,
               exptime=C.exptime, coef=coef, fitmask=fm,
               pos={k: (float(r_.x), float(r_.y)) for k, r_ in (("Mira", mira), ("M3V", m3v), ("pair", pair))})
    if return_D: out["D"] = D.reshape(len(t), C.ny, C.nx)
    return out


# ------------------------------------------------------------------------------------------------ light curves
def detrend(t, f, ok, mask_tr, window=0.75):
    m = ok & np.isfinite(f) & ~mask_tr
    tr = np.full(len(t), np.nan)
    if m.sum() < 20: return tr
    trm = flatten(t[m], f[m], method="biweight", window_length=window, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
    k = np.isfinite(trm)
    fin = np.isfinite(t)
    tr[fin] = np.interp(t[fin], t[m][k], trm[k])
    tk = t[m][k]
    idx = np.clip(np.searchsorted(tk, t[fin]), 1, len(tk) - 1)
    dist = np.minimum(np.abs(t[fin] - tk[idx - 1]), np.abs(t[fin] - tk[idx]))
    tr_f = tr[fin]; tr_f[dist > 0.5] = np.nan; tr[fin] = tr_f
    return tr


_TC_CACHE = {}
def t_tdb_of(sector, t, timecorr):
    """Target-barycentred BTJD_TDB (phot p00/p06 correction, JPL Horizons cache in ../phot/data/horizons)."""
    key = (sector, len(t))
    if key not in _TC_CACHE:
        add = np.zeros(len(t)); m = np.isfinite(t) & np.isfinite(timecorr)
        add[m] = roemer_target(sector, t[m], timecorr[m])
        _TC_CACHE[key] = add
    return t + _TC_CACHE[key] / 86400.0


def in_transit_windows(t, P=P1, T0=T01, half=0.1):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    return np.abs(ph * P) < half


def make_lc(sector, R, extra_mask=None, models=(("B", "M3V"), ("A", "M3V"), ("B", "Mira"), ("A", "pair"), ("A", "Mira"))):
    """p06.process equivalent. quality_ok: finite & QUALITY in {0, 2048} & bkg excess < 60 & pointing ok & chi2_B < 3x median.
    Normalisation: M3V 979 e-/s (PRF units); detrending wotan biweight 0.75 d with +-0.1 d transit windows masked."""
    t = R["t"]; t_tdb = t_tdb_of(sector, t, R["timecorr"])
    chiB = R["FIT"]["B"]["chi2"]
    bexc = R["bkg"] - np.median(R["bkg"][R["q0"]])
    base_ok = R["qok"] & (bexc < 60.0) & R["jit_ok"] & np.isfinite(chiB)
    base_ok &= chiB < 3 * np.nanmedian(chiB[R["q0"]])
    mask_tr = in_transit_windows(t_tdb)
    if extra_mask is not None: mask_tr = mask_tr | extra_mask
    fref = dict(FREF_INV, Mira=float(R["mu"]["F"] - 0.9 * (tmag_flux(13.0135) + FREF_INV["pair"])))
    out = dict(t_tdb=t_tdb, t_tess=t, q=R["q"], quality_ok=base_ok, bkg=R["bkg"], dx=R["dx"], dy=R["dy"],
               intransit_window=in_transit_windows(t_tdb), exptime=R["exptime"])
    for mk, nm in models:
        if mk not in R["FIT"]: continue
        names = list(R["FIT"][mk]["names"])
        dF = R["FIT"][mk]["F"][:, names.index(nm)]
        raw = 1.0 + dF / fref[nm]
        tr = detrend(t_tdb, raw, base_ok, mask_tr)
        det = raw / tr
        ok = base_ok & np.isfinite(det)
        sd = robust_std(det[ok & ~mask_tr])
        chi = R["FIT"][mk]["chi2"]; scale = np.sqrt(np.clip(chi / np.nanmedian(chi[ok]), 1.0, None))
        key = f"{nm}_{mk}"
        out[f"dF_{key}"] = dF; out[f"raw_{key}"] = raw; out[f"trend_{key}"] = tr; out[f"det_{key}"] = det
        out[f"err_{key}"] = sd * scale; out[f"rstd_{key}"] = sd
    return out


# ------------------------------------------------------------------------------------------------ analysis helpers
def harm_fit(t, y, w, Pp, T0=TMAX_RR, nh=3, extra=None):
    """Weighted fit of nh harmonics; returns first-harmonic semi-amplitude, projection on cos (max at T0), all coefs."""
    ph = 2 * np.pi * (t - T0) / Pp
    cols = [np.ones(len(t))] + [f(h * ph) for h in range(1, nh + 1) for f in (np.cos, np.sin)]
    if extra is not None: cols += list(extra)
    X = np.column_stack(cols) * np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(X, y * np.sqrt(w), rcond=None)
    return float(np.hypot(c[1], c[2])), float(c[1]), c


def trap_frac(t, tc, T14, fl, exptime, nsub=None):
    if nsub is None: nsub = 1 if exptime < 0.003 else (5 if exptime < 0.008 else 11)
    y = np.zeros(np.shape(t)); T23 = fl * T14
    for o in (np.arange(nsub) + 0.5) / nsub - 0.5:
        a = np.abs(t + o * exptime - tc)
        y = y + np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-9)))
    return y / nsub


def jd(obj, fn):
    save_json(obj, fn)
