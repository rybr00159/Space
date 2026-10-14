"""Per-sector stacked eclipse times (primary and secondary separately) with EMPIRICAL errors from bootstrapping events.
For each event: chi2(shift) curve of the best-fit physical template (cadence-integrated) with free depth scale and quadratic
baseline, on a 0.05-min grid over +-15 min.  Sector time = argmin of the summed curve (parabolic refinement); error = robust
sd of argmins over 2000 bootstrap resamplings of the sector's events.  Then linear / quadratic ephemeris fits to the 12 sector
times, and a quadratic positive control at this level (inject c2*tyr^2 into the sector times with bootstrap-error noise).
This is an independent aggregation of the same photometry (no per-event errors, no z-scale calibration involved).
Usage: s16_sector_times.py <params.json>"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import lclib as L, etvlib as E
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
pb = json.load(open(sys.argv[1])); P, T0 = pb['P'], pb['T0']
SH = np.arange(-15, 15.0001, 0.05)
def curve(t, y, tc, texp, nsub, scale):
    dt = t - tc
    tf = np.arange(dt.min() - 16 / 1440 - 0.01, dt.max() + 16 / 1440 + 0.01, 10 / 86400.0)
    dfl = 1 - E.template(pb, tf + tc, texp, nsub, scale)
    X0 = np.vstack([np.ones_like(dt), dt, dt * dt]).T
    out = np.empty(SH.size)
    for i, s in enumerate(SH / 1440):
        A = np.column_stack([X0, -np.interp(dt - s, tf, dfl)])
        c, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ c; out[i] = r @ r
    return out
def argmin_par(c):
    j = int(np.argmin(c))
    if j < 3 or j > c.size - 4: return np.nan
    sl = slice(j - 6, j + 7); cc = np.polyfit(SH[sl], c[sl], 2)
    return -cc[1] / (2 * cc[0]) if cc[0] > 0 else np.nan
rng = np.random.default_rng(17)
rows = []
for lc in L.default_set():
    sec = lc['sector']; lc = L.clip_outliers(lc, P=P, T0=T0); scale = pb.get(f's{sec}', 1.0)
    ws = L.windows(lc, P=P, T0=T0, hw=0.20)
    res = []
    for w in ws:
        i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
        cc = np.polyfit(dt[oot], y[oot], 2); res.append(y[oot] / np.polyval(cc, dt[oot]) - 1)
    r = np.concatenate(res); sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    for ph in (0.0, 0.5):
        C = []; nn = []
        for w in ws:
            if w['ph'] != ph: continue
            i = w['idx']; tcm = T0 + (w['n'] + ph) * P
            C.append(curve(lc['t'][i], lc['f'][i], tcm, lc['texp'], lc['nsub'], scale) / sig ** 2); nn.append(w['n'])
        C = np.array(C); x0 = argmin_par(C.sum(0))
        bs = []
        for b in range(2000):
            k = rng.integers(0, len(C), len(C)); bs.append(argmin_par(C[k].sum(0)))
        bs = np.array(bs); bs = bs[np.isfinite(bs)]
        err = 1.4826 * np.median(np.abs(bs - np.median(bs)))
        cs = C.sum(0); j = np.argmin(cs); m = np.abs(SH - x0) < 3
        cc = np.polyfit(SH[m], cs[m], 2); ef = 1 / np.sqrt(cc[0])
        Emid = float(np.mean(nn)) + ph
        rows.append(dict(sector=sec, ph=ph, n_ev=len(C), E=Emid, t=float(T0 + Emid * P), dt_min=float(x0), err_boot_min=float(err),
                         err_formal_min=float(ef), err_min=float(max(err, ef))))
        print(f'S{sec:3d} {"prim" if ph == 0 else "sec "}: N={len(C):2d} O-C {x0:+.3f} min; bootstrap err {err:.3f}; formal {ef:.3f}', flush=True)
E_ = np.array([r['E'] for r in rows]); oc = np.array([r['dt_min'] for r in rows]); es = np.array([r['err_min'] for r in rows])
ph = np.array([r['ph'] for r in rows]); tyr = (np.array([r['t'] for r in rows]) - 2900) / 365.25
def wls(A, y, e):
    Aw = A / e[:, None]; c, *_ = np.linalg.lstsq(Aw, y / e, rcond=None); rr = (y - A @ c) / e
    return c, float(rr @ rr), np.linalg.inv(Aw.T @ Aw)
Al = np.array([np.ones_like(E_), E_ / 1000, (ph == 0.5).astype(float)]).T; Aq = np.column_stack([Al, tyr ** 2])
cl, chl, covl = wls(Al, oc, es); cq, chq, covq = wls(Aq, oc, es)
out = dict(rows=rows, linear=dict(chi2=chl, dof=len(oc) - 3, dP_d=float(cl[1] / 1000 / 1440), dP_err_d=float(np.sqrt(covl[1, 1]) / 1000 / 1440),
                                  dT0_min=float(cl[0]), sec_offset_min=float(cl[2]), sec_offset_err=float(np.sqrt(covl[2, 2]))),
           quad=dict(chi2=chq, c2=float(cq[3]), c2_err=float(np.sqrt(covq[3, 3])), dchi2=chl - chq))
print(f'linear: chi2 {chl:.2f}/{len(oc)-3}; P = {P + cl[1]/1000/1440:.8f} +- {np.sqrt(covl[1,1])/1000/1440:.8f} d; secondary offset {cl[2]:+.2f} +- {np.sqrt(covl[2,2]):.2f} min')
print(f'quadratic: c2 = {cq[3]:+.4f} +- {np.sqrt(covq[3,3]):.4f} min/yr^2; dchi2 {chl-chq:.2f}')
# positive control at sector level
pc = []
for c2 in (0.0, 0.05, 0.1, 0.2):
    rec = []
    for b in range(2000):
        y = Al @ cl + c2 * tyr ** 2 + rng.normal(0, es)
        cqb, _, _ = wls(Aq, y, es); rec.append(cqb[3])
    rec = np.array(rec); pc.append(dict(c2_inj=c2, c2_rec_med=float(np.median(rec)), c2_rec_sd=float(np.std(rec))))
    print(f'  sector-level control: c2 inj {c2:.2f} -> {np.median(rec):+.4f} +- {np.std(rec):.4f}')
out['positive_control'] = pc
json.dump(out, open(f'{FD}/results/sector_times.json', 'w'), indent=1)
