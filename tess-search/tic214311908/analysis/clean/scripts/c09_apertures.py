"""c09: transit depth per sector under different photometric weightings of the SAME cleaned pixel data:
PRF-weighted fits (INV model B, REF 'M') versus small apertures of 1, 2, 4 pixels (the M3V-PRF-brightest pixels) and
r <= 1.5 px around the M3V.

INV apertures: difference cube D of the cleaned INV pipeline (mode 'all'), minus the per-cadence model-A Mira (b image),
pair (PRF) and plane terms -> M3V + residual; aperture sum divided by the M3V PRF fraction inside the aperture ->
M3V flux change in PRF units; normalised by 979 e-/s; then the same +-0.3 d window fits (fixed shape, free depth at the
predicted mid-time, linear baseline) as for the PRF light curve.
REF apertures: residual frames Rz of the referee-method windows (per-pixel baselines 1, tau, tau^2, B, dx, dy fitted out
of transit: the Mira, the pair and the background are removed per pixel and per window), summed over the aperture,
divided by the PRF fraction, normalised by 979 e-/s.
Errors: the same measurement at fake epochs (INV: 60 per sector; REF: 60 per sector), std of single-window depths ->
per-sector error = that / sqrt(n_real windows) (empirical).
Outputs: out/apertures.json, out/apertures_per_transit.csv
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from refclean import SectorC
from c05_inject import fit_window, accept
import pandas as pd
import r05_local as R5

rng = np.random.default_rng(909)
EPH = json.load(open(f"{CL}/out/ephemeris_final.json"))["combined"]
P, T0 = EPH["P_d"], EPH["T0_epoch0_BTJD_TDB"]


def apertures(T3):
    order = np.argsort(T3.ravel())[::-1]
    ap = {}
    for n in (1, 2, 4):
        m = np.zeros(T3.size, bool); m[order[:n]] = True; ap[f"px{n}"] = m.reshape(T3.shape)
    return ap


def real_and_fake(t, n_fake=60):
    tf = t[np.isfinite(t)]
    E = np.arange(np.ceil((tf.min() - T0) / P - 0.1), np.floor((tf.max() - T0) / P + 0.1) + 1).astype(int)
    real = [(T0 + e * P, int(e)) for e in E]
    fake = []
    tries = 0
    while len(fake) < n_fake and tries < 20000:
        tries += 1
        tc = rng.uniform(tf.min() + 0.3, tf.max() - 0.3)
        if abs((((tc - T0) / P + 0.5) % 1 - 0.5) * P) < 0.45: continue
        if fake and min(abs(f[0] - tc) for f in fake) < 0.15: continue
        fake.append((tc, -9999))
    return real, fake


def inv_sector(s):
    S = setup(s, "all")
    R = photometry(s, S, return_D=True, fit_models=("A", "B"))
    L = make_lc(s, R, models=(("B", "M3V"),))
    D = R["D"]; F = R["FIT"]["A"]["F"]          # names: Mira, M3V, pair, bkg0, bkgx, bkgy
    yy, xx = np.mgrid[0:21, 0:21]
    mx, my = R["pos"]["Mira"]; px, py = R["pos"]["pair"]; m3x, m3y = R["pos"]["M3V"]
    u = (xx - mx) / 3; v = (yy - my) / 3
    prf = S["prf"]; T3 = prf.image(m3x, m3y); TP = prf.image(px, py)
    other = (F[:, 0, None, None] * R["b_img"][None] + F[:, 2, None, None] * TP[None] + F[:, 3, None, None]
             + F[:, 4, None, None] * u[None] + F[:, 5, None, None] * v[None])
    Rz = D - other
    aps = apertures(T3); aps["r1.5"] = np.hypot(xx - m3x, yy - m3y) <= 1.5
    lcs = {"prf_B": L["raw_M3V_B"], "prf_A": 1 + F[:, 1] / F_M3V}
    for k, m in aps.items():
        frac = T3[m].sum()
        lcs[k] = 1 + Rz[:, m].sum(1) / frac / F_M3V
    t = L["t_tdb"]; ok = L["quality_ok"]
    real, fake = real_and_fake(t[ok])
    rows = []
    for kind, cen in (("real", real), ("fake", fake)):
        for tc, ep in cen:
            m = ok & (np.abs(t - tc) < 0.3)
            if m.sum() < 8 or not accept(t[m], tc, 0.3, R["exptime"]): continue
            for k, f in lcs.items():
                if not np.all(np.isfinite(f[m])): continue
                r = fit_window(t[m], f[m], tc, R["exptime"], free=False)
                if r: rows.append(dict(pipe="inv", sector=s, kind=kind, epoch=ep, tc=tc, weighting=k, depth=r["depth_fix"]))
    return rows, {k: float(T3[m].sum()) for k, m in aps.items()}


def ref_sector(s):
    S = SectorC(s, "all")
    m3x, m3y = S.pos["M3V"]
    T3 = S.P_M3V
    yy, xx = np.mgrid[0:21, 0:21]
    aps = apertures(T3); aps["r1.5"] = np.hypot(xx - m3x, yy - m3y) <= 1.5
    real, fake = real_and_fake(S.ttdb[S.qok])
    rows = []
    regidx = np.where(S.reg)[0]
    for kind, cen in (("real", real), ("fake", fake)):
        for tc, ep in cen:
            w = S.window(tc, keep_cube=True)
            if w is None or not accept(w["t"], tc, {12: 0.6, 39: 0.4}.get(s, 0.35), S.exptime): continue
            lcs = {"prf_M": 1 + w["aM"] / F_M3V, "prf_T": 1 + w["aT"][:, 0] / F_M3V}
            for k, m in aps.items():
                mm = m.ravel()[regidx]
                if mm.sum() < m.sum(): continue          # aperture must lie in the referee's pixel set
                lcs[k] = 1 + w["Rz"][:, mm].sum(1) / T3[m].sum() / F_M3V
            for k, f in lcs.items():
                r = fit_window(w["t"], f, tc, S.exptime, free=False)
                if r: rows.append(dict(pipe="ref", sector=s, kind=kind, epoch=ep, tc=tc, weighting=k, depth=r["depth_fix"]))
    return rows, {k: float(T3[m].sum()) for k, m in aps.items()}


def main():
    if "--from-csv" in sys.argv:
        d = pd.read_csv(f"{CL}/out/apertures_per_transit.csv"); fr = json.load(open(f"{CL}/out/aperture_prf_fractions.json"))
    else:
        rows = []; fr = {}
        for s in SECTORS:
            r1, f1 = inv_sector(s); r2, f2 = ref_sector(s)
            rows += r1 + r2; fr[s] = dict(inv=f1, ref=f2)
            print(f"S{s} done: PRF fractions INV {f1}, REF {f2}", flush=True)
        d = pd.DataFrame(rows)
        d.to_csv(f"{CL}/out/apertures_per_transit.csv", index=False, float_format="%.4f")
        jd(fr, f"{CL}/out/aperture_prf_fractions.json")
    out = dict(prf_fraction=fr, table={})
    for (pipe, wgt), g in d.groupby(["pipe", "weighting"]):
        res = {}
        for s in SECTORS:
            gr = g[(g.sector == s) & (g.kind == "real")]; gf = g[(g.sector == s) & (g.kind == "fake")]
            if len(gr) == 0 or len(gf) < 10: continue
            e1 = float(np.std(gf.depth)); res[s] = dict(depth=float(gr.depth.mean()), err=e1 / np.sqrt(len(gr)), n=len(gr), single_err=e1,
                                                        fake_mean=float(gf.depth.mean()))
        secs = sorted(res)
        v = np.array([res[s]["depth"] for s in secs]); e = np.array([res[s]["err"] for s in secs])
        w = 1 / e ** 2; mu = float(np.sum(w * v) / np.sum(w))
        early = [s for s in secs if s in (39, 66)]; late = [s for s in secs if s in (103, 104)]
        def wm(ss):
            vv = np.array([res[s]["depth"] for s in ss]); ee = np.array([res[s]["err"] for s in ss]); ww = 1 / ee ** 2
            return float(np.sum(ww * vv) / np.sum(ww)), float(1 / np.sqrt(np.sum(ww)))
        a, ea = wm(early); b, eb = wm(late)
        out["table"][f"{pipe}:{wgt}"] = dict(per_sector=res, mean=mu, mean_err=float(1 / np.sqrt(np.sum(w))),
                                             chi2=float(np.sum(w * (v - mu) ** 2)), late_minus_early=[b - a, float(np.hypot(ea, eb))])
        print(f"{pipe}:{wgt:6s} " + " ".join(f"S{s} {res[s]['depth']:5.1f}+-{res[s]['err']:.1f}" for s in secs) +
              f" | mean {mu:.2f}+-{1/np.sqrt(np.sum(w)):.2f} chi2 {np.sum(w*(v-mu)**2):.1f}; late-early {b-a:+.1f}+-{np.hypot(ea,eb):.1f}", flush=True)
    # paired differences aperture - PRF per sector (same events: error from fake-epoch paired differences)
    out["paired"] = {}
    for pipe, ref_w in (("inv", "prf_B"), ("ref", "prf_M")):
        g = d[d["pipe"] == pipe]
        for wgt in g.weighting.unique():
            if wgt == ref_w: continue
            res = {}
            for s in SECTORS:
                def diffs(kind):
                    a = g[(g.sector == s) & (g.kind == kind) & (g.weighting == wgt)].set_index("tc").depth
                    b = g[(g.sector == s) & (g.kind == kind) & (g.weighting == ref_w)].set_index("tc").depth
                    j = pd.concat([a, b], axis=1, join="inner"); return (j.iloc[:, 0] - j.iloc[:, 1]).values
                dr = diffs("real"); df = diffs("fake")
                if len(dr) == 0 or len(df) < 10: continue
                res[s] = dict(diff=float(dr.mean()), err=float(np.std(df) / np.sqrt(len(dr))), ratio_fake_mean=float(df.mean()))
            secs = sorted(res)
            v = np.array([res[s]["diff"] for s in secs]); e = np.array([res[s]["err"] for s in secs]); w = 1 / e ** 2
            mu = float(np.sum(w * v) / np.sum(w))
            out["paired"][f"{pipe}:{wgt}-{ref_w}"] = dict(per_sector=res, mean=mu, mean_err=float(1 / np.sqrt(w.sum())), chi2=float(np.sum(w * (v - mu) ** 2)), dof=len(secs) - 1)
            print(f"paired {pipe}:{wgt}-{ref_w}: " + " ".join(f"S{s} {res[s]['diff']:+.1f}+-{res[s]['err']:.1f}" for s in secs) + f" | mean {mu:+.2f}, chi2 {np.sum(w*(v-mu)**2):.1f}/{len(secs)-1}", flush=True)
    jd(out, f"{CL}/out/apertures.json")


if __name__ == "__main__":
    main()
