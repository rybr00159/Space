"""Referee r16: systematics of the referee flux loss and odd/even difference from method choices in r05:
fit radius around the M3V (1.5 / 2.5 / 3.5 px), second-order pointing terms, window half-width (0.2 / 0.25 / 0.35 / 0.5 d),
out-of-transit mask (0.06 / 0.07 / 0.09 d), linear instead of quadratic time baseline (linbase), and the PRF template displaced
by 0.1 px. Per variant: per-sector flux loss at the r06 global shape and the all-sector weighted mean (errors from 60
fake windows per sector), and the odd-even difference. Output: out/variants.csv, out/variants.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
import r05_local as L
from r05_local import Sector, trap_shape, real_centres, fake_centres
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]
VAR = {"default": {}, "rM1.5": dict(rM=1.5), "rM3.5": dict(rM=3.5), "order2": dict(order2=True), "half0.25": dict(half=0.25), "half0.5": dict(half=0.5),
       "mask0.06": dict(mask=0.06), "mask0.09": dict(mask=0.09), "half0.2": dict(half=0.2), "linbase": dict(linbase=True), "shift0.1": dict(shift=(0.1, -0.1)), "shift-0.1": dict(shift=(-0.1, 0.1))}


def depth(r, tc, ex):
    y = trap_shape(r["t"], tc, T14, FL, ex); sg = rstd(r["aM"][r["oot"]])
    A = np.column_stack([np.ones(len(y)), -y]); c, *_ = np.linalg.lstsq(A / sg, r["aM"] / sg, rcond=None); return c[1]


rows = []
for s in SECTORS:
    S0 = Sector(s); cen, ks = real_centres(S0, P, T0); fk = fake_centres(S0, 60, np.random.default_rng(50 + s), P, T0)
    for vn, kw in VAR.items():
        S = Sector(s)
        if "rM" in kw:
            m3x, m3y = S.pos["M3V"]; S.regM = (np.hypot(S.xx - m3x, S.yy - m3y) <= kw["rM"]).ravel()
            S.reg = S.reg | S.regM
        L.MASK = kw.get("mask", 0.07)
        tm = None
        if "shift" in kw: tm = S.C.psf(S.pos["M3V"][0] + kw["shift"][0], S.pos["M3V"][1] + kw["shift"][1])
        half = kw.get("half", None)
        dd = []; ee = []
        for tc, k in zip(cen, ks):
            if s == 12 and k == -388: continue
            r = S.window(tc, half=half, order2=kw.get("order2", False), tmpl_M=tm, bmax=60.0, linbase=kw.get("linbase", False))
            if r is None: continue
            dd.append((k, depth(r, tc, S.exptime)))
        nul = []
        for tc in fk:
            r = S.window(tc, half=half, order2=kw.get("order2", False), tmpl_M=tm, linbase=kw.get("linbase", False))
            if r is None: continue
            nul.append(depth(r, tc, S.exptime))
        e1 = float(np.std(nul))
        dd = np.array(dd)
        rows.append(dict(sector=s, variant=vn, n=len(dd), dF=float(dd[:, 1].mean()), err=e1 / np.sqrt(len(dd)), single_err=e1,
                         odd=float(dd[dd[:, 0] % 2 != 0, 1].mean()) if np.any(dd[:, 0] % 2 != 0) else np.nan, n_odd=int(np.sum(dd[:, 0] % 2 != 0)),
                         even=float(dd[dd[:, 0] % 2 == 0, 1].mean()) if np.any(dd[:, 0] % 2 == 0) else np.nan, n_even=int(np.sum(dd[:, 0] % 2 == 0))))
    L.MASK = 0.07
    print(f"S{s}: " + ", ".join(f"{r['variant']} {r['dF']:.1f}+-{r['err']:.1f}" for r in rows if r["sector"] == s), flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{R}/out/variants.csv", index=False)
summ = {}
for vn in VAR:
    d = D[D.variant == vn]; w = 1 / d.err ** 2
    mu = float(np.sum(w * d.dF) / w.sum()); er = float(1 / np.sqrt(w.sum()))
    # odd-even with per-window error = single_err
    no = d.n_odd.values; ne = d.n_even.values; se = d.single_err.values
    ok = (no > 0); oke = (ne > 0)
    wo = no[ok] / se[ok] ** 2; we = ne[oke] / se[oke] ** 2
    mo = float(np.sum(wo * d.odd.values[ok]) / wo.sum()); me = float(np.sum(we * d.even.values[oke]) / we.sum())
    eoe = float(np.hypot(1 / np.sqrt(wo.sum()), 1 / np.sqrt(we.sum())))
    summ[vn] = dict(dF=mu, err=er, chi2_sectors=float(np.sum(w * (d.dF - mu) ** 2)), odd=mo, even=me, odd_minus_even=mo - me, oe_err=eoe)
    print(f"ALL {vn:10s}: dF {mu:.2f}+-{er:.2f} e/s (sector chi2 {summ[vn]['chi2_sectors']:.1f}/5); odd-even {mo-me:+.2f}+-{eoe:.2f} e/s")
vals = np.array([v["dF"] for v in summ.values()])
summ["spread"] = dict(min=float(vals.min()), max=float(vals.max()), std=float(vals.std()))
save_json(summ, f"{R}/out/variants.json")
