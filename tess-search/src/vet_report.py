"""Full vetting dossier for one candidate: writes plots, a JSON of numbers and a Markdown summary."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
import stage2, vet, known
from centroid import localise
from quickplot import quick

W = "/home/user/work"


def j(x):
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, tuple): return list(x)
    return x


def run(tic, P, t0, dur, depth, star, outdir):
    os.makedirs(outdir, exist_ok=True)
    rep = dict(tic=tic, star=star)
    secs = stage2.load_all(tic)
    rep["sectors"] = [(s["sector"], s["prov"], s["exptime"], float(s["sigma"])) for s in secs]
    jf = stage2.joint_fit(secs, P, dur, rel=0.003)
    rep["joint_bls"] = jf
    f = vet.fit(secs, jf["P"], jf["t0"], jf["dur"], max(jf["depth"], depth / 2))
    rep["fit"] = f
    rad, mass = star.get("rad"), star.get("mass")
    rep["rp_earth"] = f["rp"] * (rad or np.nan) * 109.1
    rep["rho_fit"] = vet.rho_star(f["a"], f["P"])
    rep["rho_tic"] = (mass / rad ** 3) if (rad and mass) else np.nan
    # per-sector depth at the final ephemeris
    per = []
    from bls import diagnostics
    for s in secs:
        d = diagnostics(s["t"], s["f"], f["P"], f["t0"], jf["dur"], s["sigma"])
        per.append(dict(sector=s["sector"], prov=s["prov"], depth_ppm=d["depth"] * 1e6, err_ppm=d["depth_err"] * 1e6,
                        n_tr=d["n_tr"], oddeven_sig=d["oddeven_sig"], sec_sig=d["sec_sig"], oot_rms=d["oot_rms_ratio"]))
    rep["per_sector"] = per
    allt = np.concatenate([s["t"] for s in secs]); allf = np.concatenate([s["f"] for s in secs])
    sig = np.concatenate([np.full(len(s["t"]), s["sigma"]) for s in secs])
    d = diagnostics(allt, allf, f["P"], f["t0"], jf["dur"], float(np.median(sig)))
    rep["combined"] = d
    rep["held_out"] = vet.held_out(secs, f["P"], jf["dur"])
    # neighbours
    try:
        g = vet.gaia_neighbours(star["ra"], star["dec"], star["Tmag"], max(f["rp"] ** 2, 1e-5))
        rep["gaia"] = g[["source_id", "sep", "phot_g_mean_mag", "phot_rp_mean_mag", "tmag_est", "ruwe", "parallax", "can_mimic"]].head(25).to_dict("records")
    except Exception as e:
        rep["gaia_error"] = str(e)[:200]
    try:
        rep["vizier"] = known.vizier_cone(star["ra"], star["dec"], 150)
    except Exception as e:
        rep["vizier_error"] = str(e)[:200]
    rep["known_on_star"] = known.known_signals(tic, f["P"])[0]
    # localisation in up to three sectors with the most transits (FFI cutouts)
    loc = []
    for s in sorted(secs, key=lambda s: -s["sector"])[:3]:
        try:
            L = localise(star["ra"], star["dec"], s["sector"], f["P"], f["t0"], jf["dur"])
        except Exception as e:
            L = None
        if L:
            fig, ax = plt.subplots(1, 2, figsize=(8, 4))
            for a, img, ttl in ((ax[0], L["direct_img"], "direct"), (ax[1], L["diff_img"], "difference (out - in)")):
                a.imshow(img, origin="lower", cmap="viridis"); a.plot(*L["target_px"], "r+", ms=12); a.set_title(f"S{s['sector']} {ttl}", fontsize=9)
            ax[1].plot(*L["diff_centroid"], "wx", ms=10)
            if rep.get("gaia"):
                for gg in rep["gaia"][:15]:
                    pass
            fig.tight_layout(); fig.savefig(f"{outdir}/diffimg_s{s['sector']}.png", dpi=90); plt.close(fig)
            loc.append({k: v for k, v in L.items() if k not in ("diff_img", "direct_img", "wcs")})
    rep["localisation"] = loc
    vet.plot(tic, secs, f["P"], f["t0"], jf["dur"], f, f"{outdir}/transits.png",
             title=f"TIC {tic}  P={f['P']:.5f} d  depth={f['rp']**2*1e6:.0f} ppm  Rp={rep['rp_earth']:.2f} Re")
    quick(tic, f["P"], f["t0"], jf["dur"], f"{outdir}/phase_periodogram.png", secs=secs)
    json.dump(rep, open(f"{outdir}/report.json", "w"), default=j, indent=1)
    return rep


if __name__ == "__main__":
    tic = int(sys.argv[1]); sector = sys.argv[2]  # a sector number, or "multi" for the combined S104+S105 search
    c = pd.read_csv(f"{W}/results/new_multi.csv" if sector == "multi" else f"{W}/results/new_s{int(sector):04d}.csv")
    r = c[c.tic == tic].iloc[0].to_dict()
    star = {k: (float(r[k]) if pd.notna(r[k]) else None) for k in ("Tmag", "Teff", "rad", "mass", "d", "ra", "dec", "contratio")}
    rep = run(tic, r["P"], r["t0"], r["dur"], r["m_depth"], star, f"{W}/vet/TIC{tic}")
    print(json.dumps({k: rep[k] for k in ("fit", "rp_earth", "rho_fit", "rho_tic", "combined")}, default=j, indent=1))
