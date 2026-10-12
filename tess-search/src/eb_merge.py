"""Add a further batch of selected EBs (eb_selected.csv + an eb_final.py output) to the final table, using the same
recipe as the first 30. Usage: eb_merge.py vet/eb_final_b3.json vet/eb_neighbours_b3.json"""
import sys, json
import numpy as np, pandas as pd

W = "/home/user/work"
NOTES = {
    217276354: "near-twin; eccentric (secondary at phase 0.526)",
    406279367: "near-equal eclipses; odd and even depths differ in all 8 sectors, so the orbit is twice the search period",
    6598941: "F star + low-mass companion (EBLM-like): flat-bottomed eclipse, companion ~0.25 solar radii",
    63594264: "eclipse 4.3 % deep in S13/S39 but 2.9 % in S104: dilution differs between sectors",
    412733784: "slightly eccentric (secondary at phase 0.517)",
    90919279: "slightly eccentric (secondary at phase 0.513)",
    76440895: "well-separated eclipses; grouped with close binaries only because P < 0.8 d",
    33445288: "eclipses on top of ellipsoidal variation",
    # combined S104+S105 search
    100378891: "eccentric (secondary at phase 0.357, e cos w ~ -0.22)",
    139256801: "highly eccentric (secondary at phase 0.722, e cos w ~ +0.35)",
    404124262: "flat-bottomed (total) primary eclipse; cool companion (secondary 8 % of primary)",
    281740175: "flat-bottomed (total) primary eclipse; cool companion; secondary at phase 0.514",
    # final Sector 104 triage pass
    212554303: "M-dwarf pair; odd and even eclipses differ in all 6 sectors, so the orbit is twice the search period",
    316096209: "M-dwarf pair (3,100 K); odd and even eclipses differ in all 3 sectors, so the orbit is twice the search period",
    217228307: "out-of-eclipse modulation (spots or reflection)",
    379508029: "no secondary detected (< 0.04 %); the orbit could be twice this period if the two eclipses are identical",
    304770656: "flat-bottomed (total) primary eclipse; secondary at phase 0.489",
    343895671: "M-dwarf pair, near-twin",
    351474229: "eclipses on top of ellipsoidal variation",
    363482012: "M dwarf (3,800 K)",
}
# finer phase scan (0.0025 steps, wider window): secondary centred at 0.500, 570 +- 47 ppm
# 63594264 and 33445288: the automatic search returned the primary eclipse wing; the secondaries sit at phase 0.50
PHASE_FIX = {124358211: (0.0040, 0.503), 217228307: (0.0033, 0.487), 304770656: (0.0073, 0.489), 343895671: (0.2765, 0.501),
             363482012: (0.0673, 0.500), 379508029: (float("nan"), float("nan")),
             217276354: (0.287, 0.526), 412733784: (0.0051, 0.517), 6598941: (0.057e-2, 0.5), 63594264: (0.6e-2, 0.5), 33445288: (1.85e-2, 0.5)}


def main(final_json, neigh_json):
    sel = pd.read_csv(f"{W}/vet/eb_selected.csv")
    old = pd.read_csv(f"{W}/vet/new_ebs_final.csv")
    rows = {r["tic"]: r for r in json.load(open(final_json))}
    nb = json.load(open(neigh_json)) if neigh_json else {}
    out = []
    for x in sel[~sel.tic.isin(old.tic)].itertuples():
        r = rows.get(int(x.tic))
        if r is None: print("missing", x.tic); continue
        g = r["gaia"]
        fixed = abs(float(x.P) / r["P"] - 2) < 1e-3  # eb_select doubled the period
        P, Perr = (float(x.P), 2 * float(r["P_err"] or np.nan)) if fixed else (r["P"], r["P_err"])
        prim, sec, ph = x.prim, x.sec, x.phase_sec
        if fixed:  # near-equal eclipses: the deeper alternate eclipse is the primary
            e = json.load(open(f"{W}/vet/EB_TIC{int(x.tic)}/eb.json"))["alt_depths"]
            prim, sec, ph = max(e), min(e), 0.5
        if int(x.tic) in PHASE_FIX: sec, ph = PHASE_FIX[int(x.tic)]
        t0 = r["T0_BTJD"]
        if sec == sec and sec > prim:  # the eclipse at phase ph is deeper: make it the primary
            t0 = t0 + ph * P; prim, sec, ph = sec, prim, (1 - ph) % 1
        name = r["simbad"][0].split(" (")[0] if r["simbad"] else ""
        num = lambda k: round(float(g[k]), 3) if g.get(k) not in (None, "None", "nan") and g.get(k) == g.get(k) else None
        detached = x.kind == "detached"
        out.append(dict(tic=int(x.tic), name=name, ra=round(r["ra"], 5), dec=round(r["dec"], 5), Tmag=round(r["Tmag"], 2),
                        G=num("phot_g_mean_mag"), plx=num("parallax"), ruwe=num("ruwe"), Teff_TIC=r["Teff_TIC"],
                        P=P, P_err=Perr, T0_BTJD=t0, primary_pct=round(prim * 100, 2), secondary_pct=round(sec * 100, 3),
                        phase_sec=round(ph, 3),
                        ecosw=round(np.pi / 2 * (ph - 0.5), 3) if detached and not fixed and x.sec_sig > 5 else None,
                        kind=x.kind, n_sectors=x.nsec, worst_offset_sigma=round(x.worst_loc, 1), kostov_unvetted=r["kostov_unvetted"],
                        rv_amp=g.get("rv_amplitude_robust"), nss=g.get("non_single_star"), note=NOTES.get(int(x.tic), ""),
                        neighbours_able_to_mimic="; ".join(f"{n['sep']}\" T~{n['Test']}" for n in nb.get(str(int(x.tic)), []))
                        if isinstance(nb.get(str(int(x.tic))), list) else nb.get(str(int(x.tic)), "")))
    new = pd.DataFrame(out)
    return old, new


if __name__ == "__main__":
    old, new = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    dest = sys.argv[3] if len(sys.argv) > 3 else f"{W}/vet/new_ebs_batch.csv"
    new.to_csv(dest, index=False)
    print(len(new), "rows ->", dest)
