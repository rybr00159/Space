"""Triage EB-track candidates: broad novelty (SIMBAD types + all VizieR catalogues), then full EB analysis."""
import os, sys, json, subprocess
import pandas as pd
import known

W = "/home/user/work"
VAR_OTYPES = ("EB*", "SB*", "V*", "Ro*", "RR*", "dS*", "BY*", "RS*", "Pu*", "LP*", "Mi*", "sr*", "Ce*", "WU*", "Al*", "bL*", "El*",
              "**", "ev*", "gD*", "Ir*", "Er*", "Or*", "TT*", "Y*O", "a2*", "PM*?")
VARIAB_SKIP = ("IV/38", "IV/39", "J/AJ/155/39/table5", "J/AJ/155/39/table6", "I/355/epphot", "I/358/varisum")


def novelty(r):
    try:
        viz, ncat = known.vizier_all(r["ra"], r["dec"], 15)
        sim = known.simbad(r["ra"], r["dec"], 20)
    except Exception as e:
        return dict(tic=int(r["tic"]), nov_error=str(e)[:120])
    viz = [v for v in viz if not v.startswith(VARIAB_SKIP)]
    ot = [f"{x['main_id']}:{x['otype']}" for x in sim if float(x["sep"]) < 10]
    flagged = any(o.split(":")[-1] in VAR_OTYPES for o in ot)
    return dict(tic=int(r["tic"]), n_catalogs=ncat, viz_variability=";".join(viz), simbad=";".join(ot),
                novel=(not viz) and (not flagged))


if __name__ == "__main__":
    for sector in sys.argv[1:]:
        s = "multi" if sector == "multi" else f"s{int(sector):04d}"
        f = f"{W}/results/new_{s}.csv"
        if not os.path.exists(f): continue
        c = pd.read_csv(f); c = c[(c.track == "eb") & ~c.viz_known]
        cache = f"{W}/results/eb_novelty_{s}.csv"
        done = pd.read_csv(cache) if os.path.exists(cache) else pd.DataFrame(columns=["tic"])
        todo = c[~c.tic.isin(done.tic)]
        rows = [novelty(r) for r in todo.to_dict("records")]
        if rows:
            done = pd.concat([done, pd.DataFrame(rows)], ignore_index=True); done.to_csv(cache, index=False)
        nov = done[done.novel == True]
        print(f"{s}: {len(c)} uncatalogued-by-period EB candidates; {len(nov)} pass broad novelty")
        for t in nov.tic:
            if not os.path.exists(f"{W}/vet/EB_TIC{int(t)}/eb.json"):
                subprocess.run([f"{W}/venv/bin/python", "eb_analysis.py", str(int(t)), sector], cwd=f"{W}/src",
                               capture_output=True, timeout=1800)
