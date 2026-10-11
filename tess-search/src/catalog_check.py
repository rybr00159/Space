"""VizieR variable-star / EB catalogue check for candidates (period-matched, within 150 arcsec)."""
import json
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
import known


def check_one(r, radius=150.0):
    try:
        hits = known.vizier_cone(float(r["ra"]), float(r["dec"]), radius)
    except Exception as e:
        return dict(tic=int(r["tic"]), viz_error=str(e)[:100], viz_match="", viz_near="")
    match, near = [], []
    for h in hits:
        tag = f"{h['cat']}|{h['id']}|{h['type']}|P={h['period']}|{h['sep']:.0f}\""
        try:
            per = float(h["period"]) if h["period"] not in (None, "") else np.nan
        except (TypeError, ValueError):
            per = np.nan
        eb_cat = any(k in str(h["cat"]) for k in ("EB", "Kostov", "Rowan", "IJspeert", "Prsa"))
        if np.isfinite(per) and known.period_match(float(r["P"]), per, tol=0.005, kmax=4):
            match.append(tag)
        elif eb_cat and h["sep"] <= 21:
            match.append(tag + "(EB catalogue, within 1 px)")
        elif h["sep"] <= 25:
            near.append(tag)
    return dict(tic=int(r["tic"]), viz_match=";".join(match), viz_near=";".join(near))


def check(df, workers=4):
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(check_one, df.to_dict("records")))
    return df.merge(pd.DataFrame(res), on="tic", how="left")
