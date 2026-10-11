"""Downstream pipeline for the combined S104+S105 search: rescore -> known -> classify -> catalogue -> stage 2."""
import os, sys, glob, json, subprocess, warnings
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
from bls import detrend, clean, diagnostics
from run_search_multi import index, get
from classify import classify
from catalog_check import check
import known

W = "/home/user/work"; PY = f"{W}/venv/bin/python"


def main():
    res = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{W}/results/multi/group_*.parquet"))], ignore_index=True)
    if "error" in res: res = res[res.error.isna()]
    c = res[(res.snr_red >= 8.5) & res.passed].copy()
    WA = {int(k): v["windows"] for k, v in json.load(open(f"{W}/meta/badwin_s0104.json")).items()}
    WB = {int(k): v["windows"] for k, v in json.load(open(f"{W}/meta/badwin_s0105.json")).items()}
    iA, iB = index(f"{W}/data/s0104"), index(f"{W}/data/s0105")
    cache = {}
    def z_of(p):
        if p not in cache: cache[p] = {k: v for k, v in np.load(p).items()}
        return cache[p]
    rows = []
    for r in c.itertuples():
        ta, fa = get(z_of(iA[r.tic][0]), iA[r.tic][1], WA); tb, fb = get(z_of(iB[r.tic][0]), iB[r.tic][1], WB)
        t, f = np.concatenate([ta, tb]), np.concatenate([fa, fb]); t, f, s = clean(t, f)
        d = diagnostics(t, f, r.P, r.t0, r.dur, s)
        rows.append(dict(tic=r.tic, **{f"m_{k}": v for k, v in d.items()}))
    c = c.merge(pd.DataFrame(rows), on="tic")
    c = c[(c.m_snr_red >= 9) & (c.m_n_tr >= 3) & (c.m_max_event_frac <= 0.5) & (c.m_depth > 0)].copy()
    ks = c.apply(lambda r: known.known_signals(int(r.tic), r.P), axis=1)
    c["known_on_star"] = ks.map(lambda x: ";".join(x[0][:6])); c["known_match"] = ks.map(lambda x: ";".join(x[1][:6]))
    c["is_known"] = c.known_match.str.len() > 0; c["sector"] = "104,105"
    tic = pd.concat([pd.read_parquet(f"{W}/meta/tic_s0104.parquet"), pd.read_parquet(f"{W}/meta/tic_s0105.parquet")]).drop_duplicates("ID")
    c = classify(c, tic)
    c.to_csv(f"{W}/results/cands_multi_classified.csv", index=False)
    new = c[~c.is_known & c.track.isin(["planet", "eb"])]
    cc = check(new)[["tic", "viz_match", "viz_near"]]
    new = new.merge(cc, on="tic", how="left"); new["viz_known"] = new.viz_match.fillna("").str.len() > 0
    new.to_csv(f"{W}/results/new_multi.csv", index=False)
    print(f"multi: {len(res)} searched, {len(c)} candidates {c.track.value_counts().to_dict()}, new not in VizieR: "
          f"{(~new.viz_known).sum()} ({new[~new.viz_known].track.value_counts().to_dict()})", flush=True)
    inp = f"{W}/results/stage2_in_multi.csv"
    new[~new.viz_known & (new.track == "planet")].to_csv(inp, index=False)
    subprocess.run([PY, "stage2.py", inp, f"{W}/results/stage2_multi.jsonl", "3"], cwd=f"{W}/src")


if __name__ == "__main__":
    main()
