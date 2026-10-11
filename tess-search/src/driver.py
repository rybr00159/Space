"""Incremental downstream pipeline: select -> classify -> catalogue check -> stage 2, per sector."""
import os, sys, json, subprocess, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from classify import classify

W = "/home/user/work"; PY = f"{W}/venv/bin/python"


def step(sector):
    s = f"s{sector:04d}"
    subprocess.run([PY, "select_cands.py", str(sector), f"{W}/results/{s}", f"{W}/data/{s}", f"{W}/results/cands_{s}.csv",
                    f"{W}/meta/badwin_{s}.json"], cwd=f"{W}/src", check=True, capture_output=True)
    c = classify(pd.read_csv(f"{W}/results/cands_{s}.csv"), pd.read_parquet(f"{W}/meta/tic_{s}.parquet"))
    c.to_csv(f"{W}/results/cands_{s}_classified.csv", index=False)
    new = c[~c.is_known & c.track.isin(["planet", "eb"])]
    from catalog_check import check
    prev = f"{W}/results/catcheck_{s}.csv"
    done = pd.read_csv(prev) if os.path.exists(prev) else pd.DataFrame(columns=["tic"])
    todo = new[~new.tic.isin(done.tic)]
    if len(todo):
        cc = check(todo)[["tic", "viz_match", "viz_near"]]
        done = pd.concat([done, cc], ignore_index=True); done.to_csv(prev, index=False)
    new = new.merge(done[["tic", "viz_match", "viz_near"]], on="tic", how="left")
    new["viz_known"] = new.viz_match.fillna("").str.len() > 0
    new.to_csv(f"{W}/results/new_{s}.csv", index=False)
    return c, new[~new.viz_known]


if __name__ == "__main__":
    sector = int(sys.argv[1])
    c, new = step(sector)
    print(f"S{sector}: {len(c)} candidates; tracks {c.track.value_counts().to_dict()}; new planet/eb not in VizieR: {len(new)} ({new.track.value_counts().to_dict()})")
    inp = f"{W}/results/stage2_in_s{sector:04d}.csv"
    new[new.track == "planet"].to_csv(inp, index=False)
    subprocess.run([PY, "stage2.py", inp, f"{W}/results/stage2_s{sector:04d}.jsonl", "3"], cwd=f"{W}/src")
