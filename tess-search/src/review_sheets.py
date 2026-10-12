"""Visual-review sheets: each candidate folded at P (full phase) and at 2P, plus a zoom on the primary."""
import sys, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import stage2
from eb_analysis import fold_profile

sel = pd.read_csv(sys.argv[1]); tics = [int(x) for x in open(sys.argv[2]).read().split()]
out = sys.argv[3]
sel = sel.set_index("tic").loc[tics].reset_index()
per = 8
for page in range(0, len(sel), per):
    sub = sel.iloc[page:page + per]
    fig, axes = plt.subplots(len(sub), 3, figsize=(15, 2.3 * len(sub)), squeeze=False)
    for ax, r in zip(axes, sub.itertuples()):
        secs = stage2.load_all(int(r.tic))
        t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
        for a, P, lab in ((ax[0], r.P, "P"), (ax[1], 2 * r.P, "2P")):
            ph = ((t - r.t0) / P + 0.25) % 1.0 - 0.25
            pb, prof, _ = fold_profile(t, f, P, r.t0 - 0.25 * P, nb=400)
            a.plot(ph, f, ",", color="0.6"); a.plot(pb - 0.25, prof, "-", lw=0.8, color="C0")
            a.set_title(f"TIC {r.tic}  {lab} = {P:.5f} d  sec@{r.phase_sec:.3f}", fontsize=8); a.tick_params(labelsize=7)
        x = (((t - r.t0) / r.P + 0.5) % 1.0 - 0.5) * r.P * 24
        m = np.abs(x) < max(3 * r.dur_h, 2)
        ax[2].plot(x[m], f[m], ".", ms=1, color="0.6"); ax[2].set_title(f"primary zoom (h), {len(secs)} sectors", fontsize=8)
        ax[2].tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(f"{out}_{page // per + 1}.png", dpi=75); plt.close(fig)
print("done")
