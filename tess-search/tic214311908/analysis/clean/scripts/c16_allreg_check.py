"""c16: is the disappearance of the late-early depth contrast in the 'allreg' variant (variables' aperture light curves
as explicit pixel regressors) significant? Paired per-transit differences (allreg - all, raw flux) at the real epochs
compared with the same differences at the 60 fake epochs per sector (c11 table). Output: out/allreg_check.json"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import CL, SECTORS, jd
import numpy as np, pandas as pd
D = pd.read_csv(f"{CL}/out/modes_per_transit.csv")
def paired(kind):
    a = D[(D.version == "inv_allreg_raw") & (D.kind == kind)].set_index(["sector", "tc"]).depth
    b = D[(D.version == "inv_all_raw") & (D.kind == kind)].set_index(["sector", "tc"]).depth
    j = pd.concat([a.rename("reg"), b.rename("all")], axis=1, join="inner").reset_index()
    j["diff"] = j.reg - j["all"]; return j
R = paired("real"); F = paired("fake")
sd = {s: float(F[F.sector == s]["diff"].std()) for s in SECTORS}
R["err"] = R.sector.map(sd)
def wm(m):
    x = R[m]; w = 1 / x.err ** 2; return float((w * x["diff"]).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
a = wm(R.sector.isin([39, 66])); b = wm(R.sector.isin([103, 104]))
out = dict(null_sd_by_sector=sd, real=R[["sector", "tc", "reg", "all", "diff"]].round(3).values.tolist(),
           S39S66_mean_diff=a, S103S104_mean_diff=b, contrast_change=[b[0] - a[0], float(np.hypot(a[1], b[1]))],
           chi2_real_diffs=float(((R["diff"] / R.err) ** 2).sum()), n=len(R))
z = (b[0] - a[0]) / np.hypot(a[1], b[1])
print(f"null sd of paired diff per sector {sd}")
print(f"S39+66 mean diff {a[0]:+.2f}+-{a[1]:.2f}; S103+104 {b[0]:+.2f}+-{b[1]:.2f}; contrast change {b[0]-a[0]:+.2f}+-{np.hypot(a[1],b[1]):.2f} (z {z:.2f}); "
      f"chi2 of real paired diffs {out['chi2_real_diffs']:.1f}/{len(R)}")
jd(out, f"{CL}/out/allreg_check.json")
