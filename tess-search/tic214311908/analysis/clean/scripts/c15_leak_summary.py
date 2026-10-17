"""c15: summary table of the field-variable leak before/after cleaning (from out/leak_inv.json, out/leak_ref.json,
out/tile_transfer.json). REF tiled amplitudes are divided by the mean tiled transfer (c12) to give the equivalent
amplitude in a contiguous light curve. Output: out/leak_summary.json"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import CL, SECTORS, jd
import json, numpy as np
LI = json.load(open(f"{CL}/out/leak_inv.json")); LR = json.load(open(f"{CL}/out/leak_ref.json"))
tr = float(np.mean([v["transfer"] for v in json.load(open(f"{CL}/out/tile_transfer.json")).values()]))
out = dict(tiled_transfer=tr, inv={}, ref={})
for s in SECTORS:
    s = str(s)
    out["inv"][s] = {m: dict(RR_amp1_ppt=LI[s][m]["M3V_B"]["RR"]["amp1"], RR_p=LI[s][m]["M3V_B"]["RR"]["p_amp"], RR_UL95_ppt=LI[s][m]["M3V_B"]["RR"]["ul95_amp"],
                         EB075_dchi2=LI[s][m]["M3V_B"]["EB075"]["dchi2"], EB075_p=LI[s][m]["M3V_B"]["EB075"]["p_dchi2"],
                         EB159_dchi2=LI[s][m]["M3V_B"]["EB159"]["dchi2"], EB159_p=LI[s][m]["M3V_B"]["EB159"]["p_dchi2"],
                         dx_RR_amp_px=LI[s][m]["dx"]["RR"]["amp1"], dx_RR_p=LI[s][m]["dx"]["RR"]["p_amp"], rstd_ppt=LI[s][m]["M3V_B"]["rstd_ppt"],
                         modelA_RR_amp1_ppt=LI[s][m]["M3V_A"]["RR"]["amp1"], modelA_RR_p=LI[s][m]["M3V_A"]["RR"]["p_amp"])
                     for m in ("none", "rr3", "all", "allreg")}
    out["ref"][s] = {m: dict(RR_amp1_tiled_ppt=LR[s][m]["RR"]["amp1"], RR_p=LR[s][m]["RR"]["p_amp"],
                         RR_amp1_equiv_ppt=LR[s][m]["RR"]["amp1"] / tr, RR_UL95_equiv_ppt=LR[s][m]["RR"]["ul95_amp"] / tr,
                         EB075_p=LR[s][m]["EB075"]["p_dchi2"], EB159_p=LR[s][m]["EB159"]["p_dchi2"], rstd_ppt=LR[s][m]["rstd_ppt"])
                     for m in ("orig", "all")}
    i0, i1 = out["inv"][s]["none"], out["inv"][s]["all"]; r0, r1 = out["ref"][s]["orig"], out["ref"][s]["all"]
    print(f"S{s}: INV RR {i0['RR_amp1_ppt']:.2f} (p {i0['RR_p']:.3f}) -> {i1['RR_amp1_ppt']:.2f} (p {i1['RR_p']:.3f}, UL95 {i1['RR_UL95_ppt']:.2f}); "
          f"dx RR p {i0['dx_RR_p']:.2f} -> {i1['dx_RR_p']:.2f}; allreg {out['inv'][s]['allreg']['RR_amp1_ppt']:.2f} (p {out['inv'][s]['allreg']['RR_p']:.3f}); "
          f"REF equiv {r0['RR_amp1_equiv_ppt']:.2f} (p {r0['RR_p']:.2f}) -> {r1['RR_amp1_equiv_ppt']:.2f} (p {r1['RR_p']:.2f}, UL95 {r1['RR_UL95_equiv_ppt']:.2f}); "
          f"INV EB075 p {i0['EB075_p']:.2f}->{i1['EB075_p']:.2f}, EB159 p {i0['EB159_p']:.2f}->{i1['EB159_p']:.2f}")
jd(out, f"{CL}/out/leak_summary.json")
