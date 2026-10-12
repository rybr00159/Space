"""Light curve for TRICERATOPS with TIC 87055315 as the target (runner: /home/user/work/src/run_tri.py, unmodified).

Choice: the runner calls target.calc_depths(tdepth) and target.calc_probs(...) with TRICERATOPS' default
dilution_corrected=True, i.e. flux_0 must be normalised so that the out-of-transit flux of the TARGET is 1 (on-target,
dilution-corrected light curve, like PDCSAP) and tdepth must be the on-target depth; the dip in the aperture is then
computed internally as tdepth x fluxratio_target and each nearby star's light curve as
1 + (fluxratio_target/fluxratio_i)(flux_0 - 1).  Feeding the blended (SAP / pair) light curve to this runner would
make TRICERATOPS treat a 3.4-ppt dip as the neighbour's own depth (6x too shallow).  We therefore give the
neighbour's undiluted light curve: the PSF pair light curve (F_t + F_n; independent of the registration-degenerate
flux split) rescaled with D = F_t/F_n = 1/r:  f_nb = 1 + (f_pair - 1)(1 + D).
Sectors 66, 93, 104 only (200-s cadence): TRICERATOPS models a 2-min exposure, so the 30-min (S13) and 10-min (S39)
data would smear ingress/egress and bias the shape towards V-shaped (EB) scenarios.
TRI_DEPTH = undiluted model depth from mcmc/mcmc_pair_Dfree.json (the runner's automatic estimate, 1 - min of the
3-bin running mean of the 5-min binned fold, is biased deep by noise for this light curve).
Output: data/tri_lc_nb_200s.npz (t, f, sigma), data/tri_inputs.json
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbtools as nt
ROOT = nt.ROOT
st = json.load(open(os.path.join(ROOT, "results/stellar.json")))["adopted"]
mc = json.load(open(os.path.join(ROOT, "mcmc/mcmc_pair_Dfree.json")))
p = mc["params"]
D = st["D_pair"]
lc = nt.load_lc("pair")
T, F = [], []
for s in (66, 93, 104):
    v = lc[s]
    T.append(v["t"]); F.append(1 + (v["f"] - 1) * (1 + D))
t = np.concatenate(T); f = np.concatenate(F)
P, T0 = p["P"]["med"], p["T0"]["med"]
ph = ((t - T0) / P + 0.5) % 1 - 0.5
oot = np.abs(ph * P) > 0.15
sig = float(1.4826 * np.median(np.abs(f[oot] - np.median(f[oot]))))
np.savez(os.path.join(ROOT, "data/tri_lc_nb_200s.npz"), t=t, f=f, sigma=sig)
# runner folds with t0 and P; give a T0 near the middle of the 200-s data
n_mid = np.round((np.median(t) - T0) / P)
inp = dict(P=P, T0_ref=T0 + n_mid * P, dur=p["T14_h"]["med"] / 24, depth_undiluted=p["depth_undiluted_ppm"]["med"] * 1e-6,
           depth_undiluted_err=0.5 * (p["depth_undiluted_ppm"]["lo"] + p["depth_undiluted_ppm"]["hi"]) * 1e-6,
           sigma_point=sig, n_points=int(len(t)), D=D, sectors=[66, 93, 104],
           tri_mass_adopted=st["M"]["med"], tri_rad_adopted=st["R"]["med"],
           tri_plx_adopted=json.load(open(os.path.join(ROOT, "results/stellar.json")))["distance"]["parallax_corr"])
json.dump(inp, open(os.path.join(ROOT, "data/tri_inputs.json"), "w"), indent=1)
print(json.dumps(inp, indent=1))
