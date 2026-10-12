"""Limb-darkening / flat-bottom test: does the transit bottom show the curvature expected from stellar limb
darkening (planet transit) or is it flatter (e.g. a total occultation of a self-luminous body)?
Maximum-likelihood batman fits (same data / fixed per-sector sigma / profiled offsets as vshape.py):
  LDprior : Claret prior (as the MCMC)
  LDfree  : u1, u2 free within the physical region, no prior
  LD0     : u1 = u2 = 0 (uniform disc: flat bottom, no curvature)
  LDtic / LDgsp : u1, u2 fixed at the Claret values for the TIC (7276 K, logg 4.2) / GSP-Phot (5696 K, logg 3.5)
Writes ld_test.json."""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import tmodel as tm
import vshape as vs

HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, "ld_prior.json")))


def chi2_fixed_ld(x5, u1, u2):
    return vs.chi2_batman(np.r_[x5, u1, u2], "free") - (((u1 - tm.U1_MU) / tm.U1_SD) ** 2 + ((u2 - tm.U2_MU) / tm.U2_SD) ** 2)


def chi2_free_ld(x7):
    u1, u2 = x7[5], x7[6]
    c = vs.chi2_batman(x7, "free")
    if c >= 1e29: return c
    return c - (((u1 - tm.U1_MU) / tm.U1_SD) ** 2 + ((u2 - tm.U2_MU) / tm.U2_SD) ** 2)


xb = np.array([vs.mc["best"][n] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
out = {}
r = vs.nm(lambda x: vs.chi2_batman(x, "free"), xb, [1e-4, 1e-7, 1e-3, 0.02, 0.005, 0.02, 0.02])
out["LDprior"] = dict(chi2=float(vs.chi2_batman(r.x, "free") - (((r.x[5] - tm.U1_MU) / tm.U1_SD) ** 2 + ((r.x[6] - tm.U2_MU) / tm.U2_SD) ** 2)), x=r.x.tolist())
r = vs.nm(chi2_free_ld, xb, [1e-4, 1e-7, 1e-3, 0.02, 0.005, 0.05, 0.05], n=6)
out["LDfree"] = dict(chi2=float(r.fun), x=r.x.tolist())
for lab, (u1, u2) in (("LD0", (1e-6, 0.0)), ("LDtic", (LD["TIC_7276_logg4.2"]["u1"], LD["TIC_7276_logg4.2"]["u2"])),
                      ("LDgsp", (LD["GSPPhot_5696_logg3.5"]["u1"], LD["GSPPhot_5696_logg3.5"]["u2"]))):
    if lab == "LD0":
        f = lambda x5: vs.chi2_from_models(vs.eval_models(vs.integ(vs.g.profile(x5[1], x5[2], tm.a_from_rho(x5[4], x5[1]), x5[3], 0.0, 0.0)), x5[0], x5[1])) \
            if (0.01 < x5[2] < 0.6 and 0 <= x5[3] < 1 + x5[2] and 0.002 < x5[4] < 10 and tm.a_from_rho(x5[4], x5[1]) > 1.5) else 1e30
    else:
        f = lambda x5, u1=u1, u2=u2: chi2_fixed_ld(x5, u1, u2)
    r = vs.nm(f, xb[:5], [1e-4, 1e-7, 1e-3, 0.02, 0.005], n=4)
    out[lab] = dict(chi2=float(r.fun), x=r.x.tolist(), u=[u1, u2])
base = out["LDprior"]["chi2"]
for k in out:
    out[k]["dchi2_vs_LDprior"] = out[k]["chi2"] - base
    x = out[k]["x"]
    out[k]["k"], out[k]["b"], out[k]["rho"] = x[2], x[3], x[4]
json.dump(out, open(os.path.join(HERE, "ld_test.json"), "w"), indent=1)
print(json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "x"} for k, v in out.items()}, indent=1))
print("LDfree u1,u2:", out["LDfree"]["x"][5:7])
