"""Validation of the dF-map PRF localisation with injections into the real pixel data.

Per sector, fake transits are injected at the real transit epochs shifted by +P/4 (windows do not overlap the real
transits; same template shape; total flux drop 13 e-/s ~ the observed value) with three spatial profiles:
  inj_O   : dip proportional to the real out-of-transit image O (i.e. the real PSF of the blend; tests PRF-model
            mismatch with the real PSF).  Expected position = light centre of O, measured by a single-PRF fit to O.
  inj_T   : dip with the PRF model at the target position (code / coordinate-convention check -> expect (0, 0)).
  inj_N   : dip with the PRF model at the 9" neighbour (expect the neighbour position).
  inj_Treal: dip proportional to O minus the PRF model of all stars except the target (an empirical estimate of the
            target's own PSF), expected at the target.
The recovered positions use exactly the pipeline of prf_localize.py (dF map, single-PRF fit, 9x9 region).
Writes localize_validate.json.
"""
import os, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import prf_localize as pl
from prf_localize import gaia, TARGET_ID, NB_ID

HERE = os.path.dirname(os.path.abspath(__file__))
tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]; nbr = gaia[gaia.source_id == NB_ID].iloc[0]
DF_INJ = 13.0


def fit_pos(M, img, err, x0, y0, H=pl.HALF, half=4):
    sl = slice(H - half, H + half + 1)
    wgt = 1 / err[sl, sl].ravel() ** 2
    def chi(p):
        m = M.img(*p)[sl, sl].ravel()
        a = np.sum(m * img[sl, sl].ravel() * wgt) / np.sum(m * m * wgt)
        return np.sum((img[sl, sl].ravel() - a * m) ** 2 * wgt)
    return minimize(chi, [x0, y0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-5)).x


def main():
    out = {}
    for fpath in sorted(glob.glob("/home/user/work/tic87/data/tesscut/*.fits")):
        S = pl.load_sector(fpath); M = pl.PRFModel(S)
        xt, yt = map(float, S["w"].world_to_pixel_values(tgt.ra, tgt.dec)); xn, yn = map(float, S["w"].world_to_pixel_values(nbr.ra, nbr.dec))
        dF, edF, O = pl.dF_map(S, S["blocks"])
        # full-frame O (median of all cadences) for injection profiles
        Ofull = np.median(S["flux"], axis=0)
        # single-PRF fit to O (light centre), uniform weights
        xo, yo = fit_pos(M, O, np.ones_like(O), xt, yt)
        r = dict(O_light_centre_EN=list(pl.to_EN(S["w"], xo, yo, xt, yt)), neighbour_EN=list(pl.to_EN(S["w"], xn, yn, xt, yt)))
        xr, yr = fit_pos(M, dF, edF, xt, yt)
        r["real_dF_EN"] = list(pl.to_EN(S["w"], xr, yr, xt, yt))
        # expected O light centre from the Gaia + PRF model
        model = np.zeros(S["shape"])
        for row in M.stars.itertuples():
            model += M.prf.locate(float(row.px), float(row.py), S["shape"]) * 10 ** (-0.4 * (row.Tmag - 10))
        sub = model[S["ys"], S["xs"]]
        xm, ym = fit_pos(M, sub, np.ones_like(sub), xt, yt)
        r["model_light_centre_EN"] = list(pl.to_EN(S["w"], xm, ym, xt, yt))
        # empirical target PSF: O minus PRF model of the other stars (fluxes from Gaia, scaled to O's total)
        others = np.zeros(S["shape"])
        for row in M.stars.itertuples():
            if row.source_id == TARGET_ID: continue
            others += M.prf.locate(float(row.px), float(row.py), S["shape"]) * 10 ** (-0.4 * (row.Tmag - 10))
        scale = np.sum(Ofull) / np.sum(model)
        Treal = Ofull - scale * others
        profiles = dict(inj_O=Ofull / Ofull[S["ys"], S["xs"]].sum(),
                        inj_T=M.prf.locate(xt, yt, S["shape"]), inj_N=M.prf.locate(xn, yn, S["shape"]),
                        inj_Treal=Treal / Treal[S["ys"], S["xs"]].sum())
        for k in ("inj_T", "inj_N"):
            profiles[k] = profiles[k] / profiles[k][S["ys"], S["xs"]].sum()
        # fake epochs: shift by P/4
        t = S["t"]; E = S["E"]
        fake_blocks = []
        for e, ii, Tc in S["blocks"]:
            Tf = Tc + pl.P / 4
            m = np.where(np.abs(t - Tf) < pl.W)[0]
            if len(m) < 20: continue
            inn = np.abs(t[m] - Tf) < pl.tmpl.T14 / 2
            pre = np.sum(t[m] < Tf - pl.tmpl.T14 / 2); post = np.sum(t[m] > Tf + pl.tmpl.T14 / 2)
            if inn.sum() * E / 86400 / pl.tmpl.T14 >= 0.7 and pre > 5 and post > 5:
                fake_blocks.append((e, m, Tf))
        flux0 = S["flux"].copy()
        for k, prof in profiles.items():
            dip = np.zeros(len(t))
            for e, m, Tf in fake_blocks:
                dip[m] = (1 - pl.tmpl(t[m] - Tf, E)) / pl.tmpl.depth * DF_INJ
            S["flux"] = flux0 - dip[:, None, None] * prof[None, :, :]
            d2, e2, _ = pl.dF_map(S, fake_blocks)
            x2, y2 = fit_pos(M, d2, e2, xt, yt)
            # also the dF map without injection at the fake epochs (null)
            r[k + "_EN"] = list(pl.to_EN(S["w"], x2, y2, xt, yt))
            r[k + "_dF_sum"] = float(d2.sum())
        S["flux"] = flux0
        d0, e0, _ = pl.dF_map(S, fake_blocks)
        r["null_dF_sum_9x9"] = float(d0[1:-1, 1:-1].sum()); r["n_fake"] = len(fake_blocks)
        out[S["s"]] = r
        print(S["s"], {k: (np.round(v, 2).tolist() if isinstance(v, list) else round(v, 2)) for k, v in r.items()}, flush=True)
    json.dump(out, open(os.path.join(HERE, "localize_validate.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
