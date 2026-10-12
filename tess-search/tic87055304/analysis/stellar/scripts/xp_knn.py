"""Empirical, parallax-free absolute magnitude of TIC 87055304 from its Gaia DR3 XP spectrum (shape only):
k-nearest-neighbour regression of M_G' = G + 5 log10(plx[mas]/100) (not extinction corrected; the dust along this
sightline is foreground to ~300 pc, see reddening_distance.py) on XP coefficients of 3745 stars within 1.7 deg with
parallax_over_error > 10, RUWE < 1.4, 10.5 < G < 14, 0.60 < BP-RP < 1.15 (xp_fetch.py).
Features: first NB BP and NR RP continuous coefficients divided by phot_g_mean_flux (removes the distance scale,
keeps the spectral shape and colour), standardised, PCA-reduced. Validation: 10-fold cross-validation (residual
scatter, bias, and dwarf/subgiant confusion); gradient-boosting regression as a second method.
Also a class likelihood ratio: dwarf (M_G' > 3.6) vs subgiant (1.4 < M_G' < 3.0) from the class-normalised
neighbour counts among the K nearest training stars.
Output: out/xp/knn_result.json, fig/xp_knn.png
"""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsRegressor, NearestNeighbors
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = "/home/user/work/tic87/stellar"; OUT = os.path.join(ROOT, "out/xp")
TID = 6720986204578421888
meta = pd.read_csv(os.path.join(OUT, "train_meta.csv"))
z = np.load(os.path.join(OUT, "xp_coeffs.npz"))
co = pd.DataFrame({"source_id": z["sid"]}); co["i"] = np.arange(len(co))
m = meta.merge(co, on="source_id")
print("rows with XP", len(m))
NB, NR, NPC, K = 25, 25, 20, 15
X = np.hstack([z["bp"][m.i, :NB], z["rp"][m.i, :NR]]) / m.phot_g_mean_flux.values[:, None]
m["MG"] = m.phot_g_mean_mag + 5 * np.log10(m.parallax / 100.)
it = np.where(m.source_id == TID)[0][0]
tr = np.where(m.source_id != TID)[0]
sc = StandardScaler().fit(X[tr]); Xs = sc.transform(X)
pca = PCA(n_components=NPC).fit(Xs[tr]); Xp = pca.transform(Xs)
y = m.MG.values
# 10-fold CV
kf = KFold(10, shuffle=True, random_state=1)
pred_knn = np.full(len(m), np.nan); pred_gb = np.full(len(m), np.nan)
for a, b in kf.split(tr):
    A, B = tr[a], tr[b]
    knn = KNeighborsRegressor(K, weights="distance").fit(Xp[A], y[A]); pred_knn[B] = knn.predict(Xp[B])
    gb = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=1).fit(Xp[A], y[A]); pred_gb[B] = gb.predict(Xp[B])
knn = KNeighborsRegressor(K, weights="distance").fit(Xp[tr], y[tr])
gb = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=1).fit(Xp[tr], y[tr])
t_knn = float(knn.predict(Xp[[it]])[0]); t_gb = float(gb.predict(Xp[[it]])[0])
# residual scatter for stars similar to the target (BP-RP within 0.06) and in the predicted-M_G neighbourhood
sim = tr[np.abs(m.bp_rp.values[tr] - m.bp_rp.values[it]) < 0.06]
res_k = y[sim] - pred_knn[sim]; res_g = y[sim] - pred_gb[sim]
def rob(x): return float(1.4826 * np.median(np.abs(x - np.median(x))))
near = sim[np.abs(pred_knn[sim] - t_knn) < 0.5]
dw = y > 3.6; sg = (y > 1.4) & (y < 3.0)
conf = dict(subgiants_predicted_as_dwarf=float(np.mean(pred_knn[tr][sg[tr]] > 3.6)), n_sub=int(sg[tr].sum()),
            dwarfs_predicted_as_subgiant=float(np.mean(pred_knn[tr][dw[tr]] < 3.0)), n_dwarf=int(dw[tr].sum()))
# class likelihood ratio from KK nearest neighbours
out_lr = {}
for KK in (15, 30, 60):
    nn = NearestNeighbors(n_neighbors=KK).fit(Xp[tr]); dist, ind = nn.kneighbors(Xp[[it]])
    nb = tr[ind[0]]
    nd, ns = int(dw[nb].sum()), int(sg[nb].sum())
    Nd, Ns = int(dw[tr].sum()), int(sg[tr].sum())
    lr = ((nd + 0.5) / Nd) / ((ns + 0.5) / Ns)
    out_lr[str(KK)] = dict(n_dwarf=nd, n_sub=ns, N_dwarf=Nd, N_sub=Ns, likelihood_ratio_dwarf_vs_subgiant=float(lr),
                           MG_neighbours_med=float(np.median(y[nb])), MG_neighbours_p16=float(np.percentile(y[nb], 16)),
                           MG_neighbours_p84=float(np.percentile(y[nb], 84)), median_nn_dist=float(np.median(dist)))
# distance from predicted M_G' (extinction is in M_G' implicitly, so d = 10^((G - M_G' + 5)/5))
G = float(m.phot_g_mean_mag.values[it])
sig = float(np.hypot(rob(y[near] - pred_knn[near]) if len(near) > 10 else rob(res_k), 0.0))
def d_of(MG): return 10 ** ((G - MG + 5) / 5)
res = dict(n_train=int(len(tr)), K=K, NPC=NPC, NB=NB, NR=NR, target_MG_knn=t_knn, target_MG_gb=t_gb,
           cv_knn_scatter_similar=rob(res_k), cv_knn_bias_similar=float(np.median(res_k)),
           cv_gb_scatter_similar=rob(res_g), cv_gb_bias_similar=float(np.median(res_g)),
           cv_knn_scatter_near_pred=rob(y[near] - pred_knn[near]) if len(near) > 10 else None, n_near=int(len(near)),
           cv_knn_frac_outliers_gt1mag=float(np.mean(np.abs(res_k) > 1.0)),
           confusion=conf, class_LR=out_lr, sigma_MG_adopted=sig,
           MG_if_dr3_plx=float(G + 5 * np.log10(2.832 / 100)), MG_gspphot_marcs=2.0406, MG_gspphot_phoenix=4.2757,
           d_med=float(d_of(t_knn)), d_p16=float(d_of(t_knn + sig)), d_p84=float(d_of(t_knn - sig)),
           d_gb=float(d_of(t_gb)))
# how many sigma is the GSP-Phot MARCS M_G (2.04) from the XP prediction?
res["MARCS_offset_sigma"] = float((t_knn - 2.0406) / sig)
res["plx_implied_mas"] = float(1000 / res["d_med"])
json.dump(res, open(os.path.join(OUT, "knn_result.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].scatter(y[tr], pred_knn[tr], s=2, c="0.5"); ax[0].plot([0, 7], [0, 7], "k-")
ax[0].axhline(t_knn, color="r", label=f"target pred {t_knn:.2f}"); ax[0].set_xlabel("M_G' from parallax"); ax[0].set_ylabel("M_G' predicted from XP (CV)")
ax[0].legend(); ax[0].set_xlim(-1, 7); ax[0].set_ylim(-1, 7)
ax[1].scatter(m.bp_rp[tr], y[tr], s=2, c="0.6", label="training (parallax M_G')")
ax[1].plot(m.bp_rp[it], t_knn, "r*", ms=14, label="target, XP k-NN")
ax[1].plot(m.bp_rp[it], res["MG_if_dr3_plx"], "bs", ms=8, label="target, DR3 parallax")
ax[1].plot(m.bp_rp[it], 2.04 + 0.139, "m^", ms=8, label="target, GSP-Phot MARCS (998 pc)")
ax[1].invert_yaxis(); ax[1].set_xlabel("BP-RP"); ax[1].set_ylabel("M_G'"); ax[1].legend(fontsize=7)
ax[2].hist(res_k, bins=60, range=(-2, 2), histtype="step", label="kNN CV residual (similar colour)")
ax[2].hist(res_g, bins=60, range=(-2, 2), histtype="step", label="GB CV residual")
ax[2].set_xlabel("M_G' true - pred"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/xp_knn.png"), dpi=110)
