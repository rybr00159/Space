"""Combine the per-sector free-position fits of the in-transit difference images (p14) into one dimming centroid
relative to the M3V, and compare it with the Gaia positions of nearby sources (statistical errors only; add ~1-2"
for PRF/WCS systematics: field-star WCS residuals 0.01-0.04 px, M-dwarf lateral colour 0.02-0.1 px)."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
V = json.load(open(f"{ROOT}/out/validate_diffimg.json"))
E, N, eE, eN = [], [], [], []
for s in SECTORS:
    f = V[str(s)]["fits"]["free"]; E.append(f["offset_from_M3V_EN"][0]); N.append(f["offset_from_M3V_EN"][1]); eE.append(f["offset_err"][0]); eN.append(f["offset_err"][1])
E, N, eE, eN = map(np.array, (E, N, eE, eN))
wE, wN = 1 / eE ** 2, 1 / eN ** 2
mE, mN = np.sum(wE * E) / wE.sum(), np.sum(wN * N) / wN.sum(); sE, sN = 1 / np.sqrt(wE.sum()), 1 / np.sqrt(wN.sum())
chiE, chiN = np.sum(wE * (E - mE) ** 2), np.sum(wN * (N - mN) ** 2)
g = pd.read_csv(GAIA); m3 = g[g.source_id == IDS["M3V"]].iloc[0]
g["oE"] = (g.ra - m3.ra) * np.cos(np.radians(m3.dec)) * 3600; g["oN"] = (g.dec - m3.dec) * 3600   # epoch 2016 relative offsets
# offsets at the precision-weighted mean epoch of the sectors (the M3V moves (-114, +86) mas/yr)
ep = {12: 2019.4, 39: 2021.5, 66: 2023.5, 93: 2025.5, 103: 2026.2, 104: 2026.3}
wt = 1 / (eE ** 2 + eN ** 2); epoch = float(np.sum(wt * np.array([ep[s] for s in SECTORS])) / wt.sum()); dt = epoch - 2016.0
g["oE"] = g.oE + (np.nan_to_num(g.pmra) - m3.pmra) * dt / 1000; g["oN"] = g.oN + (np.nan_to_num(g.pmdec) - m3.pmdec) * dt / 1000
print(f"offsets evaluated at epoch {epoch:.1f}")
near = g[(np.hypot(g.oE, g.oN) < 25) & (g["T"] < 19)].copy()
out = dict(epoch=epoch, centroid_E=mE, centroid_E_err=sE, centroid_N=mN, centroid_N_err=sN, chi2_E=chiE, chi2_N=chiN, dof=len(E) - 1, sources=[])
for sys_ in (0.0, 1.5):
    for r in near.itertuples():
        dE_, dN_ = mE - r.oE, mN - r.oN
        c2 = dE_ ** 2 / (sE ** 2 + sys_ ** 2) + dN_ ** 2 / (sN ** 2 + sys_ ** 2)
        out["sources"].append(dict(source_id=int(r.source_id), G=float(r.G), T=float(r.T), oE=float(r.oE), oN=float(r.oN), sys_arcsec=sys_, chi2_2dof=float(c2)))
S = pd.DataFrame(out["sources"])
print(f"dimming centroid relative to the M3V (Gaia 2016 frame): E {mE:+.2f} +- {sE:.2f}\", N {mN:+.2f} +- {sN:.2f}\" (sector scatter chi2 E {chiE:.1f}, N {chiN:.1f} for {len(E)-1} dof)")
print(S.sort_values(["sys_arcsec", "chi2_2dof"]).round(2).to_string())
save_json(out, f"{ROOT}/out/centroid_combined.json")
