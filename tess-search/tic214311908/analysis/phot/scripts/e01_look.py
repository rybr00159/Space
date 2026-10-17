import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
for ax, s in zip(axs.ravel(), SECTORS):
    C = Cutout(s); g = C.gaia(17)
    good = C.fin & (C.q == 0)
    med = np.median(C.flux[good], 0)
    ax.imshow(np.log10(np.clip(med - np.percentile(med, 5) + 10, 1, None)), origin="lower", cmap="gray_r")
    m = (g.x > -0.5) & (g.x < 20.5) & (g.y > -0.5) & (g.y < 20.5)
    ax.scatter(g.x[m], g.y[m], s=np.clip(5 * (17 - g["T"][m]) ** 1.5, 2, 200), facecolors="none", edgecolors="c", lw=0.7)
    for k, c in (("Mira", "r"), ("M3V", "y"), ("pairA", "m"), ("ref", "g")):
        r = g[g.source_id == IDS[k]].iloc[0]; ax.plot(r.x, r.y, "+", color=c, ms=10, mew=2)
    ax.set_title(f"S{s} cam{C.camera} ccd{C.ccd} t={C.exptime*86400:.0f}s")
    print(s, "Mira xy", g[g.source_id == IDS["Mira"]][["x", "y"]].values, "M3V", g[g.source_id == IDS["M3V"]][["x", "y"]].values,
          "J", np.round(C.jacobian(10, 10), 2).tolist())
fig.tight_layout(); fig.savefig(f"{ROOT}/explore/e01_images.png", dpi=70)
