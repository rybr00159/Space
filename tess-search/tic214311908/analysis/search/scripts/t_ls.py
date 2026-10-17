import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from astropy.timeseries import LombScargle
for star, model in [("M3V","B"),("M3V","A"),("pair","A"),("pair","Bp"),("Mira","B"),("Mira","A")]:
    d = load_star(star, model if model!="Bp" else "B")
    m = ~d.intr.values & np.isfinite(d.det.values)
    tb, fb, nb_, sb, cb = bin_lc(d.t.values[m], d.det.values[m], d.sector.values[m], d.cad.values[m], 30.0)
    y = fb.copy()
    for s in SECTORS:
        k = sb==s; y[k] = (y[k]-np.median(y[k]))
    fr = np.linspace(1/5, 1/0.25, 200000)
    L = LombScargle(tb, y); pw = L.power(fr)
    j = np.argmax(pw); 
    # amplitude at 0.6485
    f0 = 1/0.64853
    jj = np.argmin(np.abs(fr - f0)); win = (np.abs(fr-f0) < 0.002)
    print(star, model, "best P %.5f pow %.4f" % (1/fr[j], pw[j]), "near 0.6485: P %.5f pow %.4f" % (1/fr[win][np.argmax(pw[win])], pw[win].max()))
    # per sector
    out=[]
    for s in SECTORS:
        k = sb==s
        L2 = LombScargle(tb[k], y[k]); f2 = np.linspace(1/5, 1/0.25, 20000); p2 = L2.power(f2)
        win2 = np.abs(f2 - f0) < 0.01
        amp = np.sqrt(4*p2[win2].max()*np.var(y[k]))
        out.append("S%d best %.4f (%.3f) @0.6485: P %.4f pow %.3f A~%.1fppt" % (s, 1/f2[np.argmax(p2)], p2.max(), 1/f2[win2][np.argmax(p2[win2])], p2[win2].max(), amp*1e3))
    print("   " + " | ".join(out))
