import time, sys
sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
d = load_star("M3V", "B")
t0 = time.time()
tr = detrend(d.t.values, d.raw.values, d.intr.values)
print("detrend s", time.time() - t0)
det2 = d.raw.values / tr
print("max |det2-det| ppt", np.nanmax(np.abs(det2 - d.det.values)) * 1e3, np.nanpercentile(np.abs(det2 - d.det.values), 99) * 1e3)
m = ~d.intr.values
tb, fb, nb_, sb, cb = bin_lc(d.t.values[m], d.det.values[m], d.sector.values[m], d.cad.values[m], 10.0)
print("binned N", len(tb), {s: int((sb == s).sum()) for s in SECTORS})
y, w = prep_weights(fb - 1, sb)
for s in SECTORS: print(s, "rstd binned ppt", 1e3 / np.sqrt(w[sb == s][0]))
B = tb.max() - tb.min()
fr = freq_grid(0.3, 30.0, B)
print("nfreq", len(fr), "baseline", B)
t0 = time.time(); r = run_bls(tb[:100], y[:100], w[:100], fr[:1000]); print("compile", time.time() - t0)
t0 = time.time(); r = run_bls(tb, y, w, fr); print("full bls s", time.time() - t0)
i = np.argmax(r["dchi2"]); print("best P", 1 / fr[i], r["dchi2"][i], r["depth"][i])
