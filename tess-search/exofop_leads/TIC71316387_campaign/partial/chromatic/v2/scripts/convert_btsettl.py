"""Convert SVO BT-Settl (CIFIST2011) VOTable spectra to a compact binned npz grid.
Grid: 3000-12000 A in 2 A bins, 12000-60000 A in 20 A bins (flux density averaged in bin)."""
import re, glob, os, sys
import numpy as np
D = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/data/models/btsettl/'
edges = np.r_[np.arange(3000, 12000, 2.0), np.arange(12000, 60001, 20.0)]
cen = 0.5 * (edges[1:] + edges[:-1])
teffs, loggs, fl = [], [], []
pat = re.compile(rb'<TD>([^<]+)</TD><TD>([^<]+)</TD>')
for fn in sorted(glob.glob(D + 'bt_*.xml')):
    T, g = os.path.basename(fn)[3:-4].split('_')
    s = open(fn, 'rb').read()
    arr = np.array(pat.findall(s), dtype=float)
    w, f = arr[:, 0], arr[:, 1]
    o = np.argsort(w); w, f = w[o], f[o]
    # integrate f over each bin using trapezoid on a cumulative integral
    cum = np.r_[0, np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(w))]
    ci = np.interp(edges, w, cum)
    fb = np.diff(ci) / np.diff(edges)
    teffs.append(float(T)); loggs.append(float(g)); fl.append(fb)
    print(T, g, len(w), w.min(), w.max(), 'median dlam 4000-9000: %.3f' % np.median(np.diff(w[(w > 4000) & (w < 9000)])), flush=True)
np.savez_compressed(D + '../btsettl_grid.npz', wave=cen, edges=edges, teff=np.array(teffs), logg=np.array(loggs), flux=np.array(fl))
print('saved')
