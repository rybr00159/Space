"""Sanity check of the OOE fold: per sector, fit trend + 4 harmonics simultaneously (full LSQ), subtract only the fitted trend,
fold at P_EB and bin (0.025 in phase).  Also with 8 harmonics (P/3, P/4 included) to look for structure."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/scripts')
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import importlib.util
spec = importlib.util.spec_from_file_location('ooe', '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/scripts/s05_ooe.py')
import lclib as L
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
# re-implement minimal pieces (avoid re-running s05)
import glob, pandas as pd
from astropy.io import fits
from numpy.polynomial import legendre as LG
P, T0 = L.P0, L.T00
def spoc(sec):
    f = glob.glob(f'{L.RD}/mastDownload/HLSP/*0071316387-s{sec:04d}*/*lc.fits')[0]
    d = fits.open(f)[1].data; t = np.array(d['TIME'], float); fl = np.array(d['PDCSAP_FLUX'], float); q = np.array(d['QUALITY'])
    m = np.isfinite(t) & np.isfinite(fl) & (q == 0) & (fl > 0); return dict(t=t[m], f=fl[m] / np.median(fl[m]))
def prep(lc):
    t, f = lc['t'], lc['f']; ph = ((t - T0) / P) % 1
    dd = np.minimum(np.minimum(ph, 1 - ph), np.abs(ph - 0.5)) * P; m = dd > 0.06; t, f = t[m], f[m]
    rm = pd.Series(f, index=pd.to_datetime(t, unit='D')).rolling('0.5D', center=True, min_periods=5).median().values
    r = f - rm; s = 1.4826 * np.nanmedian(np.abs(r)); k = ~(np.abs(r) > 5 * s); t, f = t[k], f[k]
    br = np.r_[0, np.where(np.diff(t) > 0.5)[0] + 1, t.size]; seg = np.zeros(t.size, int)
    for i in range(len(br) - 1): seg[br[i]:br[i + 1]] = i
    return t, f, seg
def trend(t, seg, dper=2.5):
    cols = []
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < 20: continue
        ts = t[m]; x = 2 * (ts - ts.min()) / (ts.max() - ts.min()) - 1; deg = int(max(2, np.ceil((ts.max() - ts.min()) / dper)))
        for j in range(deg + 1):
            c = np.zeros(t.size); c[m] = LG.legval(x, np.eye(deg + 1)[j]); cols.append(c)
    return np.array(cols).T
fig, axs = plt.subplots(3, 2, figsize=(11, 9)); axs = axs.ravel()
for ax, sec in zip(axs, (1, 28, 68, 95, 105, 106)):
    lc = spoc(sec) if sec in (1, 28, 68) else L.load_tc(sec)
    t, f, seg = prep(lc); A = trend(t, seg); ph = 2 * np.pi * (t - T0) / P
    H = np.vstack([np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)]).T
    X = np.hstack([A, H]); c, *_ = np.linalg.lstsq(X, f - 1, rcond=None)
    y = f - 1 - A @ c[:A.shape[1]]
    phase = ((t - T0) / P) % 1; b = np.linspace(0, 1, 41); idx = np.digitize(phase, b)
    bx = [(b[j - 1] + b[j]) / 2 for j in range(1, 41) if (idx == j).sum() > 10]; by = [np.mean(y[idx == j]) * 1e6 for j in range(1, 41) if (idx == j).sum() > 10]
    ax.plot(bx, by, 'ko', ms=3); xx = np.linspace(0, 1, 200); ph2 = 2 * np.pi * xx
    ax.plot(xx, (np.vstack([np.cos(ph2), np.sin(ph2), np.cos(2 * ph2), np.sin(2 * ph2)]).T @ c[A.shape[1]:]) * 1e6, 'r-')
    ax.set_title(f'S{sec}: harmonics (ppm) a1={c[-4]*1e6:.0f} b1={c[-3]*1e6:.0f} a2={c[-2]*1e6:.0f} b2={c[-1]*1e6:.0f}', fontsize=9)
    ax.set_xlabel('phase'); ax.set_ylabel('ppm')
plt.tight_layout(); plt.savefig('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/plots/ooe_fold_sectors.png', dpi=100)
print('ok')
