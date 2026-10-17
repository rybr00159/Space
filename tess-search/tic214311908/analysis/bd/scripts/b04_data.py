"""b04: light-curve data sets for the shape fits and the injection calibration.
REF (independent local-window PRF photometry, ../referee_phot): real windows = lc_M3V_referee_windows.csv (strict
  selection, 21 transits); FAKE windows = the referee's own fake-epoch windows (data/win_s*.npz 'fake_*': identical
  processing, per-pixel baselines fitted out-of-'transit' |dt| > 0.07 d, 62-104 per sector), flux = 1 + aM/979,
  per-window error = robust OOT scatter /979 (as r23_export).
INV-clean (model-B deblended photometry, ../phot, with the Gaia DR3 RRab leak removed by ../referee2_search r01): the
  raw model-B flux (flux_raw_norm, rescaled to the 979 e-/s scale) is re-detrended here with wotan biweight 0.75 d
  (as phot/p06), masking +-0.1 d around the real AND (for fake sets) the fake epochs, then the referee2_search
  3-harmonic RRab model (per sector, 'model' in m3v_clean_r2.npz) is subtracted. Real windows use the same code with only
  the real epochs masked (checked against the delivered cleaned flux). Fake ephemerides: 48 phase offsets in
  0.12-0.38 and 0.62-0.88. Flares (flare_flag) and quality_ok == 0 removed. Windows +-0.30 d.
Outputs: data/ref_real.npz, data/ref_fake.npz, data/inv_real.npz, data/inv_fake.npz, out/b04_data.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
from wotan import flatten
info = {}
# ---------------- REF
r = pd.read_csv(f'{B}/referee_phot/lc_M3V_referee_windows.csv')
r = r[np.isfinite(r.flux_norm) & np.isfinite(r.flux_err)]
np.savez(f'{BD}/data/ref_real.npz', t=r.time_btjd_tdb.values, f=r.flux_norm.values, e=r.flux_err.values, sec=r.sector.values,
         cad=r.cadence_s.values.astype(float), win=(r.sector.values*1000 + r.epoch.values), tc=T0 + r.epoch.values*P, oot=r.oot.values.astype(bool))
info['ref_real'] = dict(n=len(r), n_windows=int(r.groupby(['sector', 'epoch']).ngroups))
FT, FF, FE, FS, FC, FW, FTC, FO = [], [], [], [], [], [], [], []
for s in (12, 39, 66, 93, 103, 104):
    z = np.load(f'{B}/referee_phot/data/win_s{s}.npz', allow_pickle=True)
    t = z['fake_t']; aM = z['fake_aM']; w = z['fake_win']; oot = z['fake_oot']; tc = z['fake_tc']; cad = float(z['exptime'])*DAY
    for i in range(len(tc)):
        m = w == i
        if m.sum() < 10: continue
        sg = rstd(aM[m & oot])
        FT.append(t[m]); FF.append(1 + aM[m]/FM3V); FE.append(np.full(m.sum(), sg/FM3V)); FS.append(np.full(m.sum(), s))
        FC.append(np.full(m.sum(), round(cad))); FW.append(np.full(m.sum(), s*1000 + i)); FTC.append(np.full(m.sum(), tc[i])); FO.append(oot[m])
cat = lambda L: np.concatenate(L)
np.savez(f'{BD}/data/ref_fake.npz', t=cat(FT), f=cat(FF), e=cat(FE), sec=cat(FS), cad=cat(FC).astype(float), win=cat(FW), tc=cat(FTC), oot=cat(FO))
info['ref_fake'] = dict(n=int(len(cat(FT))), n_windows=int(len(np.unique(cat(FW)))))
# ---------------- INV-clean
d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
d = d[(d.quality_ok == 1) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].copy()
z = np.load(f'{B}/referee2_search/data/m3v_clean_r2.npz')
mod = pd.Series(z['model'], index=np.round(z['t'], 7))
clean_deliv = pd.Series(z['clean'], index=np.round(z['t'], 7))
d['tk'] = np.round(d.time_btjd_tdb.values, 7)
d = d[d.tk.isin(mod.index)].copy()
d['rr'] = mod.loc[d.tk].values; d['clean_deliv'] = clean_deliv.loc[d.tk].values
SC = 935.97/979.0
t = d.time_btjd_tdb.values; raw = 1 + SC*(d.flux_raw_norm.values - 1); e = d.flux_err.values*SC
sec = d.sector.values.astype(int); cad = d.cadence_s.values.astype(float); flare = d.flare_flag.values.astype(bool)
def phdist(tt, t0):
    return (((tt - t0)/P + 0.5) % 1 - 0.5)*P
real = np.abs(phdist(t, T0)) < 0.1
def detrend(mask):
    f = np.full(len(t), np.nan)
    for s in np.unique(sec):
        m = sec == s; use = m & ~mask
        tr = flatten(t[use], raw[use], method='biweight', window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
        ok = np.isfinite(tr)
        f[m] = raw[m]/np.interp(t[m], t[use][ok], tr[ok])
    return f - d.rr.values          # RRab leak removed (additive, fractional units)
f_real = detrend(real)
dtr = phdist(t, T0); wr = (np.abs(dtr) < 0.30) & ~flare & np.isfinite(f_real)
ep = np.round((t - T0)/P).astype(int)
np.savez(f'{BD}/data/inv_real.npz', t=t[wr], f=f_real[wr], e=e[wr], sec=sec[wr], cad=cad[wr], win=sec[wr]*1000 + ep[wr] % 1000,
         tc=(T0 + ep*P)[wr], oot=(np.abs(dtr[wr]) > 0.07))
dd = f_real[wr] - d.clean_deliv.values[wr]
info['inv_real'] = dict(n=int(wr.sum()), n_windows=int(len(np.unique(sec[wr]*1000 + ep[wr] % 1000))),
                        check_vs_delivered_clean=dict(median_ppt=float(np.nanmedian(dd)*1e3), rms_ppt=float(np.nanstd(dd)*1e3)))
print(info, flush=True)
phs = np.r_[np.linspace(0.12, 0.38, 24), np.linspace(0.62, 0.88, 24)]
FT, FF, FE, FS, FC, FW, FTC, FO, FSET = [], [], [], [], [], [], [], [], []
for j, ph in enumerate(phs):
    T0f = T0 + ph*P
    fake = np.abs(phdist(t, T0f)) < 0.1
    ff = detrend(real | fake)
    dtf = phdist(t, T0f); w = (np.abs(dtf) < 0.30) & ~flare & np.isfinite(ff) & (np.abs(dtr) > 0.30)
    epf = np.round((t - T0f)/P).astype(int)
    FT.append(t[w]); FF.append(ff[w]); FE.append(e[w]); FS.append(sec[w]); FC.append(cad[w])
    FW.append(j*100000 + sec[w]*1000 + epf[w] % 1000); FTC.append((T0f + epf*P)[w]); FO.append(np.abs(dtf[w]) > 0.07); FSET.append(np.full(w.sum(), j))
    print('fake phase', round(ph, 3), 'n', int(w.sum()), flush=True)
np.savez(f'{BD}/data/inv_fake.npz', t=cat(FT), f=cat(FF), e=cat(FE), sec=cat(FS), cad=cat(FC), win=cat(FW), tc=cat(FTC), oot=cat(FO), fset=cat(FSET), phases=phs)
info['inv_fake'] = dict(n_sets=len(phs), n=int(len(cat(FT))))
jdump(dict(method=__doc__, info=info), f'{BD}/out/b04_data.json')
