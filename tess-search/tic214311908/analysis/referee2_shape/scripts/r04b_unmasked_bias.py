"""r04b: recovery of a 2-ppt secondary at phase 0.5 when the flux is detrended WITHOUT masking phase 0.5 (as in the
delivered phot/lc_M3V_all.csv flux_norm used by ../shape s10) vs with the mask (r04). Output: out/secondary_unmasked_bias.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
src = open('tic214311908/analysis/referee2_shape/scripts/r04_secondary.py').read()
src = src[:src.index('rows = []')]
exec(src)
# variant of measure() whose detrending masks only the real transits (not the tested window)
def measure_unmasked(ph, raw):
    global fake_mask_off
    T0f = T0REF + ph*PREF
    f = np.full(len(t), np.nan)
    for s in np.unique(sec):
        m = sec == s; use = m & ~real
        tr = flatten(t[use], raw[use], method='biweight', window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
        ok = np.isfinite(tr); f[m] = raw[m]/np.interp(t[m], t[use][ok], tr[ok])
    f = 1 + (f - 1)*935/979
    ep = np.round((t - T0f)/PREF).astype(int); dt = t - (T0f + ep*PREF)
    w = (np.abs(dt) < 0.3) & np.isfinite(f) & ~(np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.3)
    D = dict(t=t[w], f=f[w], sig=sigw[w], sector=sec[w], cad=cad[w], epoch=ep[w], dt=dt[w], wid=sec[w]*1000 + ep[w])
    lc = LC(D); tm = lc.transit(T0f, PREF, 0.04, 34.0, 0.0, 0.125, 0.533) - 1; tm = tm/(-tm.min())
    r = lc.f - 1; x = tm.copy(); wi = lc.wi; ww = lc.w; sw = np.bincount(wi, ww, lc.nw)
    r = r - (np.bincount(wi, ww*r, lc.nw)/sw)[wi]; x = x - (np.bincount(wi, ww*x, lc.nw)/sw)[wi]
    return float(np.sum(ww*x*r)/np.sum(ww*x*x))
D = dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), sector=sec, cad=cad, epoch=np.zeros(len(t), int), dt=((t - T0REF - 0.5*PREF)/PREF + 0.5) % 1*PREF - 0.5*PREF, wid=sec)
L = LC(D); tm = L.transit(T0REF + 0.5*PREF, PREF, 0.04, 34.0, 0.0, 0.125, 0.533) - 1; tm = tm/(-tm.min())
a0 = measure_unmasked(0.5, raw0); a1 = measure_unmasked(0.5, raw0 + 0.002*tm*979/935)
out = dict(real_unmasked_ppt=a0*1e3, recovered_fraction_unmasked=(a1 - a0)/0.002)
print(out); json.dump(out, open(f'{RD}/out/secondary_unmasked_bias.json', 'w'), indent=1)
