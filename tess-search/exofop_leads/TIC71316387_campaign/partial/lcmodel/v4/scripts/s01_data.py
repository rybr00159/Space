"""Build the eclipse-window data sets used by every fit (re-derived from the light curves, not reused):
  default : TESS-SPOC PDCSAP S1/S28/S68 (crowding-corrected by SPOC) + own TESScut aperture photometry S95/S105/S106
            (target-specific barycentric times, PRF/Gaia contamination removed).
  qlp     : SPOC S1/S28/S68 + QLP SAP S95/S105 (time-corrected) + TESScut S106 (no QLP for S106).
  tc      : own TESScut photometry for all six sectors.
Each eclipse window is +-0.20 d around the predicted conjunction; per-event multiplicative quadratic baseline is profiled
in the likelihood (fitlib.profile_chi2), so no prior detrending touches the eclipses.  Per-sector white noise from OOE
residuals, inflated by the Pont et al. time-averaging beta (30-90 min)."""
import numpy as np, sys, pickle, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import lclib as L, fitlib as F
OUT = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results'
sets = {
    'default': lambda: L.default_set(),
    'qlp': lambda: [L.load_spoc(1), L.load_spoc(28), L.load_spoc(68), L.load_qlp(95), L.load_qlp(105), L.load_tc(106)],
    'tc': lambda: [L.load_tc(s) for s in (1, 28, 68, 95, 105, 106)],
}
summ = {}
for name, fn in sets.items():
    print(name)
    d = F.prepare(fn(), hw=0.20, verbose=True)
    pickle.dump(d, open(f'{OUT}/data_{name}.pkl', 'wb'))
    ev = d['evtab']
    summ[name] = dict(npts=int(d['t'].size), nev=len(ev),
                      per_sector={str(s): dict(prod=v['prod'], n_prim=sum(1 for e in ev if e['sector'] == s and e['ph'] == 0.0),
                                               n_sec=sum(1 for e in ev if e['sector'] == s and e['ph'] == 0.5),
                                               sigma_ppm=v['sigma'] * 1e6, beta=v['beta']) for s, v in d['info'].items()})
    print('  points', d['t'].size, 'events', len(ev))
json.dump(summ, open(f'{OUT}/data_summary.json', 'w'), indent=1)
