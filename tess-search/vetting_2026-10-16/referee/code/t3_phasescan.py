"""TIC 214311908 aperture (M dwarf next to the Mira): S/N of the box depth at the
5.69149-d ephemeris vs at 200 trial phases (empirical false-alarm check incl. red noise)."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from wotan import flatten
from reflib import *
P, T0, T14 = 5.69149, 3840.913, 1.6 / 24
OUT = {}
allt, allf = [], []
for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(214311906, sec)
    g = gaia_table(214311906, btjd_to_jyear(np.median(cut['time'])))
    a = g[g.source_id == 5925373272850784768].iloc[0]
    t, f, e, ex = aperture_lc(214311906, sec, r=1.0, center=(a['ra_ep'], a['dec_ep']))
    fl = flatten(t, f, method='biweight', window_length=0.5, break_tolerance=0.5)
    m = np.isfinite(fl)
    allt.append(t[m]); allf.append(fl[m])
def snr_at(t, f, phase0):
    ph = ((t - T0) / P - phase0 + 0.5) % 1 - 0.5
    inn = np.abs(ph * P) < 0.4 * T14
    near = (np.abs(ph * P) > 0.6 * T14) & (np.abs(ph * P) < 0.6 * T14 + 0.15)
    if inn.sum() < 5:
        return np.nan
    # depth relative to local surroundings, error from per-event scatter of in-transit means
    d = np.mean(f[near]) - np.mean(f[inn])
    s = robust_std(f) / np.sqrt(inn.sum())
    return d / s
for name, (t, f) in [(f'S{s}', (tt, ff)) for s, tt, ff in zip((12, 39, 66, 93, 103, 104), allt, allf)] + [('all', (np.concatenate(allt), np.concatenate(allf)))]:
    s0 = snr_at(t, f, 0.0)
    trial = np.array([snr_at(t, f, p) for p in np.linspace(0.05, 0.95, 200)])
    trial = trial[np.isfinite(trial)]
    fap = np.mean(trial >= s0)
    OUT[name] = dict(snr=s0, trial_max=float(trial.max()), trial_std=float(trial.std()), frac_trials_ge=float(fap))
    print(f'{name:5s}: S/N at ephemeris {s0:5.1f}; trial phases: sd {trial.std():.2f}, max {trial.max():.1f}, fraction >= observed {fap:.3f}')
save_json(OUT, f'{BASE}/results/t3_phasescan.json')
