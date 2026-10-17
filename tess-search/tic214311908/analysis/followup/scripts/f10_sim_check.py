"""f10: Monte Carlo check of the per-event depth-error model used in f06.
Simulated ground light curves: 2-min bins over the full 4.9-h window, trapezoid transit (T14 1.70 h, T23/T14 0.73, depth
15.2 ppt), white noise sigma_2min, plus an Ornstein-Uhlenbeck (AR(1)) red component with tau = 30 min and stationary rms
sigma_ou, plus an airmass-like quadratic trend. Fit: depth x fixed template + quadratic baseline (linear least squares),
mid-time either known or offset by a random 2-min (1-sigma) ephemeris error that the fit does not know.
Compares the scatter of the recovered depth with the analytic f06 white term sigma_2min sqrt(1/N_in + 1/N_out), and
gives the per-event depth error as a function of sigma_ou (the 'sigma_red' of f06 is the red term of this table).
Output: out/sim_depth_error.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd
rng = np.random.default_rng(42)
dt = 2/1440; hw = (T14_H/2 + 1.5 + 0.1)/24
t = np.arange(-hw, hw + 1e-9, dt)
def trap(t, t0=0.0):
    T14 = T14_H/24; T23 = 0.73*T14; x = np.abs(t - t0)
    return np.where(x <= T23/2, 1.0, np.where(x >= T14/2, 0.0, (T14/2 - x)/((T14 - T23)/2)))
X = lambda t0: np.vstack([trap(t, t0), np.ones_like(t), t, t*t]).T
rows = []
for s2 in (0.8, 1.5, 2.6, 4.4):
    for sou in (0.0, 0.5, 1.0):
        for terr in (0.0, 2.0):
            dep = []
            for k in range(1500):
                phi = np.exp(-dt/(30/1440)); e = np.zeros_like(t); e[0] = rng.normal(0, sou*1e-3)
                for i in range(1, len(t)): e[i] = phi*e[i - 1] + rng.normal(0, sou*1e-3*np.sqrt(1 - phi**2))
                y = -DEPTH*trap(t, rng.normal(0, terr/1440)) + rng.normal(0, s2*1e-3, len(t)) + e + 1e-3*(0.5*t*24 + 0.3*(t*24)**2)
                b, *_ = np.linalg.lstsq(X(0.0), y, rcond=None)
                dep.append(-b[0])
            dep = np.array(dep)
            nin = np.sum(trap(t) > 0.5); nout = np.sum(trap(t) == 0)
            white = s2*np.sqrt(1/nin + 1/nout)
            rows.append(dict(sigma_2min_ppt=s2, sigma_ou_ppt=sou, timing_err_min=terr, depth_bias_ppt=round((dep.mean() - DEPTH)*1e3, 3),
                             depth_sd_ppt=round(dep.std()*1e3, 3), analytic_white_ppt=round(white, 3),
                             ratio_sd_to_white=round(dep.std()*1e3/white, 2),
                             red_term_ppt=round(float(np.sqrt(max((dep.std()*1e3)**2 - (white*dep.std()*1e3/white if sou == 0 else white)**2, 0))), 3)))
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f'{OUT}/sim_depth_error.csv', index=False)
