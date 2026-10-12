"""Convert MIST v1.2 (vvcrit 0.4) basic isochrones for [Fe/H] = -1.0 ... +0.5 into a compact npz table
(log_age 8.5-10.15, phases 0 = MS and 2 = SGB/RGB, 0.5 < M_init < 3.0).
Output: cache/mist_compact.npz with columns feh, logage, eep, mini, mass, logL, logT, logR, logg, phase."""
import glob, os
import numpy as np
D = "/home/user/work/tic87/stellar/cache/mist/MIST_v1.2_vvcrit0.4_basic_isos"
rows = []
for fn in sorted(glob.glob(os.path.join(D, "*.iso"))):
    tag = os.path.basename(fn).split("_")[3]
    feh = (-1 if tag[0] == "m" else 1) * float(tag[1:])
    a = np.loadtxt(fn, comments="#")
    # columns: 0 EEP 1 log_age 2 initial_mass 3 star_mass 7 log_L 10 log_Teff 11 log_R 12 log_g 24 phase
    m = (a[:, 1] >= 8.5) & (a[:, 1] <= 10.15) & np.isin(a[:, 24], [0, 2]) & (a[:, 2] > 0.5) & (a[:, 2] < 3.0)
    b = a[m]
    rows.append(np.column_stack([np.full(len(b), feh), b[:, 1], b[:, 0], b[:, 2], b[:, 3], b[:, 7], b[:, 10], b[:, 11], b[:, 12], b[:, 24]]))
    print(fn, feh, len(b))
t = np.vstack(rows)
np.savez_compressed("/home/user/work/tic87/stellar/cache/mist_compact.npz", tab=t,
                    cols=np.array(["feh", "logage", "eep", "mini", "mass", "logL", "logT", "logR", "logg", "phase"]))
