"""eem: Mamajek 'Modern Mean Dwarf Stellar Color and Effective Temperature Sequence' (copy downloaded by the
referee_star agent, ../referee_star/data/EEM_dwarf.txt) as interpolation functions of mass for main-sequence
K/M dwarfs: Teff, R, M_G, Bp-Rp, G-Rp, M_T (TESS: G - T from Stassun+19 for Bp-Rp <= 3, linear extrapolation of
its slope (0.24 mag/mag) beyond), M_V, V-Ic, g-r / i-z (SDSS-like) and BCv. Brown dwarfs: L_T = 0, radius given."""
import numpy as np, pandas as pd
EEMF = 'tic214311908/analysis/referee_star/data/EEM_dwarf.txt'
def load():
    lines = open(EEMF).read().splitlines()
    hdr = lines[22].lstrip('#').split()
    rows = []
    for l in lines[23:140]:
        if l.startswith('#'): break
        rows.append(l.split()[:len(hdr)])
    d = pd.DataFrame(rows, columns=hdr[:len(rows[0])])
    for c in d.columns[1:]:
        d[c] = pd.to_numeric(d[c].astype(str).str.rstrip(':'), errors='coerce')
    return d
EEM = load()
def g_minus_t(bprp):
    x = np.minimum(bprp, 3.0)
    gt = 0.00522555*x**3 - 0.0891337*x**2 + 0.633923*x - 0.0324473
    return gt + 0.24*np.clip(bprp - 3.0, 0, None)
MS = EEM[(EEM.Msun.notna()) & (EEM.Msun >= 0.078) & (EEM.Msun <= 1.6) & EEM.M_G.notna() & EEM['Bp-Rp'].notna()].copy()
MS = MS.sort_values('Msun').drop_duplicates('Msun')
MS['M_T'] = MS.M_G - g_minus_t(MS['Bp-Rp'])
def ms(m, col):
    return np.interp(m, MS.Msun.values, MS[col].values)
def ms_T(m): return ms(m, 'M_T')
def ms_R(m): return ms(m, 'R_Rsun')
def ms_Teff(m): return ms(m, 'Teff')
# ---- Mann et al. (2015, 2019) relations for M dwarfs (metallicity-insensitive Ks-band luminosity budget)
MANN19 = (-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4)   # log10 M(MKs), n = 5, no [Fe/H] term
def mann_mass(MKs):
    x = MKs - 7.5
    return 10**sum(c*x**i for i, c in enumerate(MANN19))
def mann_radius(MKs):
    return 1.9515 - 0.3520*MKs + 0.01680*MKs**2
_mk = np.linspace(4.0, 11.0, 2000); _mm = mann_mass(_mk)
def mann_MKs(M):
    """inverse of Mann+19 (valid 0.075-0.70 Msun)"""
    return float(np.interp(M, _mm[::-1], _mk[::-1]))
MS['T_minus_Ks'] = MS.M_T - MS.M_Ks
def t_minus_ks(m): return ms(m, 'T_minus_Ks')
