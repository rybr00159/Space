"""rlib: shared helpers for the referee-2 re-derivation of the FPP of the 5.69-d signal on TIC 214311908.
Own code. Reads (read-only) phase-1/phase-2 data products; writes only into referee2_fpp/."""
import numpy as np, json, os
B = 'tic214311908/analysis'
R2 = f'{B}/referee2_fpp'
G = 6.674e-8; RSUN = 6.957e10; MSUN = 1.989e33; AU = 1.496e13; DAY = 86400.0; REARTH = 6.371e8
P = 5.6915087
T0 = 3840.9141          # phase-1 ephemeris (phot/out/timing.json); shape free fits 3840.9150 / 3840.9160
DIST = 64.31
MKS_OBS = 6.4406


def jdump(obj, path):
    json.dump(obj, open(path if path.startswith('/') else f'{R2}/{path}', 'w'), indent=1,
              default=lambda o: o.tolist() if hasattr(o, 'tolist') else (float(o) if isinstance(o, np.floating) else str(o)))


def jload(path):
    return json.load(open(path if path.startswith('/') else f'{B}/{path}'))


# ---------------- Mann et al. 2019 (mass-MKs, 5th order, zp 7.5, no [Fe/H]) and Mann et al. 2015 (R-MKs)
_M19 = (-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4)
def mass_of_mks(mk):
    x = np.asarray(mk, float) - 7.5
    return 10**sum(c*x**i for i, c in enumerate(_M19))
def radius_of_mks(mk):
    mk = np.asarray(mk, float)
    return 1.9515 - 0.3520*mk + 0.01680*mk**2
_mkg = np.linspace(4.0, 11.5, 6000); _mmg = mass_of_mks(_mkg)
def mks_of_mass(m):
    return np.interp(np.asarray(m, float), _mmg[::-1], _mkg[::-1])


# ---------------- Mamajek dwarf table (data file downloaded in phase 1 by referee_star), own parser
def eem_table():
    import pandas as pd
    L = open(f'{B}/referee_star/data/EEM_dwarf.txt').read().splitlines()
    i0 = [i for i, l in enumerate(L) if l.startswith('#SpT')][0]
    hdr = L[i0].lstrip('#').split()
    rows = []
    for l in L[i0+1:]:
        if l.startswith('#') or not l.strip(): break
        rows.append(l.split()[:len(hdr)])
    d = pd.DataFrame(rows, columns=hdr[:len(rows[0])])
    for c in d.columns[1:]:
        d[c] = pd.to_numeric(d[c].astype(str).str.rstrip(':'), errors='coerce')
    d = d[(d.Msun > 0.07) & (d.Msun < 0.75) & d.M_G.notna() & d['Bp-Rp'].notna() & d.M_Ks.notna()].copy()
    d = d.sort_values('Msun').drop_duplicates('Msun')
    # Gaia BP, RP absolute mags; TESS via Stassun+19 G-T(Bp-Rp) (valid to Bp-Rp ~3; linear slope beyond)
    d['M_RP'] = d.M_G - d['G-Rp']
    d['M_BP'] = d.M_RP + d['Bp-Rp']
    x = np.minimum(d['Bp-Rp'], 3.0)
    gt = 0.00522555*x**3 - 0.0891337*x**2 + 0.633923*x - 0.0324473 + 0.24*np.clip(d['Bp-Rp'] - 3.0, 0, None)
    d['M_T'] = d.M_G - gt
    d['M_I'] = d.Mv - d['V-Ic']
    return d
_E = None
def eem(m, col):
    """interpolate a Mamajek column at mass m, in the space of MKs (smoother than mass)."""
    global _E
    if _E is None: _E = eem_table()
    return np.interp(np.asarray(m, float), _E.Msun.values, _E[col].values)
def colour_offsets(m, col):
    """col - Ks for mass m from the Mamajek table, applied to Mann MKs (keeps the Ks luminosity budget of Mann+19)"""
    return eem(m, col) - eem(m, 'M_Ks')
