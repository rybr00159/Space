"""Main-sequence relations (v6, SMOOTHED).  Source: Pecaut & Mamajek dwarf sequence (EEM_dwarf_UBVIJHK_colors_Teff,
version 2022.04.16), rows 0.09-1.1 Msun.  The raw table is lumpy (R jumps 0.36 -> 0.42 Rsun between M3V and M2.5V,
i.e. 0.37 -> 0.40 Msun), which would imprint artificial edges on a mass posterior, so every quantity is a
least-squares degree-5 polynomial in x = log10(M/Msun):  log R (rms 0.011 dex), log Teff (0.005 dex), M_T (0.064 mag),
M_Ks (0.057 mag), M_G (0.067 mag), Bp-Rp.  The scatter is absorbed by the dR (7%) and dJ (10%) nuisance priors of the
MS-constrained fit.  TESS absolute magnitude from G and Bp-Rp with the TIC-8 relation (Stassun et al. 2019); checked
against M_Ic from the same table: M_T - M_Ic = +0.25 +- 0.05 mag over M0-M5 (smooth, so ratios are unaffected).
TESS-band disc-averaged surface brightness SB_T = 10^(-0.4 M_T) / R^2 (arbitrary units).  Valid 0.09-0.9 Msun."""
import numpy as np
F = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/ext/mamajek.txt'

def _load():
    rows = []
    started = False
    for line in open(F):
        if line.startswith('#SpT'):
            if started: break
            started = True; hdr = line[1:].split(); continue
        if not started or line.startswith('#') or not line.strip(): continue
        v = line.split()
        if len(v) < len(hdr): continue
        rows.append(v[:len(hdr)])
    out = {}
    for j, h in enumerate(hdr):
        col = []
        for r in rows:
            try: col.append(float(r[j]))
            except ValueError: col.append(np.nan)
        out[h] = np.array(col)
    out['SpT'] = np.array([r[0] for r in rows])
    return out

T = _load()
sel = np.isfinite(T['Msun']) & np.isfinite(T['M_G']) & np.isfinite(T['Bp-Rp']) & np.isfinite(T['R_Rsun']) & (T['Msun'] <= 1.1) & (T['Msun'] >= 0.09)
_m = T['Msun'][sel]; o = np.argsort(_m); _m = _m[o]
SPT = T['SpT'][sel][o]
_x = np.log10(_m)
_c = T['Bp-Rp'][sel][o]
_MTraw = T['M_G'][sel][o] - 0.00522555 * _c ** 3 + 0.0891337 * _c ** 2 - 0.633923 * _c + 0.0324473
DEG = 5
FITS = {}; SCAT = {}
for name, y in [('logR', np.log10(T['R_Rsun'][sel][o])), ('logT', np.log10(T['Teff'][sel][o])), ('MT', _MTraw),
                ('MKs', T['M_Ks'][sel][o]), ('MG', T['M_G'][sel][o]), ('bprp', _c)]:
    ok = np.isfinite(y)
    cc = np.polyfit(_x[ok], y[ok], DEG); FITS[name] = cc
    SCAT[name] = float(np.std(y[ok] - np.polyval(cc, _x[ok])))
MMIN, MMAX = 0.09, 0.90

def _p(name, m): return np.polyval(FITS[name], np.log10(np.asarray(m, float)))
def teff(m): return 10 ** _p('logT', m)
def radius(m): return 10 ** _p('logR', m)
def MT(m): return _p('MT', m)
def MKs(m): return _p('MKs', m)
def MGmag(m): return _p('MG', m)
def bprp(m): return _p('bprp', m)
def SB_T(m): return 10 ** (-0.4 * MT(m)) / radius(m) ** 2
def spt(m):
    return SPT[np.argmin(np.abs(_x - np.log10(m)))]
def density(m):   # solar units
    return np.asarray(m, float) / radius(m) ** 3
_mg = np.exp(np.linspace(np.log(MMIN), np.log(MMAX), 3000)); _rho = density(_mg)
assert np.all(np.diff(_rho) < 0), 'density not monotonic in mass'
def mass_from_density(rho):
    return np.interp(rho, _rho[::-1], _mg[::-1], left=np.nan, right=np.nan)
_sb = SB_T(_mg)
assert np.all(np.diff(_sb) > 0), 'surface brightness not monotonic in mass'
def m2_from_J(J, m1):
    """secondary mass whose TESS surface brightness ratio to a primary of mass m1 equals J"""
    return np.interp(np.asarray(J) * SB_T(m1), _sb, _mg, left=np.nan, right=np.nan)

if __name__ == '__main__':
    print('rms about smooth fits:', {k: round(v, 4) for k, v in SCAT.items()})
    for m in (0.12, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9):
        print(m, spt(m), round(float(teff(m))), round(float(radius(m)), 3), round(float(MT(m)), 2), round(float(MKs(m)), 2), round(float(density(m)), 2))
